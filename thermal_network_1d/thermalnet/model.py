"""モデル定義層 --- YAML/JSON に書いた「物理量」から熱回路を組む。

lessons/ では抵抗値を Python で直接計算していたが、実務では
「熱伝導率・寸法・発熱量・温度指定・流体への排熱」を編集して回したい。
このモジュールはその編集レイヤを提供する。

    python -m thermalnet.model models/cpu_heatsink.yaml
    python -m thermalnet.model models/cpu_heatsink.yaml --set nodes.die.heat=120
    python -m thermalnet.model models/cpu_heatsink.yaml --sweep nodes.die.heat=50:150:11
    python -m thermalnet.model models/cpu_heatsink.yaml --transient --t-end 300

書式は docs/model_format.md を参照。要点だけ:

    nodes:      発熱 heat [W] / 温度指定 fixed [degC] / 熱容量 capacitance [J/K]
    elements:   ノード間の熱抵抗。type ごとに必要な物理量を書く
                conduction / cylinder / interface / contact / convection /
                radiation / fin_array / spreading / flow / resistance

温度に依存する要素(自然対流の相関式・輻射)は自動的に非線形の枝になり、
solve_steady が反復して解く。
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from typing import Any, Callable, Dict, List

import numpy as np

from . import elements as el
from . import fluids, materials
from .network import ThermalNetwork


# =====================================================================
# ヘルパ
# =====================================================================
def _get(spec: dict, *names, default=None, required=False, where=""):
    """spec から最初に見つかったキーの値を返す(別名を許すため)。"""
    for name in names:
        if name in spec and spec[name] is not None:
            return spec[name]
    if required:
        raise KeyError(f"{where}: '{names[0]}' が必要です (指定: {sorted(spec)})")
    return default


def _area(spec: dict, where: str) -> float:
    """面積を得る。area [m^2] か、size(正方形の一辺)か、width x length で指定。"""
    if "area" in spec:
        return float(spec["area"])
    if "size" in spec:
        return float(spec["size"]) ** 2
    if "width" in spec and "length" in spec:
        return float(spec["width"]) * float(spec["length"])
    if "diameter" in spec:
        return math.pi * (float(spec["diameter"]) / 2.0) ** 2
    raise KeyError(f"{where}: 面積が要ります (area / size / width+length / diameter)")


def _conductivity(spec: dict, where: str, local_materials: dict) -> float:
    """熱伝導率 k [W/(m K)] を、直接指定か材料名から得る。"""
    if "k" in spec:
        return float(spec["k"])
    name = _get(spec, "material", "mat")
    if name is None:
        raise KeyError(f"{where}: 熱伝導率が要ります (k: 値 か material: 名前)")
    if str(name).lower() in local_materials:
        return float(local_materials[str(name).lower()]["k"])
    return materials.property_of(name, "k")


def _material_props(name: str, local_materials: dict) -> dict:
    if str(name).lower() in local_materials:
        return local_materials[str(name).lower()]
    return materials.get(name)


def _equivalent_radius(spec: dict, prefix: str, where: str) -> float:
    """円形換算半径を得る。radius / diameter / size(正方形の一辺)/ area から。"""
    if f"{prefix}_radius" in spec:
        return float(spec[f"{prefix}_radius"])
    if f"{prefix}_diameter" in spec:
        return float(spec[f"{prefix}_diameter"]) / 2.0
    if f"{prefix}_size" in spec:       # 正方形 -> 同じ面積の円に換算
        return float(spec[f"{prefix}_size"]) / math.sqrt(math.pi)
    if f"{prefix}_area" in spec:
        return math.sqrt(float(spec[f"{prefix}_area"]) / math.pi)
    raise KeyError(f"{where}: {prefix}_radius / {prefix}_diameter / "
                   f"{prefix}_size / {prefix}_area のいずれかが要ります")


# =====================================================================
# 熱伝達率(定数 or 相関式)
# =====================================================================
def _h_provider(spec: dict, where: str) -> tuple:
    """熱伝達率を返す関数と、非線形かどうかのフラグを作る。

    h: 数値                       -> 定数(線形)
    correlation: 相関式名 + 引数  -> (T_surface, T_fluid) の関数(非線形)

    使える相関式:
      natural_vertical      : 垂直平板の自然対流 (height)
      natural_horizontal_up : 上向き高温水平面   (char_length または area+perimeter)
      forced_plate          : 平板の強制対流     (velocity, length)
      duct                  : ダクト内強制対流   (velocity, hydraulic_diameter)
    """
    if "h" in spec and spec["h"] is not None:
        h_value = float(spec["h"])
        return (lambda ts, tf: h_value), False

    corr = _get(spec, "correlation", "corr")
    if corr is None:
        raise KeyError(f"{where}: h か correlation のどちらかが必要です")
    corr = str(corr).lower()

    if corr in ("natural_vertical", "vertical_plate"):
        height = float(_get(spec, "height", "length", required=True, where=where))
        return (lambda ts, tf: fluids.h_natural_vertical_plate(ts, tf, height)), True

    if corr in ("natural_horizontal_up", "horizontal_plate_up"):
        if "char_length" in spec:
            lc = float(spec["char_length"])
        else:
            lc = _area(spec, where) / float(
                _get(spec, "perimeter", required=True, where=where))
        return (lambda ts, tf: fluids.h_natural_horizontal_plate_up(ts, tf, lc)), True

    if corr in ("forced_plate", "flat_plate"):
        vel = float(_get(spec, "velocity", required=True, where=where))
        length = float(_get(spec, "length", required=True, where=where))
        return (lambda ts, tf: fluids.h_forced_flat_plate(vel, length,
                                                          0.5 * (ts + tf))), True

    if corr == "duct":
        vel = float(_get(spec, "velocity", required=True, where=where))
        d_h = float(_get(spec, "hydraulic_diameter", "d_h", required=True,
                         where=where))
        uniform = bool(spec.get("uniform_flux", True))
        return (lambda ts, tf: fluids.h_duct(vel, d_h, 0.5 * (ts + tf),
                                             uniform)), True

    raise ValueError(f"{where}: 未知の相関式 '{corr}'")


# =====================================================================
# 要素 -> 熱抵抗
# =====================================================================
def _build_element(spec: dict, local_materials: dict):
    """要素の仕様から (熱抵抗, 説明文) を作る。抵抗は数値か (Ta,Tb)->R の関数。"""
    etype = str(_get(spec, "type", required=True, where="element")).lower()
    name = spec.get("name", f"{spec.get('from')}->{spec.get('to')}")
    where = f"element '{name}'"

    # --- 直接指定 ---------------------------------------------------
    if etype in ("resistance", "r"):
        value = float(_get(spec, "value", "resistance", "R", required=True,
                           where=where))
        return value, f"R = {value:.5g} K/W (直接指定)"

    # --- 伝導 -------------------------------------------------------
    if etype in ("conduction", "cond", "plane"):
        k = _conductivity(spec, where, local_materials)
        thickness = float(_get(spec, "thickness", "length", required=True,
                               where=where))
        area = _area(spec, where)
        return (el.r_cond_plane(thickness, k, area),
                f"L/(kA): L={thickness:.4g} m, k={k:.4g}, A={area:.4g} m^2")

    if etype in ("cylinder", "pipe"):
        k = _conductivity(spec, where, local_materials)
        r_in = float(_get(spec, "r_in", "inner_radius", required=True, where=where))
        r_out = float(_get(spec, "r_out", "outer_radius", required=True, where=where))
        length = float(_get(spec, "length", default=1.0))
        return (el.r_cond_cylinder(r_in, r_out, k, length),
                f"ln(r2/r1)/(2 pi k L): r1={r_in:.4g}, r2={r_out:.4g}, k={k:.4g}")

    if etype == "sphere":
        k = _conductivity(spec, where, local_materials)
        r_in = float(_get(spec, "r_in", "inner_radius", required=True, where=where))
        r_out = float(_get(spec, "r_out", "outer_radius", required=True, where=where))
        return (el.r_cond_sphere(r_in, r_out, k),
                f"(1/r1-1/r2)/(4 pi k): r1={r_in:.4g}, r2={r_out:.4g}")

    # --- 界面 -------------------------------------------------------
    if etype in ("interface", "tim"):
        area = _area(spec, where)
        if "resistivity" in spec:
            rho_th = float(spec["resistivity"])
            return (el.r_interface(rho_th, area),
                    f"熱抵抗率/A: {rho_th:.4g} m^2K/W / {area:.4g} m^2")
        # 厚みと k から作る場合(グリース層など)
        k = _conductivity(spec, where, local_materials)
        thickness = float(_get(spec, "thickness", required=True, where=where))
        return (el.r_cond_plane(thickness, k, area),
                f"L/(kA): L={thickness:.4g} m, k={k:.4g}, A={area:.4g} m^2")

    if etype == "contact":
        area = _area(spec, where)
        h_c = float(_get(spec, "conductance", "h_c", required=True, where=where))
        return (el.r_contact(h_c, area), f"1/(h_c A): h_c={h_c:.4g}, A={area:.4g}")

    # --- 対流 -------------------------------------------------------
    if etype in ("convection", "conv"):
        area = _area(spec, where)
        h_fun, nonlinear = _h_provider(spec, where)
        if not nonlinear:
            h = h_fun(0.0, 0.0)
            return el.r_conv(h, area), f"1/(hA): h={h:.4g}, A={area:.4g} m^2"
        desc = f"1/(hA) 相関式 {spec.get('correlation')}, A={area:.4g} m^2"
        return (lambda ts, tf: el.r_conv(max(h_fun(ts, tf), 1e-9), area)), desc

    # --- 輻射 -------------------------------------------------------
    if etype in ("radiation", "rad"):
        area = _area(spec, where)
        eps = float(_get(spec, "emissivity", "eps", required=True, where=where))
        view = float(_get(spec, "view_factor", default=1.0))
        desc = f"1/(h_r A): eps={eps:.3g}, F={view:.3g}, A={area:.4g} m^2"
        return (lambda ts, tsur: el.r_conv(
            max(el.h_radiation(eps, ts, tsur, view), 1e-9), area)), desc

    # --- フィンアレイ -----------------------------------------------
    if etype in ("fin_array", "fins"):
        k = _conductivity(spec, where, local_materials)
        fin_t = float(_get(spec, "fin_thickness", required=True, where=where))
        fin_l = float(_get(spec, "fin_height", "fin_length", required=True,
                           where=where))
        fin_w = float(_get(spec, "fin_width", required=True, where=where))
        n_fins = float(_get(spec, "n_fins", "count", required=True, where=where))
        base_area = _get(spec, "base_area")
        base_area = float(base_area) if base_area is not None else _area(spec, where)
        h_fun, nonlinear = _h_provider(spec, where)

        def _r_fin(h):
            return el.r_fin_array(max(h, 1e-9), k, fin_t, fin_l, fin_w, n_fins,
                                  base_area)

        if not nonlinear:
            h = h_fun(0.0, 0.0)
            eta = el.fin_efficiency(h, k, fin_t, fin_l)
            return _r_fin(h), (f"フィン {n_fins:.0f} 枚, h={h:.4g}, "
                               f"効率 eta={100 * eta:.1f} %")
        return (lambda ts, tf: _r_fin(h_fun(ts, tf))), \
               f"フィン {n_fins:.0f} 枚, 相関式 {spec.get('correlation')}"

    # --- 拡がり抵抗 -------------------------------------------------
    if etype in ("spreading", "spread"):
        k = _conductivity(spec, where, local_materials)
        a = _equivalent_radius(spec, "source", where)
        if "plate_radius" in spec or "plate_diameter" in spec \
                or "plate_size" in spec or "plate_area" in spec:
            b = _equivalent_radius(spec, "plate", where)
            thickness = float(_get(spec, "thickness", required=True, where=where))
            h_raw = _get(spec, "h", "sink_h", required=True, where=where)
            if isinstance(h_raw, str) and h_raw.strip().lower() == "auto":
                # 下流の合成抵抗から h を自己整合的に決める(build_network が反復)
                return (_AutoSpreading(a, b, thickness, k, math.pi * b ** 2),
                        f"Lee et al.: a={a:.4g}, b={b:.4g}, t={thickness:.4g}, "
                        f"k={k:.4g}, h=auto")
            h_sink = float(h_raw)
            return (el.r_spreading_disk(a, b, thickness, k, h_sink),
                    f"Lee et al.: a={a:.4g}, b={b:.4g}, t={thickness:.4g}, "
                    f"k={k:.4g}, h={h_sink:.4g}")
        isothermal = bool(spec.get("isothermal", False))
        return (el.r_spreading_halfspace(a, k, isothermal),
                f"半無限体: a={a:.4g}, k={k:.4g}, "
                f"{'等温' if isothermal else '等熱流束'}")

    # --- 流体の昇温(排熱) -----------------------------------------
    if etype in ("flow", "caloric", "fluid"):
        cp = _get(spec, "cp")
        fluid_name = _get(spec, "fluid", default="air")
        props = _material_props(fluid_name, local_materials)
        cp = float(cp) if cp is not None else float(props["cp"])
        if "mass_flow" in spec:
            mdot = float(spec["mass_flow"])
        elif "volume_flow" in spec or "volumetric_flow" in spec:
            vdot = float(_get(spec, "volume_flow", "volumetric_flow"))
            rho = float(_get(spec, "rho", default=props["rho"]))
            mdot = rho * vdot
        elif "lpm" in spec:            # L/min は現場でよく使う単位
            rho = float(_get(spec, "rho", default=props["rho"]))
            mdot = rho * float(spec["lpm"]) / 60000.0
        elif "cfm" in spec:            # 空冷ファンの単位
            rho = float(_get(spec, "rho", default=props["rho"]))
            mdot = rho * float(spec["cfm"]) * 0.000471947
        else:
            raise KeyError(f"{where}: mass_flow / volume_flow / lpm / cfm "
                           "のいずれかが要ります")
        mean = bool(spec.get("mean", True))
        factor = 2.0 if mean else 1.0
        return (el.r_flow(mdot, cp) / factor,
                f"1/({factor:.0f} mdot cp): mdot={mdot:.5g} kg/s, cp={cp:.4g}"
                f" -> 昇温 {1.0 / (mdot * cp):.4g} K/W"
                f" ({'平均温度基準' if mean else '出口温度基準'})")

    raise ValueError(f"{where}: 未知の要素タイプ '{etype}'")


# =====================================================================
# 拡がり抵抗の h を下流から自動決定する仕組み
# =====================================================================
class _AutoSpreading:
    """h: auto と書かれた拡がり抵抗。

    拡がり抵抗の相関式は「板の裏面から熱が抜ける速さ(= h)」を必要とする。
    ところが実際の h 相当値は、その板より下流にある全抵抗で決まる:

        h_eq = 1 / (R_downstream x A_plate)

    手で反復して合わせるのは面倒なので、回路を組んでから自動で解く。
    ネットワーク構築後に _resolve_auto_spreading() が呼ばれ、
    枝の抵抗値が数回の反復で自己整合値に置き換わる。
    """

    def __init__(self, a, b, thickness, k, plate_area):
        self.a = a
        self.b = b
        self.thickness = thickness
        self.k = k
        self.plate_area = plate_area
        self.h = 1000.0                      # 初期推定
        self.value = el.r_spreading_disk(a, b, thickness, k, self.h)

    def update(self, r_downstream: float) -> float:
        """下流の合成抵抗 [K/W] から h と拡がり抵抗を更新し、変化量を返す。"""
        self.h = max(1.0 / (max(r_downstream, 1e-12) * self.plate_area), 1e-3)
        new_value = el.r_spreading_disk(self.a, self.b, self.thickness,
                                        self.k, self.h)
        delta = abs(new_value - self.value)
        self.value = new_value
        return delta


def _downstream_resistance(net: ThermalNetwork, br_index: int,
                           sinks: List[str]) -> float:
    """指定した枝を外したとき、その下流側ノードから境界までの合成抵抗 [K/W]。

    枝を外すと上流側が回路から切り離されて行列が特異になるので、
    「下流側から到達できるノードだけ」を取り出した部分回路で計算する。
    """
    branch = net.branches[br_index]
    others = [b for i, b in enumerate(net.branches) if i != br_index]
    probe = ThermalNetwork("probe")
    probe.branches = others
    reachable = probe.connected_nodes(branch.b)
    sink = next((s for s in sinks if s in reachable), None)
    if sink is None:
        raise ValueError(
            f"要素 '{branch.label}' の下流に温度指定ノードがありません。"
            "h: auto は下流の合成抵抗から h を決めるので、"
            "その先が境界につながっている必要があります。")
    pruned = ThermalNetwork("downstream")
    for name in reachable:
        pruned.add_node(name, fixed=(net.nodes[name].fixed if name == sink else None))
    for b in others:
        if b.a in reachable and b.b in reachable:
            pruned.add_resistor(b.a, b.b, b.resistance, b.label)
    t_sink = float(net.nodes[sink].fixed)
    return pruned.total_resistance(branch.b, sink, sink_temperature=t_sink)


def _resolve_auto_spreading(net: ThermalNetwork, descriptions: dict,
                            max_iter: int = 30, tol: float = 1e-10) -> None:
    """h: auto の拡がり抵抗を、下流の合成抵抗と整合するまで反復する。

    descriptions["_auto"] は {枝のインデックス: _AutoSpreading} の対応表。
    各反復で h_eq = 1/(R_downstream x A_plate) として抵抗を作り直す。
    """
    autos = descriptions.pop("_auto", None)
    if not autos:
        return
    sinks = [n for n, node in net.nodes.items() if node.fixed is not None]
    if not sinks:
        raise ValueError("h: auto を使うには温度指定ノードが必要です")

    for iterations in range(1, max_iter + 1):
        delta_max = 0.0
        for br_index, auto in autos.items():
            branch = net.branches[br_index]
            r_down = _downstream_resistance(net, br_index, sinks)
            delta_max = max(delta_max, auto.update(r_down))
            net.branches[br_index].resistance = auto.value
            descriptions[branch.label] = (
                f"Lee et al.: a={auto.a:.4g}, b={auto.b:.4g}, "
                f"t={auto.thickness:.4g}, k={auto.k:.4g}, "
                f"h={auto.h:.5g} W/(m^2 K) <- 下流 {r_down:.5g} K/W から自動決定 "
                f"({iterations} 回反復)")
        if delta_max < tol:
            break


# =====================================================================
# ネットワーク構築
# =====================================================================
def build_network(spec: dict) -> tuple:
    """モデル仕様(dict)から ThermalNetwork と要素の説明を作る。"""
    local_materials = {str(k).lower(): v
                       for k, v in (spec.get("materials") or {}).items()}
    net = ThermalNetwork(spec.get("name", "model"))
    descriptions: Dict[str, str] = {}

    # --- ノード -----------------------------------------------------
    for name, node_spec in (spec.get("nodes") or {}).items():
        node_spec = node_spec or {}
        cap = node_spec.get("capacitance")
        if cap is None and "material" in node_spec:
            props = _material_props(node_spec["material"], local_materials)
            if "mass" in node_spec:
                cap = float(node_spec["mass"]) * props["cp"]
            elif "volume" in node_spec:
                cap = el.capacitance(props["rho"], props["cp"],
                                     float(node_spec["volume"]))
            else:
                raise KeyError(f"node '{name}': material を書いたら "
                               "volume か mass も要ります")
        heat = node_spec.get("heat", 0.0)
        if "heat_schedule" in node_spec:      # 過渡用の時間変化する発熱
            heat = _schedule(node_spec["heat_schedule"], name)
        net.add_node(name, capacitance=float(cap or 0.0),
                     fixed=(None if node_spec.get("fixed") is None
                            else float(node_spec["fixed"])),
                     heat=heat)

    # --- 要素 -------------------------------------------------------
    for i, element in enumerate(spec.get("elements") or []):
        if element.get("disabled"):
            continue          # 一時的に外して効果を見るためのスイッチ
        a = _get(element, "from", "a", required=True, where=f"element[{i}]")
        b = _get(element, "to", "b", required=True, where=f"element[{i}]")
        label = element.get("name", f"{a}->{b}")
        resistance, desc = _build_element(element, local_materials)
        for node in (a, b):
            if node not in net.nodes:
                net.add_node(node)
        if isinstance(resistance, _AutoSpreading):
            descriptions.setdefault("_auto", {})[len(net.branches)] = resistance
            net.add_resistor(a, b, resistance.value, label)
        else:
            net.add_resistor(a, b, resistance, label)
        descriptions[label] = desc

    validate(net)
    _resolve_auto_spreading(net, descriptions)
    return net, descriptions


def _schedule(rows, name: str) -> Callable[[float], float]:
    """[[t0, Q0], [t1, Q1], ...] を階段状(ゼロ次ホールド)の Q(t) にする。"""
    pairs = sorted((float(t), float(q)) for t, q in rows)
    if not pairs:
        raise ValueError(f"node '{name}': heat_schedule が空です")
    times = np.array([p[0] for p in pairs])
    powers = np.array([p[1] for p in pairs])

    def q_of_t(t: float) -> float:
        idx = int(np.searchsorted(times, t, side="right") - 1)
        return float(powers[max(idx, 0)])

    return q_of_t


def validate(net: ThermalNetwork) -> None:
    """よくある作り間違いを検出する。

    * 温度固定ノードが 1 つも無い   -> 解が定まらない
    * 温度固定ノードにつながらない孤立ノードがある -> 浮きノード
    """
    fixed = [n for n, node in net.nodes.items() if node.fixed is not None]
    if not fixed:
        raise ValueError("温度を固定したノードがありません。"
                         "外気・冷却水などの境界に fixed: 値 を設定してください。")
    adjacency: Dict[str, set] = {n: set() for n in net.nodes}
    for br in net.branches:
        adjacency[br.a].add(br.b)
        adjacency[br.b].add(br.a)
    seen, stack = set(fixed), list(fixed)
    while stack:
        current = stack.pop()
        for nb in adjacency[current]:
            if nb not in seen:
                seen.add(nb)
                stack.append(nb)
    floating = sorted(set(net.nodes) - seen)
    if floating:
        raise ValueError(
            f"温度固定ノードにつながっていないノードがあります: {floating}\n"
            "  熱の逃げ道がないので温度が決まりません(発熱があれば発散します)。")


# =====================================================================
# 読み込みと上書き
# =====================================================================
def load_spec(path: str) -> dict:
    """YAML または JSON のモデルファイルを読む。"""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    if path.lower().endswith((".yaml", ".yml")):
        import yaml
        return yaml.safe_load(text)
    return json.loads(text)


def apply_override(spec: dict, assignment: str) -> dict:
    """'nodes.die.heat=120' や 'elements.TIM2.resistivity=1e-5' を適用する。

    elements はリストなので name で引けるようにしている
    (ファイルを書き換えずに条件を振れる = 感度分析が一瞬でできる)。
    """
    if "=" not in assignment:
        raise ValueError(f"--set の書式は key.path=値 です: {assignment}")
    path, raw = assignment.split("=", 1)
    keys = path.split(".")
    try:
        value: Any = json.loads(raw)
    except json.JSONDecodeError:
        value = raw
    target: Any = spec
    for depth, key in enumerate(keys[:-1]):
        if isinstance(target, list):
            match = [x for x in target if str(x.get("name")) == key]
            if not match:
                raise KeyError(f"--set {path}: 要素 '{key}' が見つかりません")
            target = match[0]
        else:
            if key not in target:
                raise KeyError(f"--set {path}: '{key}' が見つかりません "
                               f"(深さ {depth}, 候補 {sorted(target)})")
            target = target[key]
    if isinstance(target, list):
        raise KeyError(f"--set {path}: リストの直下には代入できません")
    target[keys[-1]] = value
    return spec


# =====================================================================
# レポート
# =====================================================================
def report_steady(net: ThermalNetwork, descriptions: dict, sol,
                  source: str | None = None) -> str:
    lines = []
    lines.append(f"モデル: {net.name}")
    lines.append(f"  反復 {sol.iterations} 回 / 残差 {sol.residual:.2e} K")
    lines.append("")
    lines.append("  ノード温度")
    lines.append("    " + f"{'node':<20s}{'T [degC]':>12s}  {'条件':<28s}")
    lines.append("    " + "-" * 62)
    for name, node in net.nodes.items():
        cond = []
        if node.fixed is not None:
            cond.append("温度指定")
        if node.heat:
            cond.append("発熱 可変" if callable(node.heat)
                        else f"発熱 {node.heat:.4g} W")
        if node.capacitance:
            cond.append(f"C={node.capacitance:.4g} J/K")
        lines.append(f"    {name:<20s}{sol[name]:12.2f}  {', '.join(cond):<28s}")

    total_heat = sum((n.heat(0.0) if callable(n.heat) else n.heat)
                     for n in net.nodes.values())
    lines.append("")
    lines.append("  要素(熱抵抗の内訳)")
    lines.append("    " + f"{'element':<24s}{'R [K/W]':>10s}{'dT [K]':>9s}"
                          f"{'Q [W]':>9s}{'dT share':>10s}")
    lines.append("    " + "-" * 62)
    rows = []
    for br in net.branches:
        ta, tb = sol[br.a], sol[br.b]
        r = br.value(ta, tb)
        q = (ta - tb) / r
        rows.append((br.label or f"{br.a}->{br.b}", r, ta - tb, q))
    dt_sum = sum(abs(row[2]) for row in rows) or 1.0
    for label, r, dt, q in rows:
        lines.append(f"    {label:<24s}{r:10.4f}{dt:9.2f}{q:9.2f}"
                     f"{100 * abs(dt) / dt_sum:9.1f}%")
    if len({(br.a, br.b) for br in net.branches}) != len(net.branches):
        lines.append("    ※ 並列枝があるため dT share の合計は 100 % になりません")

    # エネルギー収支: 発熱の合計と、温度指定ノードから出入りする熱の合計が釣り合う
    lines.append("")
    lines.append("  エネルギー収支")
    lines.append(f"    {'発熱の合計':<24s}{total_heat:10.4f} W")
    boundary_total = 0.0
    for name, node in net.nodes.items():
        if node.fixed is None:
            continue
        inflow = 0.0        # 回路 -> この境界ノード に流れ込む熱 [W]
        for br in net.branches:
            if br.a == name:
                inflow -= (sol[br.a] - sol[br.b]) / br.value(sol[br.a], sol[br.b])
            elif br.b == name:
                inflow += (sol[br.a] - sol[br.b]) / br.value(sol[br.a], sol[br.b])
        boundary_total += inflow
        arrow = "回路 -> 境界" if inflow >= 0 else "境界 -> 回路"
        lines.append(f"    {name + ' (温度指定)':<24s}{inflow:10.4f} W  {arrow}")
    lines.append(f"    {'残差':<24s}{abs(total_heat - boundary_total):10.2e} W")

    if source is None:
        heats = [n for n, node in net.nodes.items() if node.heat]
        source = heats[0] if len(heats) == 1 else None
    if source:
        sinks = [n for n, node in net.nodes.items() if node.fixed is not None]
        q_src = net.nodes[source].heat
        q_src = q_src(0.0) if callable(q_src) else q_src
        if q_src:
            lines.append(f"  合成熱抵抗 R({source} -> {sinks[0]}) = "
                         f"{(sol[source] - sol[sinks[0]]) / q_src:.4f} K/W")

    if descriptions:
        lines.append("")
        lines.append("  要素の内訳(どの物理量から抵抗が出たか)")
        for label, desc in descriptions.items():
            lines.append(f"    {label:<24s} {desc}")
    return "\n".join(lines)


# =====================================================================
# CLI
# =====================================================================
def _parse_sweep(text: str) -> tuple:
    """'nodes.die.heat=50:150:11' -> (パス, [値, ...])"""
    path, spec = text.split("=", 1)
    parts = spec.split(":")
    if len(parts) == 3:
        lo, hi, num = float(parts[0]), float(parts[1]), int(parts[2])
        return path, list(np.linspace(lo, hi, num))
    return path, [float(v) for v in spec.split(",")]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="YAML/JSON のモデルを解く 1D 熱回路ソルバ",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""例:
  python -m thermalnet.model models/cpu_heatsink.yaml
  python -m thermalnet.model models/cpu_heatsink.yaml --set nodes.die.heat=120
  python -m thermalnet.model models/cpu_heatsink.yaml \\
      --sweep nodes.die.heat=50:150:11 --watch die
  python -m thermalnet.model models/enclosure.yaml --transient --t-end 1800
""")
    parser.add_argument("model", help="モデルファイル (.yaml / .json)")
    parser.add_argument("--set", action="append", default=[], metavar="PATH=VALUE",
                        help="値を上書きする (例: nodes.die.heat=120)")
    parser.add_argument("--sweep", metavar="PATH=LO:HI:N",
                        help="パラメータを振って表にする")
    parser.add_argument("--watch", action="append", default=[], metavar="NODE",
                        help="スイープ/過渡で追跡するノード")
    parser.add_argument("--transient", action="store_true", help="過渡解析を行う")
    parser.add_argument("--t-end", type=float, help="過渡の終了時刻 [s]")
    parser.add_argument("--dt", type=float, help="時間刻み [s]")
    parser.add_argument("--initial", type=float, help="初期温度 [degC]")
    parser.add_argument("--guess", type=float, default=None, help="初期推定温度")
    parser.add_argument("--plot", metavar="PNG", help="結果を図に保存")
    parser.add_argument("--quiet", action="store_true", help="要素内訳を省略")
    args = parser.parse_args(argv)

    spec = load_spec(args.model)
    for assignment in args.set:
        apply_override(spec, assignment)

    solve_opts = spec.get("solve") or {}
    guess = args.guess if args.guess is not None else solve_opts.get("guess", 25.0)

    if args.sweep:
        return _run_sweep(spec, args, guess)

    net, descriptions = build_network(spec)

    if args.transient:
        return _run_transient(net, spec, args, solve_opts)

    sol = net.solve_steady(guess=guess,
                           relaxation=float(solve_opts.get("relaxation", 0.7)))
    print(report_steady(net, {} if args.quiet else descriptions, sol,
                        solve_opts.get("source")))
    if args.plot:
        _plot_steady(net, sol, args.plot)
    return 0


def _run_sweep(spec: dict, args, guess: float) -> int:
    path, values = _parse_sweep(args.sweep)
    watch = args.watch or None
    print(f"スイープ: {path}")
    header_done = False
    for value in values:
        local = json.loads(json.dumps(spec))       # 深いコピー
        apply_override(local, f"{path}={value}")
        net, _ = build_network(local)
        sol = net.solve_steady(guess=guess)
        names = watch or [n for n, node in net.nodes.items() if node.fixed is None]
        if not header_done:
            print("    " + f"{'value':>12s}" + "".join(f"{n:>14s}" for n in names))
            print("    " + "-" * (12 + 14 * len(names)))
            header_done = True
        print("    " + f"{value:12.4g}" + "".join(f"{sol[n]:14.2f}" for n in names))
    return 0


def _run_transient(net: ThermalNetwork, spec: dict, args, solve_opts: dict) -> int:
    tr_opts = spec.get("transient") or {}
    t_end = args.t_end or tr_opts.get("t_end")
    dt = args.dt or tr_opts.get("dt")
    if t_end is None or dt is None:
        print("過渡解析には --t-end と --dt (またはモデルの transient: 節) が要ります",
              file=sys.stderr)
        return 2
    initial = args.initial if args.initial is not None else tr_opts.get("initial", 25.0)
    if not any(node.capacitance for node in net.nodes.values()):
        print("警告: 熱容量が 1 つも設定されていません。"
              "過渡応答は出ません(定常解が返ります)", file=sys.stderr)
    store = max(1, int(round((t_end / dt) / 2000)))
    tr = net.solve_transient(float(t_end), float(dt), initial=initial,
                             theta=float(tr_opts.get("theta", 1.0)),
                             store_every=store)
    names = args.watch or [n for n, node in net.nodes.items() if node.fixed is None]
    print(f"過渡解析: 0 -> {t_end} s (dt = {dt} s)")
    print("    " + f"{'t [s]':>10s}" + "".join(f"{n:>14s}" for n in names))
    print("    " + "-" * (10 + 14 * len(names)))
    for i in np.linspace(0, len(tr.time) - 1, min(11, len(tr.time))).astype(int):
        print("    " + f"{tr.time[i]:10.2f}"
              + "".join(f"{tr[n][i]:14.2f}" for n in names))
    peaks = {n: float(np.max(tr[n])) for n in names}
    print("\n    最高温度: " + ", ".join(f"{n} {v:.2f} degC"
                                          for n, v in peaks.items()))
    if args.plot:
        _plot_transient(tr, names, args.plot)
    return 0


def _plot_steady(net: ThermalNetwork, sol, path: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels, drops = [], []
    for br in net.branches:
        ta, tb = sol[br.a], sol[br.b]
        labels.append(br.label or f"{br.a}->{br.b}")
        drops.append(br.value(ta, tb))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    axes[0].barh(labels, drops, color="#4C72B0")
    axes[0].set_xlabel("thermal resistance [K/W]")
    axes[0].set_title("Resistance breakdown")
    axes[0].grid(alpha=0.3, axis="x")
    names = list(net.nodes)
    axes[1].bar(names, [sol[n] for n in names], color="#55A868")
    axes[1].set_ylabel("temperature [degC]")
    axes[1].set_xticks(range(len(names)), names, rotation=45, ha="right", fontsize=8)
    axes[1].set_title("Node temperatures")
    axes[1].grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"\n[図を保存] {path}")


def _plot_transient(tr, names, path: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.6))
    for name in names:
        ax.plot(tr.time, tr[name], lw=2, label=name)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("temperature [degC]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"\n[図を保存] {path}")


if __name__ == "__main__":
    sys.exit(main())
