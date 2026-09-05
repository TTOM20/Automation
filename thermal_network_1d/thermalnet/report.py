"""モデルファイルから熱設計レポート(HTML)を自動生成する。

    python -m thermalnet.report models/cpu_heatsink.yaml -o report.html
    python -m thermalnet.report models/power_module_liquid.yaml --limit 125

レポートに載るもの:
    判定サマリ / 解析条件 / 主熱流経路の温度分布 / 熱抵抗の内訳 /
    感度分析 / 過渡応答 / 検算 / 前提と適用限界 / モデル定義(付録)

「解いた結果」だけでなく「検算」と「適用限界」を必ず書くのがこのレポートの型。
1D モデルは数字がいくらでも出てしまうので、どこまで信じてよいかを
明示しない限りレポートとして成立しない。
"""

from __future__ import annotations

import argparse
import copy
import datetime as _dt
import html
import os
import sys
from typing import Dict, List, Tuple

import numpy as np

from . import charts
from .network import ThermalNetwork
from .model import build_network, load_spec

DEFAULT_STYLE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "report_style.css")

# レポート内の小さなツールチップ(表と併記しているので補助的な役割)
_PATH_CHAIN = object()      # 検算表で「主熱流経路」の行だけ HTML を差し込むための印

TOOLTIP_JS = """<script>
// グラフのマークにカーソルを乗せたら値を出す(表と併記しているので補助的な役割)
(function () {
  var tip = document.getElementById("tip");
  function show(event) {
    var text = event.target.getAttribute("data-tip");
    if (!text) { return; }
    tip.textContent = text;
    tip.classList.add("on");
    move(event);
  }
  function move(event) {
    var pad = 14;
    var x = Math.min(event.clientX + pad, window.innerWidth - tip.offsetWidth - 8);
    var y = Math.max(event.clientY - tip.offsetHeight - pad, 8);
    tip.style.left = x + "px";
    tip.style.top = y + "px";
  }
  function hide() { tip.classList.remove("on"); }
  document.addEventListener("mouseover", show, true);
  document.addEventListener("mousemove", function (e) {
    if (tip.classList.contains("on")) { move(e); }
  }, true);
  document.addEventListener("mouseout", hide, true);
})();
</script>"""


def _esc(text) -> str:
    return html.escape(str(text), quote=False)


# =====================================================================
# 解析
# =====================================================================
def _heat_of(node, time: float = 0.0) -> float:
    return node.heat(time) if callable(node.heat) else float(node.heat)


def main_heat_path(net: ThermalNetwork, temps: Dict[str, float],
                   start: str, sinks: List[str]) -> Tuple[List[str], List[str]]:
    """主熱流経路: 最も熱流の大きい枝をたどって境界まで降りる。

    分岐がある回路でも「どこを通って熱が逃げているか」が一本の線で読める。
    """
    path, labels = [start], []
    visited = {start}
    current = start
    for step in range(len(net.nodes) + 1):
        if current in sinks and step > 0:
            break
        best, best_q = None, 0.0
        for br in net.branches:
            if br.a == current:
                other = br.b
            elif br.b == current:
                other = br.a
            else:
                continue
            if other in visited:
                continue
            r = br.value(temps[br.a], temps[br.b])
            q = (temps[current] - temps[other]) / r
            if q > best_q:
                best, best_q = (br, other), q
        if best is None:
            break
        br, other = best
        path.append(other)
        labels.append(br.label or f"{br.a}->{br.b}")
        visited.add(other)
        current = other
    return path, labels


def _scaled(resistance, factor: float):
    """抵抗を factor 倍する(温度依存の枝でも扱えるようにする)。"""
    if callable(resistance):
        return lambda ta, tb: resistance(ta, tb) * factor
    return resistance * factor


def sensitivity(spec: dict, target: str, base_temp: float,
                guess: float) -> List[Tuple[str, float]]:
    """各要素の熱抵抗を -20 % したときの目標ノード温度の変化 [K]。

    線形な回路では内訳表と等価な情報になるが、非線形要素(自然対流・輻射)が
    あると順位が変わるので、レポートには実際に解いた値を載せる。
    """
    results = []
    net, _ = build_network(copy.deepcopy(spec))
    for index, branch in enumerate(net.branches):
        trial, _ = build_network(copy.deepcopy(spec))
        trial.branches[index].resistance = _scaled(
            trial.branches[index].resistance, 0.8)
        try:
            t = trial.solve_steady(guess=guess)[target]
        except Exception:                                   # noqa: BLE001
            continue
        label = branch.label or f"{branch.a}->{branch.b}"
        results.append((f"{label} を -20 %", t - base_temp))
    return sorted(results, key=lambda row: row[1])


def condition_sensitivity(spec: dict, target: str, base_temp: float,
                          guess: float) -> List[Tuple[str, float]]:
    """運転条件(発熱・境界温度)を振ったときの温度変化 [K]。"""
    rows = []
    net, _ = build_network(copy.deepcopy(spec))
    heat_nodes = [n for n, node in net.nodes.items() if _heat_of(node)]
    fixed_nodes = [n for n, node in net.nodes.items() if node.fixed is not None]

    for name in heat_nodes:
        trial_spec = copy.deepcopy(spec)
        node_spec = trial_spec["nodes"][name]
        if "heat" in node_spec:
            node_spec["heat"] = float(node_spec["heat"]) * 1.2
        elif "heat_schedule" in node_spec:
            node_spec["heat_schedule"] = [[t, q * 1.2]
                                          for t, q in node_spec["heat_schedule"]]
        else:
            continue
        trial, _ = build_network(trial_spec)
        rows.append((f"{name} の発熱 +20 %",
                     trial.solve_steady(guess=guess)[target] - base_temp))

    for name in fixed_nodes:
        trial_spec = copy.deepcopy(spec)
        trial_spec["nodes"][name]["fixed"] = \
            float(trial_spec["nodes"][name]["fixed"]) + 5.0
        trial, _ = build_network(trial_spec)
        rows.append((f"{name} の温度 +5 K",
                     trial.solve_steady(guess=guess)[target] - base_temp))
    return rows


def analyze(spec: dict, limit: float | None = None,
            target: str | None = None) -> dict:
    """モデルを解き、レポートに必要な情報を全部そろえる。"""
    report_opts = spec.get("report") or {}
    solve_opts = spec.get("solve") or {}
    guess = float(solve_opts.get("guess", 25.0))
    limit = limit if limit is not None else report_opts.get("limit")
    target = target or report_opts.get("target")

    net, descriptions = build_network(copy.deepcopy(spec))
    sol = net.solve_steady(guess=guess,
                           relaxation=float(solve_opts.get("relaxation", 0.7)))
    temps = dict(sol.temperatures)

    sinks = [n for n, node in net.nodes.items() if node.fixed is not None]
    heat_nodes = [n for n, node in net.nodes.items() if _heat_of(node)]
    total_heat = sum(_heat_of(node) for node in net.nodes.values())

    # 境界ノードに出入りする熱量(正 = 回路から境界へ)
    boundary_flow = {}
    for name in sinks:
        inflow = 0.0
        for br in net.branches:
            r = br.value(temps[br.a], temps[br.b])
            if br.a == name:
                inflow -= (temps[br.a] - temps[br.b]) / r
            elif br.b == name:
                inflow += (temps[br.a] - temps[br.b]) / r
        boundary_flow[name] = inflow
    boundary = sum(boundary_flow.values())

    if target is None:
        # 温度指定ノードは評価対象にしない(境界の温度は入力であって結果ではない)
        candidates = [n for n in temps if net.nodes[n].fixed is None] or list(temps)
        target = max(candidates, key=lambda n: temps[n])

    if total_heat > 1e-12:
        # 発熱で駆動される系: 発熱ノードを起点に、投入熱量を基準にする
        path_start = max(heat_nodes, key=lambda n: temps[n]) if heat_nodes else target
        reference_heat = total_heat
        heat_label = "投入熱量"
        heat_note = "、".join(heat_nodes)
        t_sink = float(net.nodes[sinks[0]].fixed)
        r_total = (temps[target] - t_sink) / total_heat
        r_note = f"{target} → {sinks[0]}"
    else:
        # 温度差で駆動される系(壁など): 高温側境界を起点に、通過熱量を基準にする
        hot = max(sinks, key=lambda n: float(net.nodes[n].fixed))
        cold = min(sinks, key=lambda n: float(net.nodes[n].fixed))
        path_start = hot
        reference_heat = abs(boundary_flow[cold])
        heat_label = "通過熱量"
        heat_note = f"{hot} → {cold}"
        t_sink = float(net.nodes[cold].fixed)
        dt_boundary = float(net.nodes[hot].fixed) - t_sink
        r_total = dt_boundary / reference_heat if reference_heat else float("nan")
        r_note = f"{hot} → {cold}"

    path, path_labels = main_heat_path(net, temps, path_start, sinks)

    rows = []
    for br in net.branches:
        r = br.value(temps[br.a], temps[br.b])
        dt = temps[br.a] - temps[br.b]
        rows.append({
            "label": br.label or f"{br.a}->{br.b}",
            "from": br.a, "to": br.b,
            "R": r, "dT": dt, "Q": dt / r,
            "on_path": (br.label or f"{br.a}->{br.b}") in set(path_labels),
            "nonlinear": callable(br.resistance),
            "note": descriptions.get(br.label or "", ""),
        })
    rows.sort(key=lambda row: row["R"], reverse=True)

    # 過渡
    transient = None
    tr_opts = spec.get("transient") or {}
    has_capacitance = any(node.capacitance for node in net.nodes.values())
    if tr_opts.get("t_end") and tr_opts.get("dt") and has_capacitance:
        tr_net, _ = build_network(copy.deepcopy(spec))
        steps = int(round(float(tr_opts["t_end"]) / float(tr_opts["dt"])))
        tr = tr_net.solve_transient(
            float(tr_opts["t_end"]), float(tr_opts["dt"]),
            initial=float(tr_opts.get("initial", guess)),
            theta=float(tr_opts.get("theta", 1.0)),
            store_every=max(1, steps // 400))
        # 時定数の違う 3 点を選ぶ: 目標ノード / 熱容量が最大のノード / 経路末端(流体側)
        watch = [target]
        capacitive = sorted((n for n, node in net.nodes.items()
                             if node.capacitance and n != target),
                            key=lambda n: net.nodes[n].capacitance, reverse=True)
        if capacitive:
            watch.append(capacitive[0])
        tail = [n for n in reversed(path) if n not in sinks and n not in watch]
        if tail:
            watch.append(tail[0])
        watch = watch[:3]
        transient = {
            "time": tr.time,
            "series": [(n, tr[n]) for n in watch],
            "peaks": {n: (float(np.max(tr[n])),
                          float(tr.time[int(np.argmax(tr[n]))])) for n in watch},
            "opts": tr_opts,
            "schedule": next((spec["nodes"][n].get("heat_schedule")
                              for n in heat_nodes
                              if spec.get("nodes", {}).get(n, {}).get(
                                  "heat_schedule")), None),
        }

    return {
        "spec": spec, "net": net, "sol": sol, "temps": temps,
        "descriptions": descriptions, "rows": rows,
        "target": target, "sinks": sinks, "heat_nodes": heat_nodes,
        "total_heat": total_heat, "t_sink": t_sink, "r_total": r_total,
        "reference_heat": reference_heat, "heat_label": heat_label,
        "heat_note": heat_note, "r_note": r_note, "path_start": path_start,
        "boundary_flow": boundary_flow,
        "limit": float(limit) if limit is not None else None,
        "path": path, "path_labels": path_labels,
        "boundary_heat": boundary, "residual": abs(total_heat - boundary),
        "iterations": sol.iterations, "nonlinear": any(r["nonlinear"] for r in rows),
        "sensitivity": sensitivity(spec, target, temps[target], guess),
        "conditions": condition_sensitivity(spec, target, temps[target], guess),
        "transient": transient,
    }


# =====================================================================
# 前提と適用限界(モデルの中身から自動生成)
# =====================================================================
def assumptions(analysis: dict) -> List[Tuple[str, str]]:
    spec = analysis["spec"]
    items: List[Tuple[str, str]] = []
    types = [str(e.get("type", "")).lower() for e in spec.get("elements") or []
             if not e.get("disabled")]

    items.append((
        "1 次元化",
        "各部位を等温の 1 ノードとして扱っている。物体内部に無視できない温度分布が"
        "ある場合(ビオ数 Bi = hL/k が 0.1 を超える場合)は、その部位を複数ノードに"
        "分割する必要がある。"))

    if any(t in ("convection", "conv", "fin_array", "fins") for t in types):
        corr = sorted({str(e.get("correlation")) for e in spec.get("elements") or []
                       if e.get("correlation")})
        if corr:
            items.append((
                "対流の相関式",
                f"熱伝達率は相関式({', '.join(corr)})から求めている。"
                "相関式の適用範囲(流れの向き・レイノルズ数/レイリー数・"
                "助走区間の有無)を外れると誤差が大きくなる。"
                "実機形状での検証には CFD か実測が要る。"))
        else:
            items.append((
                "熱伝達率",
                "熱伝達率は定数として与えている。実機では流れの偏りや"
                "バイパスによって局所的に 2 倍以上変わりうる。"
                "この値の根拠(CFD・実測・カタログ)を明記すること。"))

    if any(t in ("spreading", "spread") for t in types):
        items.append((
            "拡がり抵抗",
            "小さい熱源から広い板へ熱が広がる 3 次元効果を、相関式"
            "(Lee et al., 1995)で 1 本の抵抗に置き換えている。"
            "この相関式は軸対称・等熱流束源・裏面一様 h を前提としており、"
            "数値解との差は 10 % 程度ある。熱源が複数ある場合や"
            "偏在する場合は別途評価が要る。"))

    if any(t in ("flow", "caloric", "fluid") for t in types):
        basis = [e for e in spec.get("elements") or []
                 if str(e.get("type", "")).lower() in ("flow", "caloric", "fluid")]
        mean = all(e.get("mean", True) for e in basis)
        items.append((
            "流体の昇温",
            "流路を通る流体の昇温を 1/(2 m cp) "
            f"({'平均温度基準' if mean else '出口温度基準'})で表している。"
            "流量は入力値であり、ファン/ポンプの動作点(圧損との交点)は"
            "このモデルには含まれていない。流路形状を変えたら流量も見直すこと。"))

    if any(t in ("radiation", "rad") for t in types):
        items.append((
            "輻射",
            "輻射は等価熱伝達率 h_r で線形化し、放射率と形態係数は入力値。"
            "周囲を大空間とみなしており、近接した物体間の相互反射は含まない。"))

    if any(t in ("interface", "tim", "contact") for t in types):
        items.append((
            "界面熱抵抗",
            "TIM・接触部の熱抵抗はカタログ値(熱抵抗率)を面積で割った値。"
            "実装圧力・ボイド・経時劣化で 2 倍以上悪化しうるので、"
            "量産設計ではワーストケースを別途確認すること。"))

    if analysis["transient"] is None:
        items.append((
            "定常のみ",
            "本レポートは定常解が中心。起動時・バースト負荷のピーク温度は"
            "熱容量に依存するため、必要なら過渡解析を追加すること。"))
    return items


# =====================================================================
# HTML 生成
# =====================================================================
def _table(headers: List[str], rows: List[List[str]],
           align_right: List[int] | None = None, caption: str = "") -> str:
    """表を組む。align_right に入れた列は右寄せ + 等幅数字にする。"""
    right = set(align_right or [])
    num_attr = ' class="num"'
    head = "".join(f"<th{num_attr if i in right else ''}>{_esc(h)}</th>"
                   for i, h in enumerate(headers))
    body = []
    for row in rows:
        cells = "".join(f"<td{num_attr if i in right else ''}>{cell}</td>"
                        for i, cell in enumerate(row))
        body.append(f"<tr>{cells}</tr>")
    cap = f"<figcaption>{_esc(caption)}</figcaption>" if caption else ""
    return (f'<figure class="table-wrap">{cap}<table><thead><tr>{head}</tr>'
            f'</thead><tbody>{"".join(body)}</tbody></table></figure>')


def render_html(analysis: dict, title: str) -> str:
    spec = analysis["spec"]
    net = analysis["net"]
    temps = analysis["temps"]
    target = analysis["target"]
    limit = analysis["limit"]
    t_target = temps[target]
    generated = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")

    # --- 判定 --------------------------------------------------------
    if limit is None:
        verdict, verdict_class, margin_text = "判定基準なし", "neutral", "—"
    else:
        margin = limit - t_target
        margin_text = f"{margin:+.1f} K"
        if margin >= 10:
            verdict, verdict_class = "余裕あり", "good"
        elif margin >= 0:
            verdict, verdict_class = "余裕僅少", "warning"
        else:
            verdict, verdict_class = "上限超過", "critical"

    # --- サマリタイル ------------------------------------------------
    tiles = [
        ("最高温度", f"{t_target:.1f}", "degC", f"{target}"),
        ("合成熱抵抗", f"{analysis['r_total']:.3f}", "K/W", analysis["r_note"]),
        (analysis["heat_label"], f"{analysis['reference_heat']:.4g}", "W",
         analysis["heat_note"]),
        ("上限までの余裕", margin_text, "",
         f"上限 {limit:.0f} degC" if limit is not None else "上限未設定"),
    ]
    tile_html = "".join(
        f'<div class="tile"><div class="tile-label">{_esc(label)}</div>'
        f'<div class="tile-value">{_esc(value)}'
        f'<span class="tile-unit">{_esc(unit)}</span></div>'
        f'<div class="tile-note">{_esc(note)}</div></div>'
        for label, value, unit, note in tiles)

    # --- §1 解析条件 -------------------------------------------------
    node_rows = []
    for name, node in net.nodes.items():
        heat = _heat_of(node)
        node_rows.append([
            _esc(name),
            f"{heat:.4g}" if heat else "—",
            f"{node.fixed:.1f}" if node.fixed is not None else "—",
            f"{node.capacitance:.4g}" if node.capacitance else "—",
            f"{temps[name]:.2f}",
        ])
    conditions_table = _table(
        ["ノード", "発熱 [W]", "温度指定 [degC]", "熱容量 [J/K]", "解 [degC]"],
        node_rows, align_right=[1, 2, 3, 4],
        caption="表 1: ノードの条件と解")

    # --- §2 温度分布 -------------------------------------------------
    path = analysis["path"]
    ladder = charts.temperature_ladder(
        path, [temps[n] for n in path], analysis["path_labels"])
    path_rows = []
    for i, label in enumerate(analysis["path_labels"]):
        drop = temps[path[i]] - temps[path[i + 1]]
        share = 100 * drop / (temps[path[0]] - temps[path[-1]] or 1.0)
        path_rows.append([
            f"{i + 1}", _esc(label), f"{path[i]} → {path[i + 1]}",
            f"{drop:.2f}", f"{share:.1f}"])
    path_table = _table(
        ["#", "区間", "ノード", "温度降下 [K]", "経路内の割合 [%]"],
        path_rows, align_right=[3, 4],
        caption="表 2: 主熱流経路に沿った温度降下(# は上図の番号に対応)")

    # --- §3 熱抵抗の内訳 ---------------------------------------------
    rows = analysis["rows"]
    breakdown_chart = charts.bar_chart(
        [(r["label"], r["R"]) for r in rows], "K/W",
        highlight={r["label"] for r in rows if r["on_path"]})
    r_sum_path = sum(r["R"] for r in rows if r["on_path"]) or sum(
        r["R"] for r in rows)
    breakdown_rows = []
    for r in rows:
        share = 100 * r["R"] / r_sum_path if r["on_path"] and r_sum_path else None
        breakdown_rows.append([
            _esc(r["label"]),
            f'<span class="chip {"on-path" if r["on_path"] else "off-path"}">'
            f'{"主経路" if r["on_path"] else "分岐"}</span>'
            + ('<span class="chip nl">非線形</span>' if r["nonlinear"] else ""),
            f"{r['R']:.4f}", f"{r['dT']:.2f}", f"{r['Q']:.2f}",
            f"{share:.1f}" if share is not None else "—"])
    breakdown_table = _table(
        ["要素", "区分", "熱抵抗 [K/W]", "温度降下 [K]", "熱流 [W]", "主経路の割合 [%]"],
        breakdown_rows, align_right=[2, 3, 4, 5],
        caption="表 3: 要素ごとの熱抵抗(降順)")

    on_path_rows = [r for r in rows if r["on_path"]] or rows
    dominant = max(on_path_rows, key=lambda r: r["R"])

    # --- §4 感度分析 -------------------------------------------------
    sens = analysis["sensitivity"][:8] + analysis["conditions"]
    sens_chart = charts.diverging_bar(sens, "K")
    sens_table = _table(
        ["変更", f"{target} の温度変化 [K]"],
        [[_esc(name), f"{value:+.2f}"] for name, value in sens],
        align_right=[1], caption="表 4: 感度分析(基準からの温度変化)")

    # --- §5 過渡 -----------------------------------------------------
    transient_section = ""
    tr = analysis["transient"]
    if tr:
        shade = None
        shade_label = ""
        if tr["schedule"] and len(tr["schedule"]) >= 3:
            powers = [q for _, q in tr["schedule"]]
            peak_index = max(range(len(powers)), key=lambda i: powers[i])
            if 0 < peak_index < len(tr["schedule"]):
                start = float(tr["schedule"][peak_index][0])
                end = (float(tr["schedule"][peak_index + 1][0])
                       if peak_index + 1 < len(tr["schedule"])
                       else float(tr["time"][-1]))
                shade = (start, end)
                shade_label = f"{powers[peak_index]:.0f} W"
        line = charts.line_chart(tr["time"], tr["series"], shade=shade,
                                 shade_label=shade_label)
        peak_rows = [[_esc(name), f"{peak:.2f}", f"{when:.1f}",
                      f"{temps[name]:.2f}", f"{peak - temps[name]:+.2f}"]
                     for name, (peak, when) in tr["peaks"].items()]
        peak_table = _table(
            ["ノード", "ピーク温度 [degC]", "到達時刻 [s]", "定常温度 [degC]", "差 [K]"],
            peak_rows, align_right=[1, 2, 3, 4],
            caption="表 5: 過渡応答のピーク")
        transient_section = f"""
<section id="transient">
  <h2><span class="sec-no">5</span>過渡応答</h2>
  <p>熱容量を含めた時間応答。定常設計だけでは分からない「短時間なら許容できる
  負荷」の判断に使う。時間刻み {tr['opts'].get('dt')} s、
  終了時刻 {tr['opts'].get('t_end')} s、初期温度
  {tr['opts'].get('initial', '—')} degC。</p>
  <div class="figure">{line}</div>
  {peak_table}
</section>"""

    # --- §6 検算 -----------------------------------------------------
    checks = [
        ["エネルギー保存",
         f"発熱 {analysis['total_heat']:.4g} W / 境界の出入り合計 "
         f"{analysis['boundary_heat']:.4g} W",
         f"残差 {analysis['residual']:.2e} W",
         "good" if analysis["residual"] < 1e-6 else "critical"],
        ["非線形反復",
         ("温度依存の要素あり" if analysis["nonlinear"] else "全要素が線形"),
         (f"{analysis['iterations']} 回で収束(残差 "
          f"{analysis['sol'].residual:.1e} K)" if analysis["nonlinear"]
          else "反復不要(直接解)"),
         "good"],
        ["主熱流経路",
         _PATH_CHAIN,
         f"経路上の合計 {r_sum_path:.4f} K/W "
         f"(合成熱抵抗 {analysis['r_total']:.4f} K/W の "
         f"{100 * r_sum_path / analysis['r_total']:.0f} %)",
         "good" if r_sum_path / analysis["r_total"] > 0.9 else "warning"],
        ["支配抵抗",
         dominant["label"],
         f"{dominant['R']:.4f} K/W(経路の "
         f"{100 * dominant['R'] / r_sum_path:.0f} %)",
         "good"],
    ]
    auto_notes = [d for d in analysis["descriptions"].values() if "自動決定" in d]
    if auto_notes:
        checks.append(["拡がり抵抗の自己整合",
                       "h: auto を使用",
                       f"{len(auto_notes)} 要素が下流の合成抵抗と整合済み",
                       "good"])
    path_chain = ('<span class="path-chain">'
                  + " → ".join(_esc(n) for n in path) + "</span>")
    check_rows = [[_esc(name), (path_chain if what is _PATH_CHAIN else _esc(what)),
                   _esc(result),
                   f'<span class="chip status-{cls}">'
                   f'{"OK" if cls == "good" else "要確認"}</span>']
                  for name, what, result, cls in checks]
    checks_table = _table(["項目", "内容", "結果", "判定"], check_rows,
                          caption="表 6: 検算")

    # --- §7 前提と適用限界 -------------------------------------------
    assumption_html = "".join(
        f"<div class=\"assumption\"><h3>{_esc(name)}</h3><p>{_esc(text)}</p></div>"
        for name, text in assumptions(analysis))

    # --- 付録 --------------------------------------------------------
    element_rows = [[_esc(r["label"]), _esc(r["from"]), _esc(r["to"]),
                     _esc(r["note"])] for r in analysis["rows"]]
    element_table = _table(["要素", "from", "to", "抵抗の根拠"], element_rows,
                           caption="表 7: 各要素がどの物理量から作られたか")

    style = open(DEFAULT_STYLE, encoding="utf-8").read()
    limit_line = (f"設計上限 {limit:.0f} degC" if limit is not None
                  else "設計上限は未設定")

    return f"""<title>{_esc(title)}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans+JP:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>{style}</style>

<article class="doc">
<header class="masthead">
  <div class="eyebrow">1D 熱回路網 解析レポート</div>
  <h1>{_esc(title)}</h1>
  <div class="meta">
    <span>モデル: <code>{_esc(spec.get('name', 'model'))}</code></span>
    <span>生成: {generated}</span>
    <span>ソルバ: thermalnet(節点法)</span>
  </div>
  <div class="verdict verdict-{verdict_class}">
    <span class="verdict-mark" aria-hidden="true"></span>
    <strong>{_esc(verdict)}</strong>
    <span>{_esc(target)} = {t_target:.1f} degC / {_esc(limit_line)}</span>
  </div>
</header>

<div class="tiles">{tile_html}</div>

<nav class="toc" aria-label="目次">
  <a href="#conditions">1 解析条件</a>
  <a href="#temperature">2 温度分布</a>
  <a href="#breakdown">3 熱抵抗の内訳</a>
  <a href="#sensitivity">4 感度分析</a>
  {'<a href="#transient">5 過渡応答</a>' if tr else ''}
  <a href="#checks">{'6' if tr else '5'} 検算</a>
  <a href="#assumptions">{'7' if tr else '6'} 前提と適用限界</a>
</nav>

<section id="conditions">
  <h2><span class="sec-no">1</span>解析条件</h2>
  <p>ノード {len(net.nodes)} 個、要素 {len(net.branches)} 本の熱回路として解いた。
  {_esc(analysis['heat_label'])}は {analysis['reference_heat']:.4g} W、境界は
  {_esc('、'.join(f'{n} = {net.nodes[n].fixed:.1f} degC' for n in analysis['sinks']))}。</p>
  {conditions_table}
</section>

<section id="temperature">
  <h2><span class="sec-no">2</span>温度分布</h2>
  <p>発熱部から境界まで、熱流が最も大きい枝をたどった<strong>主熱流経路</strong>に
  沿って温度を並べたもの。段差の大きい区間がその設計のボトルネックにあたる。</p>
  <div class="figure">{ladder}</div>
  {path_table}
</section>

<section id="breakdown">
  <h2><span class="sec-no">3</span>熱抵抗の内訳</h2>
  <p>要素ごとの熱抵抗を大きい順に並べたもの。濃い棒が主熱流経路上の要素で、
  薄い棒は分岐(並列経路)。この設計では
  <strong>{_esc(dominant['label'])}</strong> が
  {dominant['R']:.4f} K/W と最大で、主経路の
  {100 * dominant['R'] / r_sum_path:.0f} % を占める。
  対策はここから着手するのが最も効率がよい。</p>
  <div class="figure">{breakdown_chart}</div>
  {breakdown_table}
</section>

<section id="sensitivity">
  <h2><span class="sec-no">4</span>感度分析</h2>
  <p>各要素の熱抵抗を 20 % 改善した場合と、運転条件を振った場合の
  {_esc(target)} の温度変化。青(左)が改善、赤(右)が悪化。
  抵抗が線形な区間では内訳表と等価な情報になるが、温度依存の要素があると
  順位が入れ替わるため、実際に解き直した値を載せている。</p>
  <div class="figure">{sens_chart}</div>
  {sens_table}
</section>
{transient_section}
<section id="checks">
  <h2><span class="sec-no">{'6' if tr else '5'}</span>検算</h2>
  <p>1D モデルは数字がいくらでも出てしまうので、結果を出す前に必ず確認する項目。</p>
  {checks_table}
</section>

<section id="assumptions">
  <h2><span class="sec-no">{'7' if tr else '6'}</span>前提と適用限界</h2>
  <p>このモデルをどこまで信じてよいか。ここを書かないレポートは使い回されて事故になる。</p>
  <div class="assumptions">{assumption_html}</div>
</section>

<section id="appendix">
  <h2 class="appendix-head">付録　モデル定義</h2>
  {element_table}
  <details>
    <summary>モデルファイル(YAML)全文</summary>
    <pre><code>{_esc(analysis.get('source_text', ''))}</code></pre>
  </details>
</section>

<footer class="colophon">
  <p>thermalnet の節点法ソルバで生成。数値は本レポートの前提条件のもとでのみ有効。
  形状・流量・熱伝達率を変更した場合はモデルの再校正が必要。</p>
</footer>
</article>

<div id="tip" role="status" aria-live="polite"></div>
{TOOLTIP_JS}
"""


# =====================================================================
# CLI
# =====================================================================
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="モデルファイルから熱設計レポート(HTML)を生成する")
    parser.add_argument("model", help="モデルファイル (.yaml / .json)")
    parser.add_argument("-o", "--output", help="出力する HTML のパス")
    parser.add_argument("--limit", type=float, help="設計上限温度 [degC]")
    parser.add_argument("--target", help="評価対象ノード(既定は最高温度のノード)")
    parser.add_argument("--title", help="レポートの表題")
    args = parser.parse_args(argv)

    spec = load_spec(args.model)
    analysis = analyze(spec, limit=args.limit, target=args.target)
    with open(args.model, encoding="utf-8") as fh:
        analysis["source_text"] = fh.read()

    title = args.title or f"{spec.get('name', 'モデル')} 熱設計レポート"
    output = args.output or os.path.splitext(args.model)[0] + "_report.html"
    with open(output, "w", encoding="utf-8") as fh:
        fh.write(render_html(analysis, title))
    print(f"レポートを書き出しました: {output}")
    print(f"  {analysis['target']} = {analysis['temps'][analysis['target']]:.2f} degC"
          f" / 合成熱抵抗 {analysis['r_total']:.4f} K/W")
    return 0


if __name__ == "__main__":
    sys.exit(main())
