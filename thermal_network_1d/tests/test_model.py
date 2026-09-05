"""モデル定義層(thermalnet.model)の検証テスト。

実行:
    python tests/test_model.py
    pytest tests/test_model.py
"""

import copy
import glob
import io
import math
import os
import sys
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thermalnet import ThermalNetwork, elements as el, materials
from thermalnet.fluids import h_natural_vertical_plate
from thermalnet.model import (apply_override, build_network, load_spec, main,
                              report_steady)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_DIR = os.path.join(ROOT, "models")


# ------------------------------------------------------------------ 構築
def test_conduction_matches_hand_built_network():
    """モデル層で組んだ回路が、手で組んだ回路と厳密に一致すること。"""
    spec = {
        "nodes": {"hot": {"heat": 10.0}, "cold": {"fixed": 20.0}},
        "elements": [
            {"name": "wall", "type": "conduction", "from": "hot", "to": "mid",
             "k": 50.0, "thickness": 0.01, "area": 0.002},
            {"name": "film", "type": "convection", "from": "mid", "to": "cold",
             "h": 25.0, "area": 0.002},
        ],
    }
    net, _ = build_network(spec)
    got = net.solve_steady(guess=20.0)["hot"]

    hand = ThermalNetwork()
    hand.add_node("hot", heat=10.0)
    hand.add_node("mid")
    hand.add_node("cold", fixed=20.0)
    hand.add_resistor("hot", "mid", el.r_cond_plane(0.01, 50.0, 0.002))
    hand.add_resistor("mid", "cold", el.r_conv(25.0, 0.002))
    assert math.isclose(got, hand.solve_steady(guess=20.0)["hot"], rel_tol=1e-12)


def test_area_shorthands_are_equivalent():
    """area / size / width+length / diameter が同じ抵抗を与えること。"""
    def resistance(area_spec):
        spec = {
            "nodes": {"a": {"heat": 1.0}, "b": {"fixed": 0.0}},
            "elements": [dict({"name": "e", "type": "conduction", "from": "a",
                               "to": "b", "k": 10.0, "thickness": 0.001},
                              **area_spec)],
        }
        net, _ = build_network(spec)
        return net.branches[0].resistance

    base = resistance({"area": 0.0004})
    assert math.isclose(resistance({"size": 0.02}), base, rel_tol=1e-12)
    assert math.isclose(resistance({"width": 0.02, "length": 0.02}), base,
                        rel_tol=1e-12)
    circle = resistance({"diameter": 0.02})
    assert math.isclose(circle, 0.001 / (10.0 * math.pi * 0.01 ** 2), rel_tol=1e-12)


def test_material_lookup_and_local_override():
    spec = {
        "materials": {"my_alloy": {"k": 100.0, "rho": 3000.0, "cp": 800.0}},
        "nodes": {"a": {"heat": 1.0, "material": "my_alloy", "volume": 1e-4},
                  "b": {"fixed": 0.0}},
        "elements": [{"name": "e", "type": "conduction", "from": "a", "to": "b",
                      "material": "my_alloy", "thickness": 0.01, "area": 0.01}],
    }
    net, _ = build_network(spec)
    assert math.isclose(net.branches[0].resistance, 0.01 / (100.0 * 0.01))
    assert math.isclose(net.nodes["a"].capacitance, 3000.0 * 800.0 * 1e-4)
    # 未登録の材料はエラー(候補つき)
    try:
        materials.get("unobtainium")
    except KeyError as exc:
        assert "未登録" in str(exc)
    else:
        raise AssertionError("未登録材料を検出できていない")


def test_capacitance_from_mass_and_volume_agree():
    props = materials.get("aluminum")
    volume = 1e-4
    from_volume = {"material": "aluminum", "volume": volume}
    from_mass = {"material": "aluminum", "mass": props["rho"] * volume}
    caps = []
    for node_spec in (from_volume, from_mass):
        spec = {"nodes": {"a": dict(node_spec, heat=1.0), "b": {"fixed": 0.0}},
                "elements": [{"type": "resistance", "from": "a", "to": "b",
                              "value": 1.0}]}
        net, _ = build_network(spec)
        caps.append(net.nodes["a"].capacitance)
    assert math.isclose(caps[0], caps[1], rel_tol=1e-12)


# ------------------------------------------------------------------ 検証機能
def test_missing_boundary_is_rejected():
    spec = {"nodes": {"a": {"heat": 5.0}, "b": {}},
            "elements": [{"type": "resistance", "from": "a", "to": "b",
                          "value": 1.0}]}
    try:
        build_network(spec)
    except ValueError as exc:
        assert "温度を固定した" in str(exc)
    else:
        raise AssertionError("境界なしを検出できていない")


def test_floating_node_is_rejected():
    spec = {"nodes": {"a": {"heat": 5.0}, "amb": {"fixed": 25.0}, "orphan": {}},
            "elements": [{"type": "resistance", "from": "a", "to": "amb",
                          "value": 1.0}]}
    try:
        build_network(spec)
    except ValueError as exc:
        assert "orphan" in str(exc)
    else:
        raise AssertionError("浮きノードを検出できていない")


def test_unknown_element_type_is_rejected():
    spec = {"nodes": {"a": {"heat": 1.0}, "b": {"fixed": 0.0}},
            "elements": [{"type": "teleportation", "from": "a", "to": "b"}]}
    try:
        build_network(spec)
    except ValueError as exc:
        assert "未知の要素タイプ" in str(exc)
    else:
        raise AssertionError("未知の type を検出できていない")


# ------------------------------------------------------------------ 上書き
def test_apply_override_on_nodes_and_elements():
    spec = {"nodes": {"a": {"heat": 1.0}, "b": {"fixed": 0.0}},
            "elements": [{"name": "path", "type": "resistance", "from": "a",
                          "to": "b", "value": 2.0}]}
    apply_override(spec, "nodes.a.heat=10")
    apply_override(spec, "elements.path.value=0.5")
    net, _ = build_network(spec)
    assert math.isclose(net.solve_steady()["a"], 10.0 * 0.5, rel_tol=1e-12)


def test_disabled_element_is_skipped():
    spec = {"nodes": {"a": {"heat": 1.0}, "b": {"fixed": 0.0}},
            "elements": [
                {"name": "main", "type": "resistance", "from": "a", "to": "b",
                 "value": 2.0},
                {"name": "bypass", "type": "resistance", "from": "a", "to": "b",
                 "value": 2.0}]}
    both, _ = build_network(copy.deepcopy(spec))
    assert math.isclose(both.solve_steady()["a"], 1.0, rel_tol=1e-12)  # 並列で 1.0
    apply_override(spec, "elements.bypass.disabled=true")
    single, _ = build_network(spec)
    assert math.isclose(single.solve_steady()["a"], 2.0, rel_tol=1e-12)


def test_override_rejects_unknown_path():
    spec = {"nodes": {"a": {"heat": 1.0}, "b": {"fixed": 0.0}}, "elements": []}
    try:
        apply_override(spec, "nodes.nonexistent.heat=5")
    except KeyError:
        pass
    else:
        raise AssertionError("存在しないパスを検出できていない")


# ------------------------------------------------------------------ 要素
def test_flow_element_mean_factor():
    """mean: true は 1/(2 mdot cp)、false は 1/(mdot cp)。"""
    def build(mean):
        spec = {"nodes": {"a": {"heat": 1.0}, "b": {"fixed": 0.0}},
                "elements": [{"name": "f", "type": "flow", "from": "a", "to": "b",
                              "fluid": "water", "mass_flow": 0.05, "mean": mean}]}
        net, _ = build_network(spec)
        return net.branches[0].resistance

    cp = materials.get("water")["cp"]
    assert math.isclose(build(False), 1.0 / (0.05 * cp), rel_tol=1e-12)
    assert math.isclose(build(True), 0.5 / (0.05 * cp), rel_tol=1e-12)


def test_flow_units_lpm_and_cfm():
    """L/min と CFM が質量流量に正しく換算されること。"""
    rho = materials.get("water")["rho"]
    spec = {"nodes": {"a": {"heat": 1.0}, "b": {"fixed": 0.0}},
            "elements": [{"name": "f", "type": "flow", "from": "a", "to": "b",
                          "fluid": "water", "lpm": 6.0, "mean": False}]}
    net, _ = build_network(spec)
    expected_mdot = rho * 6.0 / 60000.0
    assert math.isclose(net.branches[0].resistance,
                        1.0 / (expected_mdot * materials.get("water")["cp"]),
                        rel_tol=1e-12)
    # 1 CFM = 0.000471947 m^3/s
    spec["elements"][0] = {"name": "f", "type": "flow", "from": "a", "to": "b",
                           "fluid": "air", "cfm": 10.0, "mean": False}
    net, _ = build_network(spec)
    air = materials.get("air")
    expected = 1.0 / (air["rho"] * 10.0 * 0.000471947 * air["cp"])
    assert math.isclose(net.branches[0].resistance, expected, rel_tol=1e-9)


def test_correlation_element_is_nonlinear_and_matches_direct_call():
    """correlation を使うと非線形になり、相関式を直接呼んだ値と一致すること。"""
    spec = {"nodes": {"case": {"heat": 20.0}, "amb": {"fixed": 25.0}},
            "elements": [{"name": "conv", "type": "convection", "from": "case",
                          "to": "amb", "correlation": "natural_vertical",
                          "height": 0.05, "area": 0.035}]}
    net, _ = build_network(spec)
    assert callable(net.branches[0].resistance)
    sol = net.solve_steady(guess=45.0, tol=1e-10)
    ts = sol["case"]
    h = h_natural_vertical_plate(ts, 25.0, 0.05)
    assert math.isclose(20.0, h * 0.035 * (ts - 25.0), rel_tol=1e-6)


def test_radiation_element_matches_stefan_boltzmann():
    spec = {"nodes": {"s": {"heat": 8.0}, "amb": {"fixed": 25.0}},
            "elements": [{"name": "rad", "type": "radiation", "from": "s",
                          "to": "amb", "emissivity": 0.9, "area": 0.02}]}
    net, _ = build_network(spec)
    ts = net.solve_steady(guess=80.0, tol=1e-11)["s"]
    flux = 0.9 * el.SIGMA * 0.02 * ((ts + el.T_ABS) ** 4 - (25.0 + el.T_ABS) ** 4)
    assert abs(flux - 8.0) < 1e-7


def test_interface_two_ways_agree():
    """resistivity 指定と k+thickness 指定が整合すること。"""
    area, thickness, k = 0.001, 5e-5, 3.0
    by_k = {"name": "t", "type": "interface", "from": "a", "to": "b",
            "k": k, "thickness": thickness, "area": area}
    by_rho = {"name": "t", "type": "interface", "from": "a", "to": "b",
              "resistivity": thickness / k, "area": area}
    values = []
    for element in (by_k, by_rho):
        spec = {"nodes": {"a": {"heat": 1.0}, "b": {"fixed": 0.0}},
                "elements": [element]}
        net, _ = build_network(spec)
        values.append(net.branches[0].resistance)
    assert math.isclose(values[0], values[1], rel_tol=1e-12)


def test_auto_spreading_is_self_consistent():
    """h: auto が h_eq = 1/(R_downstream x A_plate) を満たすこと。"""
    plate_area = 4.0e-3
    spec = {
        "nodes": {"src": {"heat": 100.0}, "mid": {}, "sink": {"fixed": 40.0}},
        "elements": [
            {"name": "spread", "type": "spreading", "from": "src", "to": "mid",
             "k": 170.0, "source_area": 1.0e-3, "plate_area": plate_area,
             "thickness": 0.00064, "h": "auto"},
            {"name": "down", "type": "resistance", "from": "mid", "to": "sink",
             "value": 0.025},
        ],
    }
    net, desc = build_network(spec)
    h_expected = 1.0 / (0.025 * plate_area)
    r_expected = el.r_spreading_disk(math.sqrt(1.0e-3 / math.pi),
                                     math.sqrt(plate_area / math.pi),
                                     0.00064, 170.0, h_expected)
    assert math.isclose(net.branches[0].resistance, r_expected, rel_tol=1e-9)
    assert "自動決定" in desc["spread"]


def test_heat_schedule_is_piecewise_constant():
    spec = {"nodes": {"a": {"heat_schedule": [[0, 10.0], [5, 50.0], [8, 10.0]],
                            "capacitance": 100.0},
                      "b": {"fixed": 0.0}},
            "elements": [{"type": "resistance", "from": "a", "to": "b",
                          "value": 1.0}]}
    net, _ = build_network(spec)
    q = net.nodes["a"].heat
    assert callable(q)
    assert q(0.0) == 10.0 and q(4.9) == 10.0
    assert q(5.0) == 50.0 and q(7.9) == 50.0
    assert q(8.0) == 10.0 and q(100.0) == 10.0


# ------------------------------------------------------------------ 例題モデル
def test_bundled_models_solve_and_conserve_energy():
    """models/ の全モデルが解けて、エネルギー保存を満たすこと。"""
    paths = sorted(glob.glob(os.path.join(MODEL_DIR, "*.yaml")))
    assert paths, "models/*.yaml が見つかりません"
    for path in paths:
        spec = load_spec(path)
        net, desc = build_network(spec)
        sol = net.solve_steady(guess=(spec.get("solve") or {}).get("guess", 25.0))
        supplied = sum((n.heat(0.0) if callable(n.heat) else n.heat)
                       for n in net.nodes.values())
        removed = 0.0
        for name, node in net.nodes.items():
            if node.fixed is None:
                continue
            for br in net.branches:
                r = br.value(sol[br.a], sol[br.b])
                if br.a == name:
                    removed -= (sol[br.a] - sol[br.b]) / r
                elif br.b == name:
                    removed += (sol[br.a] - sol[br.b]) / r
        assert abs(supplied - removed) < 1e-6, f"{path}: エネルギーが保存していない"
        assert report_steady(net, desc, sol)      # レポート生成が例外なく通る


def test_cpu_model_reproduces_lesson08():
    """cpu_heatsink.yaml が lesson08 の合成熱抵抗を再現すること。"""
    spec = load_spec(os.path.join(MODEL_DIR, "cpu_heatsink.yaml"))
    net, _ = build_network(spec)
    sol = net.solve_steady(guess=60.0)
    r_ja = (sol["die"] - 30.0) / 95.0
    assert abs(r_ja - 0.6403) < 0.005, f"R_ja = {r_ja}"


def test_cli_runs_on_every_model():
    """CLI が全モデルで正常終了すること(定常・スイープ・過渡)。"""
    for path in sorted(glob.glob(os.path.join(MODEL_DIR, "*.yaml"))):
        with redirect_stdout(io.StringIO()):
            assert main([path, "--quiet"]) == 0
    cpu = os.path.join(MODEL_DIR, "cpu_heatsink.yaml")
    with redirect_stdout(io.StringIO()):
        assert main([cpu, "--set", "nodes.die.heat=120", "--quiet"]) == 0
        assert main([cpu, "--sweep", "nodes.die.heat=50:150:3",
                     "--watch", "die"]) == 0
        assert main([cpu, "--transient", "--t-end", "10", "--dt", "0.1",
                     "--watch", "die"]) == 0


# ------------------------------------------------------------------ runner
def _main():
    tests = [(name, obj) for name, obj in sorted(globals().items())
             if name.startswith("test_") and callable(obj)]
    failures = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception as exc:                       # noqa: BLE001
            failures += 1
            print(f"  FAIL  {name}: {exc}")
    print(f"\n{len(tests) - failures} / {len(tests)} 件が成功")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_main())
