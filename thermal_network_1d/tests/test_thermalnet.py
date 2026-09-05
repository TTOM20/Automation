"""thermalnet の検証テスト --- 解析解・保存則との突き合わせ。

実行方法:
    python tests/test_thermalnet.py        (pytest なしで動く)
    pytest tests/test_thermalnet.py        (pytest があればこちらでも)

教材として「数値解をどう検算するか」の見本も兼ねている。
新しい素子や相関式を足したら、必ずここに検算を 1 本足すこと。
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thermalnet import (ThermalNetwork, biot_number, capacitance, fin_efficiency,
                        h_radiation, r_cond_cylinder, r_cond_plane, r_conv, r_flow,
                        r_spreading_disk, r_spreading_halfspace)
from thermalnet.elements import SIGMA, T_ABS
from thermalnet.fvm import axisymmetric_spreader


# ------------------------------------------------------------------ 素子
def test_conduction_matches_formula():
    assert math.isclose(r_cond_plane(0.01, 50.0, 0.002), 0.01 / (50.0 * 0.002))
    # 円筒: 薄肉極限では平板に一致する
    r_in, thick, k, length = 0.05, 1e-5, 20.0, 1.0
    cyl = r_cond_cylinder(r_in, r_in + thick, k, length)
    plane = r_cond_plane(thick, k, 2 * math.pi * r_in * length)
    assert abs(cyl - plane) / plane < 1e-4


def test_radiation_linearization_is_exact_at_the_operating_point():
    """h_r による線形化は、その温度で厳密に Stefan-Boltzmann と一致する。"""
    eps, ts, tsur, area = 0.8, 90.0, 25.0, 0.05
    exact = eps * SIGMA * area * ((ts + T_ABS) ** 4 - (tsur + T_ABS) ** 4)
    linear = h_radiation(eps, ts, tsur) * area * (ts - tsur)
    assert abs(exact - linear) / exact < 1e-12


def test_fin_efficiency_limits():
    # 熱伝導率が非常に大きい -> 効率 1
    assert fin_efficiency(50.0, 1e7, 0.002, 0.02) > 0.999
    # 長く薄いフィン -> 効率低下
    assert fin_efficiency(200.0, 15.0, 0.0005, 0.10) < 0.2


def test_spreading_halfspace_formulas():
    a, k = 0.004, 200.0
    assert math.isclose(r_spreading_halfspace(a, k, True), 1.0 / (4 * k * a))
    assert math.isclose(r_spreading_halfspace(a, k, False),
                        8.0 / (3 * math.pi ** 2 * k * a))
    # 等熱流束のほうが大きい(平均温度が高くなるため)
    assert r_spreading_halfspace(a, k, False) > r_spreading_halfspace(a, k, True)


def test_caloric_resistance():
    mdot, cp, q = 0.005, 1007.0, 50.0
    dt = q * r_flow(mdot, cp)
    assert math.isclose(dt, q / (mdot * cp))


# ------------------------------------------------------------------ 回路
def test_series_and_parallel():
    net = ThermalNetwork()
    net.add_node("hot", heat=10.0)
    net.add_node("mid")
    net.add_node("cold", fixed=0.0)
    net.add_resistor("hot", "mid", 2.0, "s1")
    net.add_resistor("mid", "cold", 3.0, "p1")
    net.add_resistor("mid", "cold", 6.0, "p2")
    sol = net.solve_steady()
    r_expected = 2.0 + 1.0 / (1.0 / 3.0 + 1.0 / 6.0)   # = 4.0
    assert math.isclose(sol["hot"], 10.0 * r_expected, rel_tol=1e-12)
    assert math.isclose(net.total_resistance("hot", "cold"), r_expected,
                        rel_tol=1e-12)
    # 並列枝の分流比は抵抗の逆比
    assert math.isclose(sol.flow_of("p1") / sol.flow_of("p2"), 2.0, rel_tol=1e-12)


def test_energy_is_conserved():
    net = ThermalNetwork()
    net.add_node("a", heat=7.0)
    net.add_node("b", heat=3.0)
    net.add_node("amb", fixed=20.0)
    net.add_resistor("a", "b", 1.5)
    net.add_resistor("a", "amb", 4.0, "a_out")
    net.add_resistor("b", "amb", 2.0, "b_out")
    sol = net.solve_steady()
    out = sol.flow_of("a_out") + sol.flow_of("b_out")
    assert abs(out - 10.0) < 1e-10


def test_ambiguous_parallel_branch_raises():
    net = ThermalNetwork()
    net.add_node("x", heat=1.0)
    net.add_node("y", fixed=0.0)
    net.add_resistor("x", "y", 1.0, "conv")
    net.add_resistor("x", "y", 2.0, "rad")
    sol = net.solve_steady()
    try:
        sol.flow("x", "y")
    except KeyError as exc:
        assert "並列" in str(exc)
    else:
        raise AssertionError("並列枝の曖昧さを検出できていない")
    assert math.isclose(sol.flow("x", "y", label="conv"), 1.0 / 1.5, rel_tol=1e-9)


def test_total_resistance_restores_network_state():
    net = ThermalNetwork()
    net.add_node("j", heat=5.0)
    net.add_node("amb", fixed=25.0)
    net.add_resistor("j", "amb", 2.0)
    before = net.solve_steady()["j"]
    net.total_resistance("j", "amb")
    after = net.solve_steady()["j"]
    assert math.isclose(before, after)


# ------------------------------------------------------------------ 離散化
def _slab_chain(n_cells, L, k, h, q3, t_inf, area=1.0):
    dx = L / n_cells
    net = ThermalNetwork()
    for i in range(n_cells):
        net.add_node(f"n{i}", heat=q3 * area * dx,
                     capacitance=0.0)
    net.add_node("fluid", fixed=t_inf)
    for i in range(n_cells - 1):
        net.add_resistor(f"n{i}", f"n{i+1}", dx / (k * area))
    net.add_resistor(f"n{n_cells-1}", "fluid", dx / (2 * k * area) + 1 / (h * area),
                     "out")
    return net


def test_1d_slab_converges_second_order():
    L, k, h, q3, t_inf = 0.05, 15.0, 80.0, 2.0e5, 20.0

    def exact(x):
        return t_inf + q3 * L / h + q3 * (L ** 2 - x ** 2) / (2 * k)

    errs = []
    counts = [8, 16, 32, 64]
    for n in counts:
        sol = _slab_chain(n, L, k, h, q3, t_inf).solve_steady(guess=t_inf)
        xc = (np.arange(n) + 0.5) * L / n
        num = np.array([sol[f"n{i}"] for i in range(n)])
        errs.append(np.max(np.abs(num - exact(xc))))
    order = np.polyfit(np.log(L / np.array(counts)), np.log(errs), 1)[0]
    assert 1.9 < order < 2.1, f"収束次数が 2 でない: {order}"
    # 表面熱流束は格子によらず厳密(全体のエネルギー保存)
    sol = _slab_chain(8, L, k, h, q3, t_inf).solve_steady(guess=t_inf)
    assert abs(sol.flow_of("out") - q3 * L) < 1e-9


def test_biot_number():
    assert math.isclose(biot_number(80.0, 15.0, 0.05), 80.0 * 0.05 / 15.0)


# ------------------------------------------------------------------ 過渡
def test_lumped_transient_matches_exponential():
    cap = capacitance(2700.0, 900.0, 0.03 ** 3)
    r = r_conv(25.0, 6 * 0.03 ** 2)
    tau = r * cap
    net = ThermalNetwork()
    net.add_node("block", capacitance=cap)
    net.add_node("air", fixed=25.0)
    net.add_resistor("block", "air", r)
    tr = net.solve_transient(3 * tau, tau / 500, initial={"block": 200.0}, theta=0.5)
    analytic = 25.0 + 175.0 * np.exp(-tr.time / tau)
    assert np.max(np.abs(tr["block"] - analytic)) < 1e-3


def test_crank_nicolson_beats_backward_euler():
    cap, r = 100.0, 2.0
    tau = cap * r

    def run(theta, dt):
        net = ThermalNetwork()
        net.add_node("b", capacitance=cap)
        net.add_node("a", fixed=0.0)
        net.add_resistor("b", "a", r)
        tr = net.solve_transient(2 * tau, dt, initial={"b": 100.0}, theta=theta)
        return np.max(np.abs(tr["b"] - 100.0 * np.exp(-tr.time / tau)))

    assert run(0.5, tau / 20) < run(1.0, tau / 20) / 10


def test_explicit_scheme_is_unstable_beyond_limit():
    cap, r = 100.0, 2.0
    dt_limit = 2 * r * cap
    net = ThermalNetwork()
    net.add_node("b", capacitance=cap)
    net.add_node("a", fixed=0.0)
    net.add_resistor("b", "a", r)
    stable = net.solve_transient(20 * dt_limit, 0.5 * dt_limit,
                                 initial={"b": 100.0}, theta=0.0)
    unstable = net.solve_transient(20 * dt_limit, 1.3 * dt_limit,
                                   initial={"b": 100.0}, theta=0.0)
    assert np.max(np.abs(stable["b"])) <= 100.0 + 1e-9
    assert np.max(np.abs(unstable["b"])) > 1e3


def test_zero_capacitance_node_is_algebraic():
    """熱容量 0 のノードがあっても過渡計算が破綻せず、定常解に収束する。"""
    net = ThermalNetwork()
    net.add_node("die", heat=10.0)                    # C = 0
    net.add_node("sink", capacitance=200.0)
    net.add_node("air", fixed=25.0)
    net.add_resistor("die", "sink", 0.5)
    net.add_resistor("sink", "air", 1.0)
    tr = net.solve_transient(3000.0, 1.0, initial=25.0)
    assert abs(tr.final()["die"] - (25.0 + 10.0 * 1.5)) < 1e-3


# ------------------------------------------------------------------ 非線形
def test_nonlinear_radiation_network():
    """輻射だけで放熱する板。解析的な釣り合いと一致するか。"""
    eps, area, t_amb, q = 0.9, 0.02, 25.0, 8.0
    net = ThermalNetwork()
    net.add_node("plate", heat=q)
    net.add_node("amb", fixed=t_amb)
    net.add_resistor("plate", "amb",
                     lambda ts, ta: r_conv(h_radiation(eps, ts, ta), area))
    sol = net.solve_steady(guess=t_amb + 50.0, tol=1e-10)
    ts = sol["plate"]
    residual = eps * SIGMA * area * ((ts + T_ABS) ** 4 - (t_amb + T_ABS) ** 4) - q
    assert abs(residual) < 1e-8, f"エネルギー残差 {residual}"


# ------------------------------------------------------------------ FVM
def test_fvm_energy_balance_and_decomposition():
    a, b, t, k, h, q = 0.005, 0.02, 0.003, 200.0, 500.0, 30.0
    res = axisymmetric_spreader(a, b, t, k, h, q, nr=120, nz=60, t_inf=25.0)
    # 定義どおり R_total = R_1d + R_spread
    assert math.isclose(res["R_total"], res["R_1d"] + res["R_spread"], rel_tol=1e-12)
    # 素朴な 1D より必ず高温になる(拡がり抵抗は正)
    assert res["R_spread"] > 0
    # 下面の平均温度は 1D 的な釣り合いと一致するはず: Q = h A (T_base - T_inf)
    t_base_expected = 25.0 + q / (h * math.pi * b ** 2)
    assert abs(res["T_base_mean"] - t_base_expected) < 0.05 * (t_base_expected - 25.0)


def test_fvm_grid_convergence():
    args = dict(a=0.005, b=0.02, t=0.003, k=200.0, h=500.0, q_total=1.0)
    coarse = axisymmetric_spreader(**args, nr=80, nz=40)["R_total"]
    fine = axisymmetric_spreader(**args, nr=320, nz=160)["R_total"]
    assert abs(fine - coarse) / fine < 1e-3, "格子依存性が大きすぎる"


def test_spreading_correlation_within_10_percent_of_fvm():
    """Lee et al. の相関式は FVM 解と 10 % 以内で一致する(それ以上は期待しない)。"""
    for a, b, t, k, h in [(0.005, 0.020, 0.003, 200.0, 500.0),
                          (0.010, 0.030, 0.005, 150.0, 300.0),
                          (0.008, 0.020, 0.010, 50.0, 200.0)]:
        num = axisymmetric_spreader(a, b, t, k, h, 1.0, 200, 100)["R_spread"]
        corr = r_spreading_disk(a, b, t, k, h)
        rel = abs(corr - num) / num
        assert rel < 0.10, f"相関式のずれ {100 * rel:.1f} % (a={a}, b={b}, t={t})"


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
