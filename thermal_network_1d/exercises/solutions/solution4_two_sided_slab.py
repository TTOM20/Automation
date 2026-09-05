"""Q4 解答: 両面冷却の平板を 1D 離散化する"""

import _path  # noqa: F401
import numpy as np

from thermalnet import ThermalNetwork, biot_number

L, K, H, Q3, T_INF = 0.020, 0.8, 15.0, 3.0e4, 25.0
AREA = 1.0


def exact(x):
    """中心 x=0、両面 x=+-L/2 が対流。

    k T'' + q''' = 0, 対称性から T'(0) = 0
    表面熱流束 q''' L/2 = h (T_s - T_inf) より T_s = T_inf + q''' L / (2h)
    T(x) = T_s + q'''((L/2)^2 - x^2) / (2k)
    """
    t_surf = T_INF + Q3 * L / (2 * H)
    return t_surf + Q3 * ((L / 2) ** 2 - x ** 2) / (2 * K)


def build(n_cells):
    """両端が対流のはしご回路。両端セルに『半セル伝導 + 対流』が付く。"""
    dx = L / n_cells
    net = ThermalNetwork(f"slab{n_cells}")
    for i in range(n_cells):
        net.add_node(f"n{i}", heat=Q3 * AREA * dx)
    net.add_node("air", fixed=T_INF)
    for i in range(n_cells - 1):
        net.add_resistor(f"n{i}", f"n{i+1}", dx / (K * AREA))
    r_end = dx / (2 * K * AREA) + 1.0 / (H * AREA)
    net.add_resistor("n0", "air", r_end, "left")
    net.add_resistor(f"n{n_cells-1}", "air", r_end, "right")
    return net


print("Q4: 両面冷却平板の 1D 離散化\n")
t_max_exact = exact(0.0)
t_surf_exact = T_INF + Q3 * L / (2 * H)
parabola = Q3 * L ** 2 / (8 * K)
print(f"  解析解: 表面 {t_surf_exact:.2f} degC, 中心 {t_max_exact:.2f} degC "
      f"(内部の温度差 {parabola:.2f} K)")
print(f"  Bi = h (L/2) / k = {biot_number(H, K, L / 2):.4f}"
      "  <- 0.1 を超えているので 1 ノード近似は使えない\n")

print("     N     数値の最高温度 [degC]   誤差 [K]    前段との比")
prev = None
for n in (2, 4, 8, 16, 32, 64):
    sol = build(n).solve_steady(guess=T_INF)
    temps = np.array([sol[f"n{i}"] for i in range(n)])
    xc = (np.arange(n) + 0.5) * (L / n) - L / 2
    err = float(np.max(np.abs(temps - exact(xc))))
    ratio = "-" if prev is None else f"{prev / err:.2f}"
    print(f"    {n:3d}     {temps.max():16.3f}   {err:9.4f}   {ratio:>10s}")
    prev = err
    # エネルギー保存の検算(必ずやる)
    assert abs(sol.flow_of("left") + sol.flow_of("right") - Q3 * L) < 1e-9

sol1 = build(1).solve_steady(guess=T_INF)
print(f"""
  1 ノードで済ませる 2 通りのやり方と、その誤差:
    (a) 内部温度差を完全に無視 (T = T_surf)       {t_surf_exact:.2f} degC
        -> 真の最高温度 {t_max_exact:.2f} degC を {t_max_exact - t_surf_exact:.2f} K 過小評価
    (b) N=1 のはしご (半セル抵抗 L/2k を両側に)   {sol1['n0']:.2f} degC
        -> 逆に {sol1['n0'] - t_max_exact:.2f} K 過大評価

  ずれの大きさはどちらも {parabola:.2f} K、つまり内部の放物線分布そのもの。
  (b) が過大になるのは、N=1 だと「セル中心から表面までの半セル」が
  板厚の半分になり、伝導抵抗を実際より長く見積もるため。
  N を 2 以上にすれば両側から挟まれて正しくなる(上の表)。

  樹脂のように k が小さい材料では Bi がすぐ 0.1 を超える。
  「部品は 1 ノード」と決めつけず、Bi を計算してからノード数を決めること。
  分割数は 8 もあれば誤差 0.03 K 未満で足りる。
""")
