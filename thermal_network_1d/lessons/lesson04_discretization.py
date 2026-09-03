"""lesson04: CFD の離散化はそのまま熱回路である

学ぶこと
  1. 1D 伝導方程式を有限体積法で離散化すると、抵抗のはしご回路になること
  2. 境界セルは「半セル伝導 + 対流」の直列であること(ここを間違えやすい)
  3. 分割数を増やすと解析解に 2 次精度で収束すること(格子収束の確認方法)
  4. 何ノードに割るべきかはビオ数が教えてくれること

対象問題(解析解あり):
    厚さ L の平板、内部一様発熱 q''' [W/m^3]
    x=0 : 断熱、  x=L : 熱伝達率 h で T_inf へ

    解析解:  T(x) = T_inf + q''' L / h + q''' (L^2 - x^2) / (2 k)

実行:  python lessons/lesson04_discretization.py
"""

from _util import header, section, savefig, plt
import numpy as np

from thermalnet import ThermalNetwork, biot_number

header("lesson04  1D 離散化 = はしご型熱回路")

L = 0.05          # 板厚 [m]
k = 15.0          # 熱伝導率 [W/(m K)] (ステンレス相当)
h = 80.0          # 熱伝達率 [W/(m^2 K)]
q3 = 2.0e5        # 体積発熱 [W/m^3]
t_inf = 20.0      # 周囲温度 [degC]
area = 1.0        # 単位面積 [m^2] で考える

print(f"""
  L = {L} m, k = {k} W/(m K), h = {h} W/(m^2 K),
  q''' = {q3:.0e} W/m^3 (= {q3 * L * area:.0f} W/m^2 の熱流束), T_inf = {t_inf} degC
""")


def exact(x):
    return t_inf + q3 * L / h + q3 * (L ** 2 - x ** 2) / (2.0 * k)


def build_chain(n_cells: int) -> ThermalNetwork:
    """板厚方向を n_cells に割った「はしご回路」を作る。

        [n0]--R--[n1]--R-- ... --[nN-1]--(R/2 + R_conv)--[T_inf]
         |        |                |
        Q_i      Q_i              Q_i        Q_i = q''' A dx

        R = dx / (k A)                 : セル中心間の伝導抵抗
        R/2 = dx / (2 k A)             : 最終セル中心から表面までの半セル分
        R_conv = 1 / (h A)             : 表面から流体まで
    """
    dx = L / n_cells
    net = ThermalNetwork(f"slab_{n_cells}")
    for i in range(n_cells):
        net.add_node(f"n{i}", heat=q3 * area * dx)     # 各セルの発熱を電流源に
    net.add_node("fluid", fixed=t_inf)
    for i in range(n_cells - 1):
        net.add_resistor(f"n{i}", f"n{i + 1}", dx / (k * area), f"cond {i}-{i+1}")
    net.add_resistor(f"n{n_cells - 1}", "fluid",
                     dx / (2.0 * k * area) + 1.0 / (h * area), "half cell + conv")
    return net


# ---------------------------------------------------------------- 粗い格子
section("(1) まず 4 分割で解いてみる")

net4 = build_chain(4)
print(net4.describe())
sol4 = net4.solve_steady(guess=t_inf)
dx4 = L / 4
x4 = (np.arange(4) + 0.5) * dx4
print("\n   セル   x [mm]   数値解 [degC]   解析解 [degC]   誤差 [K]")
for i in range(4):
    tn = sol4[f"n{i}"]
    print(f"    n{i}   {1000 * x4[i]:6.2f}   {tn:11.3f}   {exact(x4[i]):11.3f}"
          f"   {tn - exact(x4[i]):+7.3f}")
print(f"\n  表面温度(数値): {t_inf + q3 * L * area / (h * area):.3f} degC"
      f"   <- 全体のエネルギー保存から一意に決まる")
print("  エネルギー保存の検算: 発熱 %.1f W / 流体へ %.1f W"
      % (q3 * L * area, sol4.flow_of("half cell + conv")))
print("""
  誤差が全セルで同じ +0.26 K になっている点に注目。放物線分布の「形」は
  離散式が厳密に再現しており、ずれは境界の半セル処理に由来する定数オフセット。
  誤差の出どころを切り分けられるのが、小さいモデルで勉強する利点。""")

# ---------------------------------------------------------------- 格子収束
section("(2) 格子収束: 分割数を倍にすると誤差はどうなるか")

counts = [2, 4, 8, 16, 32, 64, 128]
errors = []
for n_cells in counts:
    net = build_chain(n_cells)
    sol = net.solve_steady(guess=t_inf)
    dx = L / n_cells
    xc = (np.arange(n_cells) + 0.5) * dx
    num = np.array([sol[f"n{i}"] for i in range(n_cells)])
    errors.append(float(np.max(np.abs(num - exact(xc)))))

print("   分割数   最大誤差 [K]   前段との比   観測次数")
prev = None
for n_cells, e in zip(counts, errors):
    if prev is None:
        print(f"    {n_cells:4d}   {e:11.3e}        -           -")
    else:
        ratio = prev / e
        print(f"    {n_cells:4d}   {e:11.3e}   {ratio:8.3f}   {np.log2(ratio):8.3f}")
    prev = e
order = np.polyfit(np.log(L / np.array(counts)), np.log(errors), 1)[0]
print(f"\n  最小二乗フィットによる収束次数 p = {order:.3f}  (理論値 2)")
print("""
  これは CFD で必ずやる「格子収束性の確認」と同じ作業。
  1D 熱回路でも分割が粗ければ誤差は出る。ただし本質的な違いがある:
  * CFD の格子誤差は流れ場・境界層の解像度に起因し、収束させるのが高価
  * 1D 回路の誤差は「モデルの立て方」に起因し、格子は数十点で十分
  1D 化の勝負どころは格子ではなく、相関式と抵抗の切り出し方(lesson06-08)。
""")

# ---------------------------------------------------------------- 集中定数
section("(3) 1 ノードで済ませてよいか? --- ビオ数の実演")

bi = biot_number(h, k, L)
print(f"  Bi = h L / k = {bi:.3f}")
t_lumped = t_inf + q3 * L * area / (h * area)   # 板内温度差を無視した場合
t_exact_max = exact(0.0)
print(f"  1 ノード近似の最高温度 {t_lumped:.2f} degC")
print(f"  真の最高温度           {t_exact_max:.2f} degC"
      f"   -> {t_exact_max - t_lumped:.2f} K の過小評価")

print("\n  Bi を変えて、1 ノード近似の誤差がどう変わるかを見る:")
print("     h [W/m^2K]     Bi      1ノード誤差 [K]   相対誤差(温度上昇比)")
for h_test in [5.0, 20.0, 80.0, 300.0, 1000.0]:
    bi_t = biot_number(h_test, k, L)
    t_l = t_inf + q3 * L / h_test
    t_e = t_inf + q3 * L / h_test + q3 * L ** 2 / (2 * k)
    rise = t_e - t_inf
    print(f"     {h_test:8.0f}   {bi_t:6.3f}   {t_e - t_l:12.2f}"
          f"        {100 * (t_e - t_l) / rise:6.1f} %")
print("""  Bi < 0.1 の行では誤差が数 % に収まる。これが集中定数近似の適用限界。
  Bi が大きいほど「物体を何ノードにも割る」必要がある = 1D 離散化が要る。""")

# ---------------------------------------------------------------- 図
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

xf = np.linspace(0, L, 200)
axes[0].plot(xf * 1000, exact(xf), "k-", lw=2, label="analytical")
for n_cells, style in [(2, "s--"), (4, "o--"), (16, "^--")]:
    net = build_chain(n_cells)
    sol = net.solve_steady(guess=t_inf)
    dx = L / n_cells
    xc = (np.arange(n_cells) + 0.5) * dx
    axes[0].plot(xc * 1000, [sol[f"n{i}"] for i in range(n_cells)], style,
                 label=f"{n_cells} nodes")
axes[0].set_xlabel("x [mm]   (x=0 adiabatic, x=L convective)")
axes[0].set_ylabel("temperature [degC]")
axes[0].set_title("1D slab with internal heat generation")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].loglog(L / np.array(counts), errors, "o-", lw=2, label="ladder network")
ref = errors[0] * (np.array(counts, dtype=float) / counts[0]) ** -2.0
axes[1].loglog(L / np.array(counts), ref, "k--", label="2nd order slope")
axes[1].set_xlabel("cell size dx [m]")
axes[1].set_ylabel("max error [K]")
axes[1].set_title(f"Grid convergence (observed p = {order:.2f})")
axes[1].legend()
axes[1].grid(alpha=0.3, which="both")

savefig(fig, "lesson04_discretization.png")

print("""
まとめ
  * 有限体積の離散式 = 節点方程式。CFD と熱回路は同じ骨格を持つ。
  * 境界セルは「半セル伝導 + 表面抵抗」の直列。ここを忘れると系統誤差になる。
  * 格子収束は 1D でも確認する癖をつける(観測次数 ~2)。
  * ノード数の必要量はビオ数で見積もれる。
次: lesson05  過渡応答(熱容量・時定数・陰解法)
""")
