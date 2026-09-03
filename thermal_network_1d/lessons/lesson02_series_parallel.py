"""lesson02: 直列・並列合成 --- 多層壁と熱橋(サーマルブリッジ)

学ぶこと
  1. 直列 R = sum R_i、並列 1/R = sum 1/R_i(電気とまったく同じ)
  2. 多層壁の熱通過率 U = 1/(R_total A) と温度分布の求め方
  3. 「熱橋」を並列枝でモデル化すると、断熱の効きが数字で見える
  4. 手計算とソルバの答えが一致することを自分で確認する

実行:  python lessons/lesson02_series_parallel.py
"""

from _util import header, section, savefig, plt
import numpy as np

from thermalnet import ThermalNetwork, r_cond_plane, r_conv

header("lesson02  直列・並列合成と多層壁")

# ---------------------------------------------------------------- 多層壁
section("(1) 多層壁: 室内 -> 石膏ボード -> グラスウール -> 合板 -> 外気")

area = 10.0            # 壁面積 [m^2]
t_in, t_out = 20.0, -5.0
layers = [
    # (名前, 厚さ[m], 熱伝導率[W/(m K)])
    ("gypsum board", 0.0125, 0.22),
    ("glass wool", 0.100, 0.040),
    ("plywood", 0.012, 0.15),
]
h_in, h_out = 9.0, 23.0     # 室内側/屋外側の総合熱伝達率 [W/(m^2 K)]

r_list = [("inside film", r_conv(h_in, area))]
r_list += [(name, r_cond_plane(t, k, area)) for name, t, k in layers]
r_list += [("outside film", r_conv(h_out, area))]

r_total = sum(r for _, r in r_list)
q_wall = (t_in - t_out) / r_total
u_value = 1.0 / (r_total * area)

print("  層ごとの熱抵抗 [K/W] と温度降下 [K]")
for name, r in r_list:
    print(f"    {name:<16s} R = {r:8.5f}   dT = {q_wall * r:6.2f}"
          f"   ({100 * r / r_total:5.1f} %)")
print(f"\n  合計 R = {r_total:.5f} K/W   熱流 Q = {q_wall:.1f} W")
print(f"  熱通過率 U = 1/(R A) = {u_value:.3f} W/(m^2 K)"
      f"   (建築で言う U 値。1/U = R 値 x 面積)")

# 同じものを回路で解いて温度分布を得る
net = ThermalNetwork("wall")
net.add_node("room", fixed=t_in)
net.add_node("outdoor", fixed=t_out)
path = ["room", "surf_in", "n1", "n2", "surf_out", "outdoor"]
for n in path[1:-1]:
    net.add_node(n)
net.add_series(path, [r for _, r in r_list], [name for name, _ in r_list])
sol = net.solve_steady(guess=0.0)
print("\n  回路ソルバによる各界面温度 [degC]")
for n in path:
    print(f"    {n:<10s} {sol[n]:7.2f}")
print(f"  手計算の熱流 {q_wall:.3f} W  vs  ソルバ {sol.flow('room', 'surf_in'):.3f} W")
assert abs(q_wall - sol.flow("room", "surf_in")) < 1e-9

print("""
  注目: 断熱材(グラスウール)が全抵抗の大半を占め、温度もそこで大きく落ちる。
  「温度降下の大きい層 = 支配抵抗」であり、そこを改善しないと効果は出ない。
  これは電子機器でも建築でも自動車でも同じ、熱回路の最も実用的な使い方。
""")

# ---------------------------------------------------------------- 熱橋
section("(2) 熱橋: 断熱層を貫く木製スタッド(並列回路)")

stud_fraction = 0.15                     # 壁面積の 15 % が木の柱
a_ins = area * (1.0 - stud_fraction)
a_stud = area * stud_fraction
k_wood = 0.13

r_ins_only = r_cond_plane(0.100, 0.040, a_ins)      # 断熱部
r_stud = r_cond_plane(0.100, k_wood, a_stud)        # 柱部(並列枝)
r_parallel = 1.0 / (1.0 / r_ins_only + 1.0 / r_stud)
r_ins_ideal = r_cond_plane(0.100, 0.040, area)

print(f"  断熱材だけの層抵抗            {r_ins_ideal:.5f} K/W")
print(f"  柱を含めた並列合成            {r_parallel:.5f} K/W"
      f"   ({100 * (1 - r_parallel / r_ins_ideal):.0f} % 悪化)")

net2 = ThermalNetwork("wall_with_stud")
net2.add_node("room", fixed=t_in)
net2.add_node("outdoor", fixed=t_out)
for n in ("surf_in", "n1", "n2", "surf_out"):
    net2.add_node(n)
net2.add_resistor("room", "surf_in", r_conv(h_in, area), "inside film")
net2.add_resistor("surf_in", "n1", r_cond_plane(0.0125, 0.22, area), "gypsum")
net2.add_resistor("n1", "n2", r_ins_only, "glass wool")     # 並列枝 1
net2.add_resistor("n1", "n2", r_stud, "wood stud")          # 並列枝 2
net2.add_resistor("n2", "surf_out", r_cond_plane(0.012, 0.15, area), "plywood")
net2.add_resistor("surf_out", "outdoor", r_conv(h_out, area), "outside film")
sol2 = net2.solve_steady(guess=0.0)
q2 = sol2.flow("room", "surf_in")
print(f"\n  熱橋なしの熱損失 {q_wall:6.1f} W")
print(f"  熱橋ありの熱損失 {q2:6.1f} W   (+{100 * (q2 / q_wall - 1):.0f} %)")
q_stud = sol2.flow_of("wood stud")
q_ins = sol2.flow_of("glass wool")
print(f"    うち断熱材を通る分 {q_ins:5.1f} W")
print(f"    うち柱を通る分     {q_stud:5.1f} W ... 面積比わずか 15 % の柱が"
      f"全体の {100 * q_stud / q2:.0f} % を運ぶ")
print("""
  並列枝を 1 本足すだけで熱橋を定量化できる。CFD を回すまでもない。
  逆に言うと、CFD の 3D 効果(拡がり・回り込み)の多くは
  「並列枝」または「追加の直列抵抗」として 1D に押し込める。これが 1D 化の本質。
""")

# ---------------------------------------------------------------- 図
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

x_pos, temps, labels = [], [], []
x = 0.0
positions = [0.0]
thicknesses = [0.02, 0.0125, 0.100, 0.012, 0.02]
for th in thicknesses:
    x += th
    positions.append(x)
axes[0].plot(np.array(positions) * 1000, [sol[n] for n in path], "o-", lw=2,
             color="#4C72B0")
for xp, name in zip(np.array(positions) * 1000, path):
    axes[0].annotate(name, (xp, sol[name]), textcoords="offset points",
                     xytext=(0, 8), fontsize=8, ha="center")
axes[0].set_xlabel("position through the wall [mm] (films drawn as 20 mm)")
axes[0].set_ylabel("temperature [degC]")
axes[0].set_title("Temperature profile: the big drop = dominant resistance")
axes[0].grid(alpha=0.3)

frac = np.linspace(0.0, 0.30, 100)
loss = []
for f in frac:
    ri = r_cond_plane(0.100, 0.040, area * (1 - f)) if f < 1 else np.inf
    rs = r_cond_plane(0.100, k_wood, area * f) if f > 0 else np.inf
    rp = 1.0 / (1.0 / ri + 1.0 / rs)
    rt = (r_conv(h_in, area) + r_cond_plane(0.0125, 0.22, area) + rp
          + r_cond_plane(0.012, 0.15, area) + r_conv(h_out, area))
    loss.append((t_in - t_out) / rt)
axes[1].plot(frac * 100, loss, lw=2, color="#C44E52")
axes[1].set_xlabel("stud area fraction [%]")
axes[1].set_ylabel("wall heat loss [W]")
axes[1].set_title("Thermal bridge = parallel branch")
axes[1].grid(alpha=0.3)

savefig(fig, "lesson02_series_parallel.png")

print("""
まとめ
  * 直列は足し算、並列は逆数の足し算。合成の順序を間違えないこと。
  * 温度降下の内訳 = 抵抗の内訳。改善対象の特定に直結する。
  * 3D の効果は「並列枝」「追加抵抗」として 1D に組み込める。
次: lesson03  節点法(行列を自分で組む)
""")
