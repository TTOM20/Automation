"""lesson03: 節点法 --- 行列を自分の手で組んでみる

学ぶこと
  1. 節点方程式 A T = b の作り方(A はコンダクタンス行列 = グラフラプラシアン)
  2. 温度固定境界の入れ方(行の置換)と発熱の入れ方(右辺)
  3. 対称・対角優位という性質、そしてエネルギー保存の検算
  4. 自作コードと thermalnet の答えが完全一致することの確認

ここを理解すると、市販の熱回路ソルバも CFD の連立方程式も
「同じものの規模違い」だと分かる。

実行:  python lessons/lesson03_nodal_method.py
"""

from _util import header, section, savefig, plt
import numpy as np

from thermalnet import ThermalNetwork

header("lesson03  節点法(nodal method)")

print("""
節点 i のエネルギー保存(定常):

    Q_i + sum_j (T_j - T_i)/R_ij = 0

移項して行列形にすると

    sum_j G_ij (T_i - T_j) = Q_i ,  G_ij = 1/R_ij

    A[i,i] = sum_j G_ij      (自分につながる全コンダクタンスの和)
    A[i,j] = -G_ij           (対称)
    b[i]   = Q_i             (発熱は右辺の「電流源」)

温度固定ノード k は、その行を  T_k = T_fixed  に置き換える。
""")

# ---------------------------------------------------------------- 対象回路
section("対象: 発熱する IC → 基板 → 筐体 → 外気(筐体からは輻射も)")

nodes = ["die", "board", "case", "ambient"]
idx = {n: i for i, n in enumerate(nodes)}

resistors = [                 # (ノード a, ノード b, R [K/W], ラベル)
    ("die", "board", 2.5, "die->board"),
    ("die", "case", 8.0, "die->case (air gap)"),
    ("board", "case", 4.0, "board->case"),
    ("board", "ambient", 30.0, "board->air"),
    ("case", "ambient", 3.0, "case->air (conv)"),
    ("case", "ambient", 12.0, "case->air (rad)"),
]
heat = {"die": 6.0}           # [W]
fixed = {"ambient": 25.0}     # [degC]

for a, b, r, label in resistors:
    print(f"    {label:<22s} R = {r:5.2f} K/W")

# ---------------------------------------------------------------- 手組み
section("(1) 行列を手で組んで numpy で解く")

n = len(nodes)
A = np.zeros((n, n))
b = np.zeros(n)

for a, bn, r, _ in resistors:
    i, j = idx[a], idx[bn]
    g = 1.0 / r
    A[i, i] += g
    A[j, j] += g
    A[i, j] -= g
    A[j, i] -= g

for name, q in heat.items():
    b[idx[name]] += q

print("  境界を入れる前の A [W/K] (対称・行和ゼロ = 純粋な拡散演算子)")
with np.printoptions(precision=4, suppress=True):
    print(A)
print(f"  行和 = {np.round(A.sum(axis=1), 12)}  <- すべて 0。"
      "これは「全体を同じ温度にすると熱は流れない」ことの表れ")

for name, t in fixed.items():
    i = idx[name]
    A[i, :] = 0.0
    A[i, i] = 1.0
    b[i] = t

t_manual = np.linalg.solve(A, b)
print("\n  解 [degC]")
for name in nodes:
    print(f"    {name:<10s} {t_manual[idx[name]]:8.3f}")

# ---------------------------------------------------------------- ライブラリ
section("(2) thermalnet で同じ回路を解く")

net = ThermalNetwork("ic_board_case")
net.add_node("die", heat=heat["die"])
net.add_node("board")
net.add_node("case")
net.add_node("ambient", fixed=fixed["ambient"])
for a, bn, r, label in resistors:
    net.add_resistor(a, bn, r, label)
sol = net.solve_steady(guess=25.0)

print("  最大差 = %.2e K  ->  自作コードとライブラリは同一のことをしている"
      % max(abs(sol[nm] - t_manual[idx[nm]]) for nm in nodes))
print()
print(sol.report("節点法の解"))

# ---------------------------------------------------------------- 検算
section("(3) エネルギー保存の検算(必ずやること)")

q_out = (sol.flow_of("board->air") + sol.flow_of("case->air (conv)")
         + sol.flow_of("case->air (rad)"))
print(f"  投入熱量 {heat['die']:.4f} W  /  外気へ出た熱量 {q_out:.4f} W"
      f"  (残差 {abs(q_out - heat['die']):.2e} W)")
print("""
  節点法は必ず保存則を満たす(離散式そのものが保存則だから)。
  逆に、保存が合わないときは回路の作り間違い(浮きノード、単位ミス)である。
""")

section("(4) 熱抵抗の「見え方」: die から見た合成抵抗")
r_ja = net.total_resistance("die", "ambient")
print(f"  R_ja = (T_die - T_amb)/Q = {r_ja:.3f} K/W")
print(f"  検算: {(sol['die'] - 25.0) / heat['die']:.3f} K/W")
print("""  データシートの R_ja / R_jc はこの合成抵抗のこと。
  ただし合成抵抗は経路の分岐比に依存する = 実装条件が変われば変わる。
  「R_ja は材料定数ではなく、系の定数」という点は覚えておく価値がある。""")

# ---------------------------------------------------------------- 図
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

A_plot = np.zeros((n, n))
for a, bn, r, _ in resistors:
    i, j = idx[a], idx[bn]
    g = 1.0 / r
    A_plot[i, i] += g
    A_plot[j, j] += g
    A_plot[i, j] -= g
    A_plot[j, i] -= g
im = axes[0].imshow(A_plot, cmap="RdBu_r",
                    vmin=-np.abs(A_plot).max(), vmax=np.abs(A_plot).max())
axes[0].set_xticks(range(n), nodes, rotation=30)
axes[0].set_yticks(range(n), nodes)
for i in range(n):
    for j in range(n):
        axes[0].text(j, i, f"{A_plot[i, j]:.3f}", ha="center", va="center", fontsize=8)
axes[0].set_title("Conductance matrix A [W/K]")
fig.colorbar(im, ax=axes[0], fraction=0.046)

paths = ["die->board", "die->case (air gap)", "board->case",
         "board->air", "case->air (conv)", "case->air (rad)"]
flows = [sol.flow_of(p) for p in paths]
axes[1].barh(paths, flows, color="#4C72B0")
axes[1].set_xlabel("heat flow [W]")
axes[1].set_title("Where does the 6 W go?")
axes[1].grid(alpha=0.3, axis="x")

savefig(fig, "lesson03_nodal.png")

print("""
まとめ
  * A は対称・対角優位・行和ゼロ。境界条件だけがこの構造を崩す。
  * 発熱は右辺、温度固定は行の置換。
  * 解いたら必ずエネルギー保存を検算する。
次: lesson04  1D 離散化 --- CFD のメッシュはそのまま熱回路になる
""")
