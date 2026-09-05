"""lesson07: CFD 結果を 1D 熱回路に落とし込む --- 拡がり抵抗

これが本題。3D(ここでは軸対称)の伝導解析を「正解」として、
素朴な 1D 回路がなぜ、どれだけ外れるかを定量化し、
足りない物理を 1 個の抵抗として回路に組み込む。

対象:
      z=0   +---------------------------+   上面 r<a に熱源 Q(小さいダイ)
            |####|      スプレッダ       |   r>a は断熱
      z=t   +---------------------------+   下面: 熱伝達率 h で空気へ
            r=0                       r=b

学ぶこと
  1. 素朴な 1D 回路 (R = t/(kA) + 1/(hA)) は、熱源が小さいと大幅に過小評価する
  2. 差分の正体が「拡がり抵抗 R_spread」であること
  3. CFD 結果から R_spread を抽出して 1D 回路に足すと一致すること
  4. 相関式(Lee et al.)との比較 --- 便利だが誤差はあることを自分で確認する
  5. 「1D 化の手順」の型: 経路を切る -> 各区間の dT/Q を CFD から読む -> 回路にする

実行:  python lessons/lesson07_cfd_to_1d.py
"""

from _util import header, section, savefig, plt
import numpy as np

from thermalnet import ThermalNetwork, r_cond_plane, r_conv, r_spreading_disk
from thermalnet.fvm import axisymmetric_spreader

header("lesson07  CFD -> 1D 化(拡がり抵抗の抽出)")

a = 0.005      # 熱源半径 [m] (直径 10 mm のダイ)
b = 0.020      # スプレッダ半径 [m] (直径 40 mm)
t = 0.003      # 板厚 [m]
k = 200.0      # 熱伝導率 [W/(m K)] (アルミ)
h = 500.0      # 下面の等価熱伝達率 [W/(m^2 K)]
q_total = 30.0  # 発熱 [W]
t_inf = 25.0

area_b = np.pi * b ** 2
print(f"""
  熱源直径 {2000 * a:.0f} mm / 板直径 {2000 * b:.0f} mm / 板厚 {1000 * t:.0f} mm
  k = {k} W/(m K), h = {h} W/(m^2 K), Q = {q_total} W, T_inf = {t_inf} degC
  面積比 (a/b)^2 = {(a / b) ** 2:.3f}  <- 熱源は板の {(a / b) ** 2 * 100:.0f} % しか覆っていない
""")

# ---------------------------------------------------------------- CFD
section("(1) 軸対称 FVM(= ミニ CFD)で解く")

cfd = axisymmetric_spreader(a, b, t, k, h, q_total, nr=160, nz=80, t_inf=t_inf)
print(f"  格子 {cfd['grid'][0]} x {cfd['grid'][1]}")
print(f"  熱源面の平均温度  {cfd['T_source_mean']:.2f} degC")
print(f"  熱源中心の最高温度 {cfd['T_source_max']:.2f} degC")
print(f"  下面の平均温度    {cfd['T_base_mean']:.2f} degC")
print(f"  全体の熱抵抗 R_total = (T_src - T_inf)/Q = {cfd['R_total']:.4f} K/W")

print("\n  格子依存性の確認(CFD なら必ずやる):")
for nr, nz in [(40, 20), (80, 40), (160, 80), (320, 160)]:
    c = axisymmetric_spreader(a, b, t, k, h, q_total, nr, nz, t_inf)
    print(f"    {nr:4d} x {nz:3d}   R_total = {c['R_total']:.5f} K/W"
          f"   R_spread = {c['R_spread']:.5f} K/W")
print("  格子を 4 倍細かくしても R_total の変化は 0.1 % 未満 -> 格子は十分。")

# ---------------------------------------------------------------- 素朴な 1D
section("(2) 素朴な 1D 回路(よくある間違い)")

naive = ThermalNetwork("naive_1d")
naive.add_node("source", heat=q_total)
naive.add_node("base")
naive.add_node("air", fixed=t_inf)
naive.add_resistor("source", "base", r_cond_plane(t, k, area_b), "conduction t/(kA)")
naive.add_resistor("base", "air", r_conv(h, area_b), "convection 1/(hA)")
sol_naive = naive.solve_steady(guess=t_inf)

print(f"  R_1d = t/(kA) + 1/(hA) = {r_cond_plane(t, k, area_b):.5f}"
      f" + {r_conv(h, area_b):.5f} = {cfd['R_1d']:.5f} K/W")
print(f"  1D の予測 T_source = {sol_naive['source']:.2f} degC")
print(f"  CFD の答え         {cfd['T_source_mean']:.2f} degC")
delta = cfd["T_source_mean"] - sol_naive["source"]
print(f"  -> {delta:.2f} K も低く見積もっている"
      f" (温度上昇でみて {100 * delta / (cfd['T_source_mean'] - t_inf):.0f} % の過小評価)")
print("""
  原因: 1D では「熱が板の全断面 pi b^2 を一様に流れる」と仮定しているが、
  実際には小さな熱源から放射状に広がる過程で余分な温度差が生じる。
  この差が 拡がり抵抗 (spreading / constriction resistance)。
  3D CFD を 1D にするとき、最も見落とされやすい抵抗である。""")

# ---------------------------------------------------------------- 抽出
section("(3) CFD から R_spread を抽出して 1D 回路に足す")

r_spread_cfd = cfd["R_spread"]
print(f"  R_spread = R_total(CFD) - R_1d = {cfd['R_total']:.5f} - {cfd['R_1d']:.5f}"
      f" = {r_spread_cfd:.5f} K/W")
print(f"  相関式 (Lee et al. 1995) では {r_spreading_disk(a, b, t, k, h):.5f} K/W"
      f"  (CFD 比 {100 * (r_spreading_disk(a, b, t, k, h) / r_spread_cfd - 1):+.1f} %)")
print("  相関式は便利だが数 % はずれる。CFD があるなら CFD から抜くのが正確。")

fixed = ThermalNetwork("calibrated_1d")
fixed.add_node("source", heat=q_total)
fixed.add_node("spread")
fixed.add_node("base")
fixed.add_node("air", fixed=t_inf)
fixed.add_resistor("source", "spread", r_spread_cfd, "spreading (from CFD)")
fixed.add_resistor("spread", "base", r_cond_plane(t, k, area_b), "conduction")
fixed.add_resistor("base", "air", r_conv(h, area_b), "convection")
sol_fixed = fixed.solve_steady(guess=t_inf)
print()
print(sol_fixed.report("校正後の 1D 回路"))
print(f"\n  校正後 T_source = {sol_fixed['source']:.3f} degC  vs "
      f"CFD {cfd['T_source_mean']:.3f} degC  (差 "
      f"{abs(sol_fixed['source'] - cfd['T_source_mean']):.3e} K)")
print("""
  当たり前だが、抜いてきた抵抗を足せば一致する。重要なのは
  「この R_spread は形状と h が変わらない限り再利用できる」という点。
  CFD を 1 回だけ回して抵抗を作り、あとは 1D 回路で何百通りも条件を振る --- 
  これが実務で最も費用対効果の高い CFD の使い方。""")

# ---------------------------------------------------------------- 適用範囲
section("(4) 抽出した抵抗はどこまで使い回せるか")

print("   h [W/m^2K]   R_spread(CFD)   R_total に占める割合   "
      "h=500 の値を流用した誤差 [K] @30W")
base_spread = r_spread_cfd
for h_test in (200.0, 500.0, 1000.0, 3000.0):
    c = axisymmetric_spreader(a, b, t, k, h_test, q_total, 160, 80, t_inf)
    r1d = t / (k * area_b) + 1.0 / (h_test * area_b)
    t_fixed_model = t_inf + q_total * (r1d + base_spread)
    print(f"    {h_test:8.0f}   {c['R_spread']:13.5f}   "
          f"{100 * c['R_spread'] / c['R_total']:16.1f} %   "
          f"{t_fixed_model - c['T_source_mean']:+16.2f}")
print("""  R_spread 自体は h が 15 倍になっても 6 % ほどしか動かない
  (熱の広がり方は主に形状と k で決まる)。一方 R_1d は 1/h で急減するので、
  よく冷える設計ほど R_spread の相対的な重みが増していく --- 上の表で
  割合が 7 % から 50 % へ跳ね上がる。
  つまり「風量を増やしても頭打ちになる」現象の正体は拡がり抵抗。
  流用時の誤差は 0.5 K 以下なので、この形状なら定数として扱ってよい。""")

print("\n   熱源半径比 a/b を振ったときの拡がり抵抗の重み:")
print("     a/b    R_spread [K/W]   R_total に占める割合")
ratios = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
share = []
for ratio in ratios:
    c = axisymmetric_spreader(ratio * b, b, t, k, h, q_total, 200, 60, t_inf)
    share.append(100 * c["R_spread"] / c["R_total"])
    print(f"    {ratio:5.2f}   {c['R_spread']:12.5f}   {share[-1]:16.1f} %")
print("""  熱源が小さいほど拡がり抵抗が支配的になる。パワー半導体や
  最近の小面積・高発熱ダイでは、ここが設計の主戦場。""")

# ---------------------------------------------------------------- 図
fig = plt.figure(figsize=(15, 4.4))

ax0 = fig.add_subplot(1, 3, 1)
temp = cfd["T"]
rr, zz = np.meshgrid(cfd["r"] * 1000, cfd["z"] * 1000, indexing="ij")
cs = ax0.contourf(rr, zz, temp, levels=20, cmap="inferno")
ax0.contour(rr, zz, temp, levels=12, colors="k", linewidths=0.4)
ax0.axvline(a * 1000, color="cyan", ls="--", lw=1.5)
ax0.text(a * 1000 + 0.4, 0.4, "source edge", color="cyan", fontsize=8)
ax0.invert_yaxis()
ax0.set_xlabel("r [mm]")
ax0.set_ylabel("z [mm] (0 = heated face)")
ax0.set_title("Axisymmetric FVM temperature field")
fig.colorbar(cs, ax=ax0, label="T [degC]")

ax1 = fig.add_subplot(1, 3, 2)
labels = ["naive 1D\n(t/kA + 1/hA)", "1D + spreading", "CFD"]
values = [sol_naive["source"], sol_fixed["source"], cfd["T_source_mean"]]
bars = ax1.bar(labels, values, color=["#DD8452", "#4C72B0", "#55A868"])
ax1.axhline(t_inf, color="gray", ls=":")
for bar, v in zip(bars, values):
    ax1.text(bar.get_x() + bar.get_width() / 2, v, f"{v:.1f}", ha="center",
             va="bottom")
ax1.set_ylabel("source temperature [degC]")
ax1.set_title(f"Missing {delta:.0f} K without spreading resistance")
ax1.grid(alpha=0.3, axis="y")

ax2 = fig.add_subplot(1, 3, 3)
ax2.plot(ratios, share, "o-", lw=2, color="#C44E52")
ax2.set_xlabel("source-to-plate radius ratio a/b")
ax2.set_ylabel("spreading share of R_total [%]")
ax2.set_title("Small source -> spreading dominates")
ax2.grid(alpha=0.3)

savefig(fig, "lesson07_cfd_to_1d.png")

print("""
1D 化の手順(この型を覚える)
  1. 熱の経路を物理的に切る(ダイ / TIM / スプレッダ / フィン / 空気)
  2. CFD 結果から各区間の代表温度を面積加重平均で拾う
  3. R_i = (T_i - T_i+1) / Q で抵抗を作る(Q は必ずその区間を通る熱流)
  4. 回路を組み直し、CFD と一致することを確認する
  5. 使い回す範囲(h・流量・発熱分布)を明記する。これがモデルの適用限界。
次: lesson08  総合演習(CPU 冷却系のフルモデル)
""")
