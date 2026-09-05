"""lesson01: 熱抵抗とは何か --- オームの法則アナロジー

学ぶこと
  1. 熱回路の 4 つの部品(伝導・対流・輻射・熱容量)と単位
  2. 「Q = dT / R」で考えると設計の見通しが一気に良くなること
  3. 輻射は非線形だが、等価熱伝達率 h_r で線形素子に化けること
  4. ビオ数 Bi で「1D 化してよいか」を判定すること

実行:  python lessons/lesson01_thermal_resistance.py
"""

from _util import header, section, savefig, plt
import numpy as np

from thermalnet import (ThermalNetwork, biot_number, h_radiation, r_cond_plane,
                        r_conv, r_radiation)

header("lesson01  熱抵抗の基礎")

print("""
熱回路網とは、伝熱経路を「抵抗のつながり」に置き換えたモデルである。
CFD が偏微分方程式を数百万セルで解くのに対し、熱回路網は数個〜数百個の
代数方程式で解く。設計初期の検討・感度分析・制御用モデルでは後者が強い。

    電気                    熱
    ------------------------------------------------
    電圧 V   [V]            温度 T   [K] または [degC]
    電流 I   [A]            熱流 Q   [W]
    抵抗 R   [ohm]          熱抵抗 R [K/W]
    容量 C   [F]            熱容量 C [J/K]
    I = V/R                 Q = dT/R
""")

# ---------------------------------------------------------------- 具体例
section("例題: 100 x 100 mm のアルミ板 (厚さ 5 mm)、裏面全面のヒータで 15 W")

area = 0.100 * 0.100          # [m^2] (熱源が板全面 = 素直な 1D。小さい熱源だと lesson07 の拡がり抵抗が要る)
thickness = 0.005             # [m]
k_al = 167.0                  # アルミ合金 [W/(m K)]
h_conv = 12.0                 # 自然対流の目安 [W/(m^2 K)]
emissivity = 0.85             # 黒アルマイト処理
q_heat = 15.0                 # [W]
t_air = 25.0                  # [degC]

r_plate = r_cond_plane(thickness, k_al, area)
r_air = r_conv(h_conv, area)
t_surface_guess = 100.0       # 輻射抵抗の評価に使う仮の表面温度(あとで検算する)
r_rad = r_radiation(emissivity, area, t_surface_guess, t_air)

print(f"  伝導抵抗   R_cond = L/(k A) = {thickness}/({k_al} x {area:.5f})"
      f" = {r_plate:.4f} K/W")
print(f"  対流抵抗   R_conv = 1/(h A) = 1/({h_conv} x {area:.5f})"
      f" = {r_air:.4f} K/W")
print(f"  輻射抵抗   R_rad  = 1/(h_r A), h_r = "
      f"{h_radiation(emissivity, t_surface_guess, t_air):.2f} W/(m^2 K)"
      f" -> {r_rad:.4f} K/W")
print("""
  ここが最初の勘所: 伝導抵抗 %.4f K/W に対して表面の抵抗は 3 桁大きい。
  つまり「板を厚くする/材質を銅に変える」より「表面積を稼ぐ/風を当てる」が効く。
  熱回路にすると、どこがボトルネックかが数字で即答できる。
""" % r_plate)

# ---------------------------------------------------------------- 回路にする
section("回路として解く(対流と輻射は並列)")

net = ThermalNetwork("plate")
net.add_node("heater", heat=q_heat)
net.add_node("surface")
net.add_node("air", fixed=t_air)
net.add_resistor("heater", "surface", r_plate, "plate conduction")
net.add_resistor("surface", "air", r_air, "convection")
net.add_resistor("surface", "air", r_rad, "radiation")
print(net.describe())

sol = net.solve_steady(guess=t_air)
print()
print(sol.report("解"))

q_conv = (sol["surface"] - t_air) / r_air
q_rad = (sol["surface"] - t_air) / r_rad
print(f"\n  対流が運ぶ熱 {q_conv:5.2f} W / 輻射が運ぶ熱 {q_rad:5.2f} W"
      f"  (輻射の分担 {100 * q_rad / q_heat:.0f} %)")
print("  自然対流の機器で輻射を無視すると、この分だけ温度を過大に見積もる。")

h_r_check = h_radiation(emissivity, sol["surface"], t_air)
print(f"\n  検算: 輻射の h_r は仮定した表面温度 {t_surface_guess} degC で "
      f"{h_radiation(emissivity, t_surface_guess, t_air):.2f} W/(m^2 K)、")
print(f"        解いた表面温度 {sol['surface']:.1f} degC では {h_r_check:.2f} W/(m^2 K)。")
print("        ずれが大きいときは h_r を更新して解き直す = 非線形反復 (lesson06)。")

# ---------------------------------------------------------------- Bi 数
section("1D 化してよいか? --- ビオ数による判定")

bi = biot_number(h_conv + 1.0 / (r_rad * area), k_al, thickness)
print(f"  Bi = h Lc / k = {bi:.5f}")
print("""  Bi < 0.1 なら物体内部の温度差は表面-流体間の温度差に比べて無視でき、
  「板全体を 1 ノード(集中定数)」として扱ってよい。
  逆に Bi が大きい(樹脂・厚肉・高 h)ときは物体内部を複数ノードに割る。
  これが lesson04 でやる「1D 離散化」に直結する。""")

# ---------------------------------------------------------------- 図
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

names = ["conduction\n(plate)", "convection", "radiation"]
values = [r_plate, r_air, r_rad]
axes[0].bar(names, values, color=["#4C72B0", "#DD8452", "#55A868"])
axes[0].set_ylabel("thermal resistance [K/W]")
axes[0].set_title("Where is the bottleneck?")
for i, v in enumerate(values):
    axes[0].text(i, v, f"{v:.3f}", ha="center", va="bottom")

h_values = np.linspace(3, 60, 200)
t_surf = []
for h in h_values:
    n2 = ThermalNetwork()
    n2.add_node("heater", heat=q_heat)
    n2.add_node("surface")
    n2.add_node("air", fixed=t_air)
    n2.add_resistor("heater", "surface", r_plate)
    n2.add_resistor("surface", "air", r_conv(h, area))
    n2.add_resistor("surface", "air", r_rad)
    t_surf.append(n2.solve_steady(guess=t_air)["heater"])
axes[1].plot(h_values, t_surf, lw=2, color="#C44E52")
axes[1].axhline(85, ls="--", color="gray")
axes[1].text(30, 87, "design limit 85 degC", color="gray")
axes[1].set_xlabel("convective coefficient h [W/(m^2 K)]")
axes[1].set_ylabel("heater temperature [degC]")
axes[1].set_title("Effect of surface cooling (15 W)")
axes[1].grid(alpha=0.3)

savefig(fig, "lesson01_resistance.png")

print("""
まとめ
  * 熱抵抗は「1 W 流したとき何 K 上がるか」。単位 K/W を体に入れる。
  * 直列に並ぶ抵抗のうち最大のものが温度を支配する = 対策すべき場所。
  * 輻射は h_r = eps sigma (Ts+Tsur)(Ts^2+Tsur^2) で線形素子にできる。
    ただし温度依存なので、厳密には反復が要る(lesson06)。
次: lesson02  直列・並列合成と多層壁
""")
