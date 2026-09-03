"""lesson06: 非線形な熱回路 --- 自然対流 + 輻射を相関式で閉じる

「CFD の 1D 化」で一番の勘所は、CFD が数値的に解いていた境界層を
相関式 h = Nu k / L で置き換えることにある。ところが自然対流も輻射も
h が温度に依存するので、回路が非線形になり反復が要る。

学ぶこと
  1. 抵抗を「温度の関数」として登録し、不動点反復で解く
  2. 緩和係数(under-relaxation)が収束に効くこと
  3. h 一定と仮定した線形モデルがどれだけずれるか
  4. 放射率の効果 --- 自然空冷では塗装するだけで温度が下がる

対象: ファンレス筐体(200 x 150 x 50 mm アルミ)内部で 20 W 発熱

実行:  python lessons/lesson06_nonlinear.py
"""

from _util import header, section, savefig, plt
import numpy as np

from thermalnet import ThermalNetwork, h_radiation, r_conv
from thermalnet.fluids import h_natural_horizontal_plate_up, h_natural_vertical_plate

header("lesson06  非線形熱回路(自然対流 + 輻射)")

w, d, hgt = 0.200, 0.150, 0.050        # 筐体寸法 [m]
a_side = 2 * (w * hgt) + 2 * (d * hgt)  # 側面 [m^2]
a_top = w * d                           # 上面 [m^2]
lc_top = a_top / (2 * (w + d))          # 上面の特性長 A/P [m]
r_internal = 0.8                        # 発熱部品 -> 筐体壁 [K/W](内部の伝導+接触)
t_amb = 25.0
power = 20.0

print(f"""
  側面積 {a_side * 1e4:.0f} cm^2 (高さ {hgt * 1000:.0f} mm)、
  上面積 {a_top * 1e4:.0f} cm^2 (特性長 Lc = A/P = {lc_top * 1000:.1f} mm)
  発熱 {power} W、室温 {t_amb} degC
""")


def build(emissivity: float, power_w: float = power) -> ThermalNetwork:
    """筐体の非線形熱回路を組む。抵抗は (T_surface, T_ambient) の関数。"""
    net = ThermalNetwork(f"enclosure_eps{emissivity}")
    net.add_node("component", heat=power_w)
    net.add_node("case")
    net.add_node("ambient", fixed=t_amb)
    net.add_resistor("component", "case", r_internal, "internal")
    net.add_resistor(
        "case", "ambient",
        lambda ts, ta: r_conv(h_natural_vertical_plate(ts, ta, hgt), a_side),
        "conv side")
    net.add_resistor(
        "case", "ambient",
        lambda ts, ta: r_conv(h_natural_horizontal_plate_up(ts, ta, lc_top), a_top),
        "conv top")
    net.add_resistor(
        "case", "ambient",
        lambda ts, ta: r_conv(h_radiation(emissivity, ts, ta), a_side + a_top),
        "radiation")
    return net


# ---------------------------------------------------------------- 反復の中身
section("(1) 不動点反復で解く(黒アルマイト eps = 0.85)")

net = build(0.85)
sol = net.solve_steady(guess=t_amb + 20.0, relaxation=0.7, tol=1e-9)
print(f"  反復回数 {sol.iterations}、最終残差 {sol.residual:.2e} K")
print()
print(sol.report("非線形解"))

h_side = h_natural_vertical_plate(sol["case"], t_amb, hgt)
h_top = h_natural_horizontal_plate_up(sol["case"], t_amb, lc_top)
h_rad = h_radiation(0.85, sol["case"], t_amb)
print(f"\n  収束時の熱伝達率: 側面 h = {h_side:.2f}, 上面 h = {h_top:.2f}, "
      f"輻射 h_r = {h_rad:.2f}  [W/(m^2 K)]")
print(f"  輻射の分担 = {100 * sol.flow_of('radiation') / power:.0f} %"
      "  <- 自然空冷では輻射が 3-5 割を占めるのが普通。無視できない")

section("(2) 緩和係数と収束の様子")

print("   緩和係数   反復回数   T_case [degC]")
for relax in (1.0, 0.9, 0.7, 0.5, 0.3):
    s = build(0.85).solve_steady(guess=t_amb + 20.0, relaxation=relax, tol=1e-9,
                                 max_iter=500)
    print(f"    {relax:6.1f}   {s.iterations:8d}   {s['case']:10.3f}")
print("""  この問題は素直なので緩和なし(1.0)でも収束するが、
  h の温度依存が強い/多ノードの系では発散しやすい。緩和はその保険。
  CFD の under-relaxation factor とまったく同じ発想。""")

section("(3) h 一定と仮定した線形モデルとの比較")

h_fixed_list = [5.0, 8.0, 11.0]
print("   仮定 h [W/m^2K]   T_component [degC]   非線形解との差 [K]")
for h_fixed in h_fixed_list:
    lin = ThermalNetwork()
    lin.add_node("component", heat=power)
    lin.add_node("case")
    lin.add_node("ambient", fixed=t_amb)
    lin.add_resistor("component", "case", r_internal)
    lin.add_resistor("case", "ambient", r_conv(h_fixed, a_side + a_top))
    lin.add_resistor("case", "ambient",
                     r_conv(h_radiation(0.85, 60.0, t_amb), a_side + a_top))
    s = lin.solve_steady(guess=t_amb)
    print(f"    {h_fixed:10.1f}        {s['component']:12.2f}"
          f"        {s['component'] - sol['component']:+8.2f}")
print("""  h を決め打ちすると簡単に 10 K 単位でずれる。
  1D 化の精度は「相関式をどれだけ的確に選ぶか」でほぼ決まる、という教訓。
  CFD を持っているなら、CFD の結果から h を逆算して校正する(lesson07)。""")

# ---------------------------------------------------------------- 図
section("(4) 発熱量スイープと放射率の効果")

powers = np.linspace(2.0, 40.0, 25)
curves = {}
for eps in (0.05, 0.30, 0.85):
    temps = []
    for p in powers:
        s = build(eps, p).solve_steady(guess=t_amb + 20.0, relaxation=0.7)
        temps.append(s["component"])
    curves[eps] = np.array(temps)
    print(f"  eps = {eps:.2f}: {power:.0f} W のとき T_component = "
          f"{np.interp(power, powers, temps):.1f} degC")
print(f"  無処理アルミ(eps=0.05)から黒色塗装(eps=0.85)にするだけで "
      f"{np.interp(power, powers, curves[0.05]) - np.interp(power, powers, curves[0.85]):.1f} K 低下")

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

for eps, temps in curves.items():
    axes[0].plot(powers, temps, lw=2, label=f"emissivity {eps:.2f}")
axes[0].plot(powers, t_amb + powers * (sol["component"] - t_amb) / power, "k--",
             lw=1.2, label="linear extrapolation")
axes[0].set_xlabel("dissipated power [W]")
axes[0].set_ylabel("component temperature [degC]")
axes[0].set_title("Nonlinear: h grows with dT, so T is sub-linear")
axes[0].legend()
axes[0].grid(alpha=0.3)

t_scan = np.linspace(30, 120, 100)
axes[1].plot(t_scan, [h_natural_vertical_plate(t, t_amb, hgt) for t in t_scan],
             lw=2, label="natural convection (side)")
axes[1].plot(t_scan, [h_natural_horizontal_plate_up(t, t_amb, lc_top) for t in t_scan],
             lw=2, label="natural convection (top)")
for eps in (0.05, 0.85):
    axes[1].plot(t_scan, [h_radiation(eps, t, t_amb) for t in t_scan], "--",
                 lw=2, label=f"radiation (eps={eps})")
axes[1].set_xlabel("surface temperature [degC]")
axes[1].set_ylabel("h [W/(m^2 K)]")
axes[1].set_title("Heat transfer coefficients are functions of T")
axes[1].legend(fontsize=8)
axes[1].grid(alpha=0.3)

savefig(fig, "lesson06_nonlinear.png")

print("""
まとめ
  * 抵抗を T の関数にすれば非線形回路も同じ枠組みで解ける(反復するだけ)。
  * 自然対流 h ~ dT^0.25 なので、T は発熱量に対して線形には上がらない。
  * 自然空冷では輻射が 3-5 割。放射率を上げるのは最も安い対策。
  * 1D 化の誤差は格子ではなく相関式の選択で決まる。
次: lesson07  CFD 結果を 1D 回路に落とす(拡がり抵抗)
""")
