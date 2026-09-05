"""Q5 解答: 表面温度 45 degC 以下に保てる最大発熱量(非線形 + 二分探索)"""

import _path  # noqa: F401
from thermalnet import ThermalNetwork, h_radiation, r_conv
from thermalnet.fluids import h_natural_horizontal_plate_up, h_natural_vertical_plate

W, D, HGT = 0.200, 0.150, 0.050
A_SIDE = 2 * (W * HGT) + 2 * (D * HGT)
A_TOP = W * D
LC_TOP = A_TOP / (2 * (W + D))
R_INTERNAL = 0.8
T_AMB, T_LIMIT = 25.0, 45.0


def case_temperature(power, emissivity, top_area_factor=1.0):
    """発熱 power [W] のときの筐体表面温度 [degC]。"""
    a_top = A_TOP * top_area_factor          # フィンで実効面積を稼いだ場合
    net = ThermalNetwork()
    net.add_node("component", heat=power)
    net.add_node("case")
    net.add_node("ambient", fixed=T_AMB)
    net.add_resistor("component", "case", R_INTERNAL)
    net.add_resistor("case", "ambient",
                     lambda ts, ta: r_conv(h_natural_vertical_plate(ts, ta, HGT),
                                           A_SIDE), "conv side")
    net.add_resistor("case", "ambient",
                     lambda ts, ta: r_conv(
                         h_natural_horizontal_plate_up(ts, ta, LC_TOP), a_top),
                     "conv top")
    net.add_resistor("case", "ambient",
                     lambda ts, ta: r_conv(h_radiation(emissivity, ts, ta),
                                           A_SIDE + A_TOP), "radiation")
    return net.solve_steady(guess=T_AMB + 15.0, relaxation=0.7)["case"]


def max_power(emissivity, top_area_factor=1.0, lo=0.1, hi=200.0, tol=1e-3):
    """T_case = 45 degC となる発熱量を二分探索で求める。

    非線形なので解析的には解けない。「解いて評価 -> 挟み込む」が定石。
    """
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if case_temperature(mid, emissivity, top_area_factor) > T_LIMIT:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


print("Q5: 表面 45 degC 以下に保てる最大発熱量\n")
print("     条件                                  上限 [W]   そのときの内部温度")
for label, eps, factor in [
        ("黒アルマイト eps=0.85", 0.85, 1.0),
        ("無処理アルミ eps=0.05", 0.05, 1.0),
        ("eps=0.85 + 上面フィン(実効 3 倍)", 0.85, 3.0),
        ("eps=0.05 + 上面フィン(実効 3 倍)", 0.05, 3.0)]:
    p = max_power(eps, factor)
    t_comp = p * R_INTERNAL + case_temperature(p, eps, factor)
    print(f"    {label:<36s} {p:7.2f}   {t_comp:8.1f} degC")

print("""
  読み取れること
  * 放射率を 0.05 -> 0.85 にするだけで許容発熱がおよそ 2 倍になる。
    塗装は最も安価な熱対策。特に自然空冷では効果が大きい。
  * 上面にフィンを足しても、側面と輻射が変わらないので効果は限定的。
    「表面積 3 倍 = 3 倍冷える」ではない --- 並列回路のうち 1 本を
    太くしても、他の枝が細いままなら合成抵抗はそこまで下がらない。
  * 触感温度の規格(例: IEC 60950 系の金属 55 degC / 樹脂 70 degC など)は
    材質と接触時間で決まる。設計値はプロジェクトの要求に合わせること。

  二分探索は非線形 1D モデルの定番の使い方。CFD でこれをやると
  1 点あたり数時間 x 十数点で現実的でない。
""")
