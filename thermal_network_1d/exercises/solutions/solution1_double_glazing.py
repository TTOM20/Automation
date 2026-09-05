"""Q1 解答: 単板ガラス vs 複層ガラス(直列合成の練習)"""

import _path  # noqa: F401
from thermalnet import ThermalNetwork, r_cond_plane, r_conv

AREA, T_IN, T_OUT = 2.0, 20.0, 0.0
H_IN, H_OUT = 9.0, 23.0


def glazing(layers, label):
    """layers = [(名前, 厚さ, 熱伝導率), ...] を直列につないで解く。"""
    net = ThermalNetwork(label)
    net.add_node("room", fixed=T_IN)
    net.add_node("outdoor", fixed=T_OUT)
    path = ["room", "s_in"]
    res = [r_conv(H_IN, AREA)]
    names = ["inside film"]
    for i, (name, th, k) in enumerate(layers):
        path.append(f"n{i}")
        res.append(r_cond_plane(th, k, AREA))
        names.append(name)
    path[-1] = "s_out"
    path.append("outdoor")
    res.append(r_conv(H_OUT, AREA))
    names.append("outside film")
    for n in path[1:-1]:
        net.add_node(n)
    net.add_series(path, res, names)
    sol = net.solve_steady(guess=10.0)
    r_total = sum(res)
    u = 1.0 / (r_total * AREA)
    q = sol.flow_of("inside film")
    print(f"\n[{label}]")
    for name, r in zip(names, res):
        print(f"    {name:<14s} R = {r:8.5f} K/W  ({100 * r / r_total:5.1f} %)")
    print(f"    合計 R = {r_total:.5f} K/W,  U = {u:.2f} W/(m^2 K),  Q = {q:.1f} W")
    print(f"    室内側ガラス表面温度 = {sol['s_in']:.2f} degC")
    return u, q, sol["s_in"]


print("Q1: 複層ガラスの効果")
u1, q1, ts1 = glazing([("glass", 0.006, 1.0)], "単板ガラス 6 mm")
u2, q2, ts2 = glazing([("glass", 0.006, 1.0),
                       ("air gap", 0.012, 0.075),
                       ("glass", 0.006, 1.0)], "複層ガラス 6-12-6")

print(f"""
  U 値      : {u1:.2f} -> {u2:.2f} W/(m^2 K)
  熱損失     : {q1:.1f} -> {q2:.1f} W  ({q1 / q2:.2f} 倍の改善)
  室内側表面温度: {ts1:.1f} -> {ts2:.1f} degC
  室内が 20 degC / 相対湿度 60 % なら露点は約 12 degC。
  単板では表面が露点を下回り結露するが、複層なら回避できる --- という判断に使う。
  なお空気層は「厚くすれば良い」わけではない: 20 mm を超えると層内で
  自然対流が立ち上がり、実効熱伝導率が上がって効果が頭打ちになる。
""")
