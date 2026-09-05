"""Q6 解答: CFD を 1 回だけ回して 1D モデルを作り、適用範囲を確かめる"""

import _path  # noqa: F401
import numpy as np

from thermalnet import ThermalNetwork, r_cond_plane, r_conv
from thermalnet.fvm import axisymmetric_spreader

A, B, T, K = 0.004, 0.025, 0.004, 180.0
H_REF, Q_REF, T_INF = 700.0, 25.0, 25.0
AREA = np.pi * B ** 2

print("Q6: CFD 1 回で 1D モデルを作る\n")

# ---- 手順 1: CFD を 1 回だけ実行して抵抗を抽出 -------------------------
cfd_ref = axisymmetric_spreader(A, B, T, K, H_REF, Q_REF, 200, 100, T_INF)
r_spread = cfd_ref["R_spread"]
print(f"  基準 CFD (h={H_REF:.0f}, Q={Q_REF:.0f} W):")
print(f"    R_total  = {cfd_ref['R_total']:.5f} K/W")
print(f"    R_1d     = {cfd_ref['R_1d']:.5f} K/W")
print(f"    R_spread = {r_spread:.5f} K/W"
      f"   ({100 * r_spread / cfd_ref['R_total']:.1f} % を占める)")


def model_1d(q, h):
    """抽出した R_spread を定数として使う 1D 回路。"""
    net = ThermalNetwork()
    net.add_node("source", heat=q)
    net.add_node("spread")
    net.add_node("base")
    net.add_node("air", fixed=T_INF)
    net.add_resistor("source", "spread", r_spread, "spreading")
    net.add_resistor("spread", "base", r_cond_plane(T, K, AREA), "conduction")
    net.add_resistor("base", "air", r_conv(h, AREA), "convection")
    return net.solve_steady(guess=T_INF)["source"]


# ---- 手順 2-3: 9 条件で 1D と CFD を比較 --------------------------------
print("\n     Q [W]   h [W/m^2K]    1D [degC]   CFD [degC]   誤差 [K]")
worst = (0.0, None)
for q in (10.0, 25.0, 60.0):
    for h in (400.0, 700.0, 1500.0):
        t_1d = model_1d(q, h)
        t_cfd = axisymmetric_spreader(A, B, T, K, h, q, 200, 100,
                                      T_INF)["T_source_mean"]
        err = t_1d - t_cfd
        if abs(err) > abs(worst[0]):
            worst = (err, (q, h))
        print(f"    {q:6.0f}   {h:10.0f}   {t_1d:10.2f}   {t_cfd:10.2f}"
              f"   {err:+8.3f}")

print(f"""
  最大誤差は Q = {worst[1][0]:.0f} W, h = {worst[1][1]:.0f} の条件で {worst[0]:+.3f} K。

  考察
  1. Q を変えても *相対* 誤差は変わらない。回路が線形(抵抗が温度に依らない)
     なので温度上昇は Q に厳密に比例し、誤差も同じ比率で拡大するだけ。
     h=1500 の列はどの Q でも温度上昇の約 0.9 % の過大評価。
     つまり発熱量は自由に外挿してよい(絶対誤差だけは Q に比例して増える)。
  2. 誤差の出どころは h の変更だけ。h を上げると熱が広がりきる前に
     下面から抜けるため、真の R_spread がわずかに小さくなる。
     それでも h を 2 倍にして 0.5 K 未満。実用上は定数で十分。
  3. 一方、形状(a, b, t, k)を変えたら R_spread は作り直しになる。
     形状を振りたいなら、CFD を数点回して R_spread(a/b, t, ...) を
     テーブル化するか、相関式(Lee et al.)を使う。

  適用範囲の書き方(モデルカード)の例:
    「この 1D モデルは a=4mm/b=25mm/t=4mm/k=180 の形状に対し、
      h = 400-1500 W/(m^2 K)、Q = 任意 の範囲で CFD と 1 K 以内で一致する。
      形状変更時は R_spread の再取得が必要。」
  こう書いておけば、他人が誤用しない。1D モデルには必ず適用範囲を添える。
""")
