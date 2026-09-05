"""Q3 解答: 実測(またはCFD)データから R と C を同定する

1D モデル作りの現場作業そのもの。ステップ応答さえあれば RC は取れる。
"""

import _path  # noqa: F401
import numpy as np

from thermalnet import ThermalNetwork

RNG = np.random.default_rng(0)
R_TRUE, C_TRUE = 1.35, 180.0        # 同定したい真値(実験では未知)
Q, T_INF = 40.0, 25.0


def measured_data(t_end=1200.0, dt=2.0, noise=0.25):
    """『実測データ』を生成する(真のモデル + 計測ノイズ)。"""
    net = ThermalNetwork()
    net.add_node("part", capacitance=C_TRUE, heat=Q)
    net.add_node("air", fixed=T_INF)
    net.add_resistor("part", "air", R_TRUE)
    tr = net.solve_transient(t_end, dt, initial=T_INF, theta=0.5)
    return tr.time, tr["part"] + RNG.normal(0.0, noise, tr["part"].shape)


def identify(time, temp):
    """R と tau を同定する。

    T(t) = T_inf + Q R (1 - exp(-t/tau))  を変形すると
        ln( (T_ss - T(t)) / (T_ss - T_inf) ) = -t / tau
    と直線になる。傾きから tau、定常値から R が出る。
    """
    t_ss = float(np.mean(temp[-10:]))               # 定常値(末尾平均)
    r_hat = (t_ss - T_INF) / Q
    y = (t_ss - temp) / (t_ss - T_INF)
    mask = (y > 0.05) & (y < 0.95)                  # 端はノイズに弱いので除く
    slope = np.polyfit(time[mask], np.log(y[mask]), 1)[0]
    tau_hat = -1.0 / slope
    return r_hat, tau_hat, tau_hat / r_hat


print("Q3: ステップ応答から R, C を同定する\n")
time, temp = measured_data()
r_hat, tau_hat, c_hat = identify(time, temp)

print(f"  真値      R = {R_TRUE:6.3f} K/W, C = {C_TRUE:6.1f} J/K, "
      f"tau = {R_TRUE * C_TRUE:6.1f} s")
print(f"  同定結果  R = {r_hat:6.3f} K/W, C = {c_hat:6.1f} J/K, tau = {tau_hat:6.1f} s")
print(f"  誤差      R {100 * (r_hat / R_TRUE - 1):+5.1f} %, "
      f"C {100 * (c_hat / C_TRUE - 1):+5.1f} %")

# 同定したモデルで再現して残差を見る
net = ThermalNetwork()
net.add_node("part", capacitance=c_hat, heat=Q)
net.add_node("air", fixed=T_INF)
net.add_resistor("part", "air", r_hat)
tr = net.solve_transient(time[-1], 2.0, initial=T_INF, theta=0.5)
resid = np.interp(time, tr.time, tr["part"]) - temp
print(f"  再現モデルの残差 RMS = {np.sqrt(np.mean(resid ** 2)):.3f} K "
      f"(計測ノイズ 0.25 K + わずかな系統誤差。妥当なフィット)")

print("\n  データ長を変えるとどうなるか:")
print("     取得時間 [s]   tau/データ長   同定 R [K/W]   同定 C [J/K]")
for t_end in (30.0, 120.0, 300.0, 600.0, 1200.0, 2400.0):
    tt, yy = measured_data(t_end=t_end)
    try:
        r_i, tau_i, c_i = identify(tt, yy)
        print(f"      {t_end:8.0f}   {R_TRUE * C_TRUE / t_end:10.2f}   "
              f"{r_i:12.3f}   {c_i:11.1f}")
    except (ValueError, np.linalg.LinAlgError):
        print(f"      {t_end:8.0f}   {R_TRUE * C_TRUE / t_end:10.2f}   同定不能")

print(f"""
  教訓: 定常値 T_ss が取れていないと R も tau も総崩れになる。
  上の表がその実例: 2.5 tau (600 s) までだと R を 9 % 過小評価し、
  5 tau (1200 s) でようやく 1 %、10 tau で誤差はほぼ消える。
  ステップ応答は最低 5 tau まで測ること。この例では tau = {R_TRUE * C_TRUE:.0f} s
  なので 1200 s 以上が目安。

  CFD で過渡解析をするときも同じ。「解析時間が足りないまま定常と見なす」のは
  1D モデル同定で最も多い失敗。lesson05 (2) で見た 96.9 % の話と同じ。
""")
