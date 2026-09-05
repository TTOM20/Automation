"""lesson05: 過渡応答 --- 熱容量・時定数・陰解法

学ぶこと
  1. 熱容量 C = rho cp V を入れると RC 回路になり、時定数 tau = R C が出ること
  2. 1 次遅れの解析解と数値解の一致確認
  3. 2 段 RC(ダイ + ヒートシンク)は「速い時定数」と「遅い時定数」を持つこと
  4. 陽解法(前進オイラー)の安定限界と、陰解法(後退オイラー)の無条件安定

実行:  python lessons/lesson05_transient.py
"""

from _util import header, section, savefig, plt
import numpy as np

from thermalnet import ThermalNetwork, capacitance, r_conv

header("lesson05  過渡応答と時定数")

# ---------------------------------------------------------------- 1 次系
section("(1) 集中定数 1 ノード: アルミ塊の冷却")

side = 0.03                       # 立方体 30 mm
volume = side ** 3
area = 6 * side ** 2
rho, cp, k_al = 2700.0, 900.0, 167.0
h = 25.0
t_air, t_init = 25.0, 200.0

cap = capacitance(rho, cp, volume)
r_surf = r_conv(h, area)
tau = r_surf * cap

print(f"  熱容量 C = rho cp V = {rho} x {cp} x {volume:.2e} = {cap:.1f} J/K")
print(f"  表面抵抗 R = 1/(hA) = {r_surf:.4f} K/W")
print(f"  時定数 tau = R C = {tau:.1f} s  (=約 {tau / 60:.1f} 分)")
print("  1 次遅れ: T(t) = T_air + (T0 - T_air) exp(-t/tau)")

net = ThermalNetwork("block")
net.add_node("block", capacitance=cap)
net.add_node("air", fixed=t_air)
net.add_resistor("block", "air", r_surf, "surface")
tr = net.solve_transient(t_end=5 * tau, dt=tau / 200, initial={"block": t_init})

t_analytic = t_air + (t_init - t_air) * np.exp(-tr.time / tau)
err = float(np.max(np.abs(tr["block"] - t_analytic)))
tr_cn = net.solve_transient(t_end=5 * tau, dt=tau / 200,
                            initial={"block": t_init}, theta=0.5)
err_cn = float(np.max(np.abs(tr_cn["block"]
                            - (t_air + (t_init - t_air) * np.exp(-tr_cn.time / tau)))))
print(f"\n  解析解との最大差: 後退オイラー(1 次精度)   {err:.4f} K")
print(f"                    クランク・ニコルソン(2 次) {err_cn:.4f} K"
      f"   <- 同じ dt で 3 桁改善")
print("   時刻 [s]   数値 [degC]   解析 [degC]")
for frac in (0.5, 1.0, 2.0, 3.0):
    i = int(np.argmin(np.abs(tr.time - frac * tau)))
    print(f"   {tr.time[i]:8.1f}   {tr['block'][i]:10.2f}   {t_analytic[i]:10.2f}"
          f"   (t = {frac} tau)")
print("  t = tau で温度差は 63 % 減る。tau は「熱の応答の速さ」そのもの。")

# ---------------------------------------------------------------- 2 次系
section("(2) 2 段 RC: 半導体ダイ + ヒートシンク")

# ダイは小さく軽い = 速い、ヒートシンクは大きく重い = 遅い
c_die = capacitance(2330.0, 700.0, 0.010 * 0.010 * 0.0007)   # Si 10x10x0.7 mm
c_sink = capacitance(2700.0, 900.0, 0.08 * 0.08 * 0.02)      # Al ブロック相当
r_die_sink = 0.35     # TIM 込み [K/W]
r_sink_air = 1.2      # ヒートシンク -> 空気 [K/W]

net2 = ThermalNetwork("die_sink")
net2.add_node("die", capacitance=c_die, heat=0.0)
net2.add_node("sink", capacitance=c_sink)
net2.add_node("air", fixed=t_air)
net2.add_resistor("die", "sink", r_die_sink, "TIM")
net2.add_resistor("sink", "air", r_sink_air, "heatsink")

print(f"  C_die  = {c_die:8.3f} J/K,  tau_die  ~ R_die_sink C_die  = "
      f"{r_die_sink * c_die:6.3f} s")
print(f"  C_sink = {c_sink:8.1f} J/K,  tau_sink ~ R_sink_air C_sink = "
      f"{r_sink_air * c_sink:6.1f} s")

power = 40.0
net2.set_heat("die", lambda t: power if t >= 0 else 0.0)
tr2 = net2.solve_transient(t_end=1200.0, dt=0.05, initial=t_air, store_every=4)

# 系の固有時定数(行列の固有値から)。CFD には無い、回路モデルならではの分析。
caps = np.array([c_die, c_sink])
a_mat = np.array([[1 / r_die_sink, -1 / r_die_sink],
                  [-1 / r_die_sink, 1 / r_die_sink + 1 / r_sink_air]])
eigs = np.linalg.eigvals(np.diag(1 / caps) @ a_mat)
print(f"  行列の固有値から出る時定数 = "
      f"{', '.join(f'{1 / e:.3f} s' for e in sorted(eigs, key=lambda x: 1 / x))}")
t_ss = t_air + power * (r_die_sink + r_sink_air)
t_1200 = tr2.final()["die"]
print(f"  t = 1200 s (= {1200 / (r_sink_air * c_sink):.1f} tau_sink) 時点: "
      f"T_die = {t_1200:.2f} degC, T_sink = {tr2.final()['sink']:.2f} degC")
print(f"  定常回路の答え: T_die = {t_ss:.2f} degC"
      f"   -> まだ {100 * (t_1200 - t_air) / (t_ss - t_air):.1f} % までしか到達していない")
print("  過渡計算は「一番遅い時定数の 4-5 倍」回さないと定常に届かない、の実例。")
print("""
  短時間の負荷変動(バースト)ではダイの小さな熱容量だけが効き、
  分単位の連続負荷ではヒートシンクの大きな熱容量が効く。
  この「時間スケールの分離」が見えるのが熱回路モデルの強み。
  制御・保護ロジック(サーマルスロットリング)の設計はここを使う。""")

section("(3) バースト負荷: 30 秒だけ 120 W")


def burst(t):
    return 120.0 if 300.0 <= t < 330.0 else 40.0


net2.set_heat("die", burst)
tr3 = net2.solve_transient(t_end=900.0, dt=0.05, initial=t_air, store_every=4)
i_peak = int(np.argmax(tr3["die"]))
print(f"  ピーク T_die = {tr3['die'][i_peak]:.2f} degC (t = {tr3.time[i_peak]:.0f} s)")
print(f"  同じ 120 W を定常で流したら "
      f"{t_air + 120 * (r_die_sink + r_sink_air):.1f} degC まで上がる")
print("  短時間なら熱容量が肩代わりしてくれる = 定常設計は過剰になりがち。")

# ---------------------------------------------------------------- 安定性
section("(4) 陽解法 vs 陰解法 --- 時間刻みの安定限界")

net3 = ThermalNetwork("stability")
net3.add_node("block", capacitance=cap)
net3.add_node("air", fixed=t_air)
net3.add_resistor("block", "air", r_surf)
dt_limit = 2.0 * tau
print(f"  この 1 ノード系の陽解法の安定限界 dt < 2 R C = {dt_limit:.1f} s")
for dt in (0.5 * dt_limit, 0.99 * dt_limit, 1.2 * dt_limit):
    tr_e = net3.solve_transient(3000.0, dt, initial={"block": t_init}, theta=0.0)
    peak = float(np.max(np.abs(tr_e["block"] - t_air)))
    verdict = "安定" if peak <= (t_init - t_air) + 1e-6 else "発散!"
    print(f"    dt = {dt:7.1f} s ({dt / dt_limit:.2f} x 限界)  最大振幅"
          f" {peak:10.2f} K  -> {verdict}")
tr_i = net3.solve_transient(6 * 5 * dt_limit, 5 * dt_limit,
                            initial={"block": t_init}, theta=1.0)
print(f"  陰解法は dt = {5 * dt_limit:.0f} s (限界の 5 倍) でも発散しない: "
      f"最終値 {tr_i['block'][-1]:.2f} degC (真値 {t_air:.2f})")
i_mid = 1
t_true_mid = t_air + (t_init - t_air) * np.exp(-tr_i.time[i_mid] / tau)
print(f"  ただし精度は別問題: t = {tr_i.time[i_mid]:.0f} s で数値 "
      f"{tr_i['block'][i_mid]:.1f} degC vs 真値 {t_true_mid:.1f} degC。")
print("  「安定 = 正しい」ではない。dt は tau/10 程度を目安に。")
print("""
  CFD の非定常計算と同じ話。陰解法は無条件安定だが各ステップで連立方程式を
  解く。熱回路は行列が小さいので陰解法が事実上ただ同然 -> 既定は theta=1。""")

# ---------------------------------------------------------------- 図
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

axes[0].plot(tr.time, tr["block"], lw=2, label="numerical")
axes[0].plot(tr.time, t_analytic, "k--", lw=1.2, label="analytical")
axes[0].axvline(tau, color="gray", ls=":")
axes[0].text(tau * 1.05, 150, "t = tau", color="gray")
axes[0].set_xlabel("time [s]")
axes[0].set_ylabel("temperature [degC]")
axes[0].set_title("First order lumped cooling")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].semilogx(tr2.time[1:], tr2["die"][1:], lw=2, label="die")
axes[1].semilogx(tr2.time[1:], tr2["sink"][1:], lw=2, label="heatsink")
axes[1].set_xlabel("time [s] (log)")
axes[1].set_ylabel("temperature [degC]")
axes[1].set_title("Two time constants (40 W step)")
axes[1].legend()
axes[1].grid(alpha=0.3, which="both")

axes[2].plot(tr3.time, tr3["die"], lw=2, label="die")
axes[2].plot(tr3.time, tr3["sink"], lw=2, label="heatsink")
axes[2].axvspan(300, 330, color="orange", alpha=0.25, label="120 W burst")
axes[2].set_xlabel("time [s]")
axes[2].set_ylabel("temperature [degC]")
axes[2].set_title("Burst load: capacitance buys headroom")
axes[2].legend()
axes[2].grid(alpha=0.3)

savefig(fig, "lesson05_transient.png")

print("""
まとめ
  * C = rho cp V、tau = R C。時定数の分離が設計の勘所。
  * 定常設計は「一番遅い時定数まで待った状態」。バーストなら緩和できる。
  * 陽解法は dt < 2 R C(一般には 2/lambda_max)。熱回路なら陰解法一択。
次: lesson06  非線形(自然対流 + 輻射)の反復解法
""")
