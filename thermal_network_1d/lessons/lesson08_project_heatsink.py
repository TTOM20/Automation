"""lesson08: 総合演習 --- 空冷 CPU 冷却系のフル 1D モデル

これまでの全要素(直列並列・相関式・拡がり抵抗・流体の昇温・熱容量・非線形)
を 1 本の回路にまとめ、実務でやることを一通りやる:
    定常設計 -> 内訳分析 -> 感度分析 -> 設計スイープ -> 過渡(バースト)検証

熱経路:
    die --R_jc--> IHS --TIM2--> sink base --spread--> fin base
        --fin(eta)--> air(平均) --caloric--> air(入口, 固定)

「air(平均)」ノードが重要。流路を通る空気は熱をもらって昇温するので、
    R_caloric = 1/(2 mdot cp)  (出口までの平均をとると係数 1/2)
を直列に入れないと冷却性能を過大評価する。CFD だと自動で出るが、
1D 化のときは自分で入れる必要がある --- 最も多い 1D 化ミス。

実行:  python lessons/lesson08_project_heatsink.py
"""

from _util import header, section, savefig, plt
import numpy as np

from thermalnet import (ThermalNetwork, capacitance, r_cond_plane, r_fin_array,
                        r_interface, r_spreading_disk)
from thermalnet.fluids import air_properties, h_duct

header("lesson08  総合演習: 空冷ヒートシンクの 1D モデル")

# ---------------------------------------------------------------- 諸元
POWER = 95.0          # TDP [W]
T_IN = 30.0           # ケース内吸気温度 [degC]

DIE = dict(size=0.012, thickness=0.0007, k=120.0)          # Si ダイ 12x12 mm
TIM1 = 0.6e-5                                              # ダイ-IHS [m^2 K/W]
IHS = dict(size=0.035, thickness=0.003, k=390.0)           # 銅 IHS 35x35 mm
TIM2 = 2.0e-5                                              # IHS-ヒートシンク
SINK = dict(base=0.060, base_t=0.006, k=200.0,             # Al ベース 60x60 mm
            fin_t=0.0012, fin_h=0.040, n_fins=22)
AIR = dict(velocity=2.5)                                   # フィン間風速 [m/s]


def air_mass_flow(velocity: float, t_film: float) -> tuple:
    """フィン間流路の質量流量と水力直径を返す。"""
    gap = (SINK["base"] - SINK["n_fins"] * SINK["fin_t"]) / (SINK["n_fins"] - 1)
    flow_area = gap * SINK["fin_h"] * (SINK["n_fins"] - 1)
    d_h = 2.0 * gap * SINK["fin_h"] / (gap + SINK["fin_h"])
    p = air_properties(t_film)
    return p["rho"] * velocity * flow_area, d_h, gap, p


# ファン特性(簡易): dp_fan = P0 (1 - V/V0) の直線近似
FAN = dict(dp_max=30.0, flow_max=6.5e-3)   # [Pa], [m^3/s]


def channel_velocity(n_fins: int, fin_h: float, fixed_pressure: bool,
                     t_film: float = 45.0) -> float:
    """フィン間流速 [m/s]。

    fixed_pressure=False なら AIR["velocity"] をそのまま使う(流速指定)。
    True なら「ファン特性と流路圧損の交点(動作点)」を解く:

        系の圧損(平行平板の層流、f Re = 24):
            dp_sys = 48 mu L u / D_h^2 = K V,   V = u A_flow
        ファン: dp_fan = P0 (1 - V/V0)
        交点:   V = P0 / (K + P0/V0)

    フィンを増やすと隙間が狭まり K が急増して流量が落ちる。
    この連成を入れないと「フィンは多いほど良い」という誤った結論になる
    --- 1D 化で最も多い失敗のひとつ。
    """
    if not fixed_pressure:
        return AIR["velocity"]
    gap = (SINK["base"] - n_fins * SINK["fin_t"]) / (n_fins - 1)
    if gap <= 1e-4:
        return 1e-3
    d_h = 2.0 * gap * fin_h / (gap + fin_h)
    a_flow = (n_fins - 1) * gap * fin_h
    mu = air_properties(t_film)["mu"]
    k_sys = 48.0 * mu * SINK["base"] / (d_h ** 2 * a_flow)   # dp = K V [Pa/(m^3/s)]
    flow = FAN["dp_max"] / (k_sys + FAN["dp_max"] / FAN["flow_max"])
    return flow / a_flow


def build_network(power: float = POWER, velocity: float | None = None,
                  tim2: float = TIM2, fin_h: float = SINK["fin_h"],
                  n_fins: int = SINK["n_fins"], k_sink: float = SINK["k"],
                  t_in: float = T_IN, with_capacitance: bool = False,
                  fixed_pressure: bool = False,
                  verbose: bool = False) -> ThermalNetwork:
    """冷却系の 1D 熱回路を構築する。"""
    if velocity is None:
        velocity = channel_velocity(n_fins, fin_h, fixed_pressure)
    a_die = DIE["size"] ** 2
    a_ihs = IHS["size"] ** 2
    a_base = SINK["base"] ** 2

    # --- 各抵抗 ---------------------------------------------------------
    r_die = r_cond_plane(DIE["thickness"], DIE["k"], a_die)
    r_tim1 = r_interface(TIM1, a_die)
    # ダイ(小) -> IHS(大) の拡がり: 等価半径に直して相関式を使う
    a_eq_die = DIE["size"] / np.sqrt(np.pi)          # 面積等価半径
    b_eq_ihs = IHS["size"] / np.sqrt(np.pi)
    r_spread_ihs = r_spreading_disk(a_eq_die, b_eq_ihs, IHS["thickness"],
                                    IHS["k"], 3000.0)
    r_ihs = r_cond_plane(IHS["thickness"], IHS["k"], a_ihs)
    r_tim2 = r_interface(tim2, a_ihs)
    b_eq_base = SINK["base"] / np.sqrt(np.pi)
    r_spread_base = r_spreading_disk(b_eq_ihs, b_eq_base, SINK["base_t"],
                                     k_sink, 800.0)
    r_base = r_cond_plane(SINK["base_t"], k_sink, a_base)

    mdot, d_h, gap, prop = air_mass_flow(velocity, 0.5 * (t_in + 55.0))
    h_fin = h_duct(velocity, d_h, 0.5 * (t_in + 55.0))
    r_fin = r_fin_array(h_fin, k_sink, SINK["fin_t"], fin_h, SINK["base"],
                        n_fins, a_base)
    # 空気の昇温。平均空気温度を代表点にするので 1/(2 mdot cp)
    r_caloric = 1.0 / (2.0 * mdot * prop["cp"])

    if verbose:
        print(f"  フィン間隙 {gap * 1000:.2f} mm, 水力直径 {d_h * 1000:.2f} mm, "
              f"h_fin = {h_fin:.1f} W/(m^2 K)")
        print(f"  質量流量 {mdot * 1000:.2f} g/s "
              f"(体積流量 {mdot / prop['rho'] * 60000:.1f} L/min), "
              f"空気の昇温 {power / (mdot * prop['cp']):.1f} K")
        from thermalnet import fin_efficiency
        print(f"  フィン効率 eta = "
              f"{100 * fin_efficiency(h_fin, k_sink, SINK['fin_t'], fin_h):.1f} %"
              f"、フィン対流抵抗 {r_fin:.4f} K/W")

    net = ThermalNetwork("cpu_cooling")
    caps = {
        "die": capacitance(2330.0, 700.0, a_die * DIE["thickness"]),
        "ihs": capacitance(8933.0, 385.0, a_ihs * IHS["thickness"]),
        "sink": capacitance(2700.0, 900.0,
                            a_base * SINK["base_t"]
                            + n_fins * SINK["fin_t"] * fin_h * SINK["base"]),
    }
    net.add_node("die", heat=power,
                 capacitance=caps["die"] if with_capacitance else 0.0)
    net.add_node("die_bottom")
    net.add_node("ihs_top", capacitance=caps["ihs"] if with_capacitance else 0.0)
    net.add_node("ihs_bottom")
    net.add_node("sink_base",
                 capacitance=caps["sink"] if with_capacitance else 0.0)
    net.add_node("fin_base")
    net.add_node("air_mean")
    net.add_node("air_in", fixed=t_in)

    net.add_resistor("die", "die_bottom", r_die, "die conduction")
    net.add_resistor("die_bottom", "ihs_top", r_tim1, "TIM1")
    net.add_resistor("ihs_top", "ihs_bottom", r_spread_ihs + r_ihs, "IHS + spreading")
    net.add_resistor("ihs_bottom", "sink_base", r_tim2, "TIM2")
    net.add_resistor("sink_base", "fin_base", r_spread_base + r_base,
                     "sink base + spreading")
    net.add_resistor("fin_base", "air_mean", r_fin, "fin convection")
    net.add_resistor("air_mean", "air_in", r_caloric, "air caloric")
    return net


# ---------------------------------------------------------------- 定常
section("(1) 定常設計点")

net = build_network(verbose=True)
sol = net.solve_steady(guess=T_IN + 30.0)
print()
print(sol.report(f"{POWER:.0f} W / 吸気 {T_IN:.0f} degC の定常解"))

r_total = (sol["die"] - T_IN) / POWER
print(f"\n  ジャンクション温度 T_j = {sol['die']:.1f} degC")
print(f"  合成熱抵抗 R_ja = {r_total:.4f} K/W")
limit = 100.0
print(f"  設計上限 {limit:.0f} degC に対する余裕 = {limit - sol['die']:+.1f} K")

section("(2) 抵抗の内訳 --- どこを叩くべきか")

labels = ["die conduction", "TIM1", "IHS + spreading", "TIM2",
          "sink base + spreading", "fin convection", "air caloric"]
res = []
for lab in labels:
    br = next(b for b in net.branches if b.label == lab)
    res.append(br.resistance)
res = np.array(res)
print("     区間                        R [K/W]    dT [K]   全体に占める割合")
for lab, r in zip(labels, res):
    print(f"    {lab:<26s} {r:8.4f}  {r * POWER:7.2f}   {100 * r / res.sum():6.1f} %")
print(f"    {'合計':<26s} {res.sum():8.4f}  {res.sum() * POWER:7.2f}")
print("""
  典型的にはフィン対流が支配的で、次に空気の昇温(caloric)。
  TIM や拡がり抵抗は小さく見えても、フィンを強化していくと相対的に効いてくる。
  「一番大きい抵抗から順に叩く」--- これが熱設計の鉄則。""")

section("(3) 感度分析 --- 各パラメータを 20 % 良くしたら何 K 下がるか")

base_tj = sol["die"]
cases = {
    "風速 +20 %": dict(velocity=AIR["velocity"] * 1.2),
    "フィン高さ +20 %": dict(fin_h=SINK["fin_h"] * 1.2),
    "フィン枚数 +20 %": dict(n_fins=int(SINK["n_fins"] * 1.2)),
    "TIM2 熱抵抗 -50 %": dict(tim2=TIM2 * 0.5),
    "ベース材を銅に (k=390)": dict(k_sink=390.0),
    "吸気温度 -5 K": dict(t_in=T_IN - 5.0),
}
deltas = {}
for name, kwargs in cases.items():
    s = build_network(**kwargs).solve_steady(guess=T_IN + 30.0)
    deltas[name] = s["die"] - base_tj
    print(f"    {name:<24s} T_j {s['die']:6.1f} degC  ({deltas[name]:+5.2f} K)")
print("""
  こういう表を 1 秒で作れるのが 1D モデルの価値。CFD で同じことをすると
  1 ケース数十分〜数時間かかる。CFD は「モデルの校正」に使い、
  探索は 1D で回す、という分業が現実的。""")

section("(4) 設計スイープ: 風速とフィン枚数")

velocities = np.linspace(0.5, 6.0, 20)
tj_v = [build_network(velocity=v).solve_steady(guess=T_IN + 30.0)["die"]
        for v in velocities]
fin_counts = np.arange(6, 45, 2)
tj_n = [build_network(n_fins=int(n), fixed_pressure=True)
        .solve_steady(guess=T_IN + 30.0)["die"] for n in fin_counts]
best = int(fin_counts[int(np.argmin(tj_n))])
print(f"  風速 0.5 -> 6.0 m/s で T_j は {tj_v[0]:.0f} -> {tj_v[-1]:.0f} degC")
print(f"  フィン枚数(ファン圧力一定の条件)の最適値 = {best} 枚 "
      f"(T_j = {min(tj_n):.1f} degC、流速 "
      f"{channel_velocity(best, SINK['fin_h'], True):.2f} m/s)")
print("""  フィンを増やすと表面積は増えるが、隙間が狭まって流速が落ちる。
  この綱引きから最適枚数が決まる。逆に「流速一定」で解くと
  枚数は多いほど良いという誤った結論になる --- 1D 化では
  ファンと流路の連成(システム抵抗曲線)を意識すること。
  この「トレードオフの山」を素早く見つけるのが 1D モデルの得意技。""")

section("(5) 過渡: ターボブースト 150 W を 20 秒")

net_t = build_network(with_capacitance=True)
net_t.set_heat("die", lambda t: 150.0 if 20.0 <= t < 40.0 else POWER)
# 定常解を初期条件にする(「定常運転中に急にブーストがかかる」状況)
tr = net_t.solve_transient(t_end=300.0, dt=0.02,
                           initial=dict(sol.temperatures), store_every=10)
i_peak = int(np.argmax(tr["die"]))
print(f"  定常 95 W  : T_j = {tr['die'][0]:.1f} degC")
print(f"  ピーク     : T_j = {tr['die'][i_peak]:.1f} degC (t = {tr.time[i_peak]:.0f} s)")
print(f"  150 W を定常で流したら {T_IN + 150 * r_total:.1f} degC まで上がる")
print(f"  -> 20 秒のバーストなら {T_IN + 150 * r_total - tr['die'][i_peak]:.1f} K"
      " 分を熱容量が肩代わりしている")
print("  ブースト時間の設計(何秒許すか)はこの曲線から決める。")

# ---------------------------------------------------------------- 図
fig, axgrid = plt.subplots(2, 3, figsize=(16, 8.2))
axes = axgrid.ravel()

axes[0].barh(labels, res, color="#4C72B0")
axes[0].set_xlabel("thermal resistance [K/W]")
axes[0].set_title("Resistance breakdown")
axes[0].grid(alpha=0.3, axis="x")

names = list(deltas)
axes[1].barh(names, [deltas[n] for n in names], color="#55A868")
axes[1].axvline(0, color="k", lw=0.8)
axes[1].set_xlabel("change in T_j [K]")
axes[1].set_yticks(range(len(names)),
                   ["velocity +20%", "fin height +20%", "fin count +20%",
                    "TIM2 -50%", "copper base", "inlet -5K"])
axes[1].set_title("Sensitivity (lower is better)")
axes[1].grid(alpha=0.3, axis="x")

axes[2].plot(velocities, tj_v, "o-", lw=2, color="#C44E52")
axes[2].axhline(limit, ls="--", color="gray")
axes[2].text(3.0, limit + 1.5, "design limit", color="gray")
axes[2].set_xlabel("fin channel velocity [m/s]")
axes[2].set_ylabel("junction temperature [degC]")
axes[2].set_title("Airflow sweep (velocity prescribed)")
axes[2].grid(alpha=0.3)

axes[3].plot(fin_counts, tj_n, "o-", lw=2, color="#8172B3")
axes[3].plot([best], [min(tj_n)], "r*", ms=16, label=f"optimum {best} fins")
axes[3].set_xlabel("number of fins")
axes[3].set_ylabel("junction temperature [degC]")
axes[3].set_ylim(min(tj_n) - 5, min(tj_n) + 60)   # 目盛りを最適点付近に寄せる
axes[3].set_title("Fin count with fan curve coupling")
axes[3].legend()
axes[3].grid(alpha=0.3)

axes[4].plot(tr.time, tr["die"], lw=2, label="junction")
axes[4].plot(tr.time, tr["sink_base"], lw=2, label="sink base")
axes[4].plot(tr.time, tr["air_mean"], lw=1.5, label="air (mean)")
axes[4].axvspan(20, 40, color="orange", alpha=0.25, label="150 W boost")
axes[4].set_xlabel("time [s]")
axes[4].set_ylabel("temperature [degC]")
axes[4].set_title("Turbo boost transient")
axes[4].legend(fontsize=8)
axes[4].grid(alpha=0.3)

# 6 枚目: 温度の「はしご」を可視化(熱回路の読み方)
chain = ["die", "die_bottom", "ihs_top", "ihs_bottom", "sink_base",
         "fin_base", "air_mean", "air_in"]
axes[5].step(range(len(chain)), [sol[n] for n in chain], where="post", lw=2,
             color="#4C72B0")
axes[5].set_xticks(range(len(chain)), chain, rotation=45, ha="right", fontsize=8)
axes[5].set_ylabel("temperature [degC]")
axes[5].set_title("Temperature ladder along the heat path")
axes[5].grid(alpha=0.3)

savefig(fig, "lesson08_heatsink.png")

print("""
まとめ(1D 熱回路網の実務的な使い方)
  * 経路を素直に直列/並列で書き下し、抵抗の内訳で優先順位を決める
  * 空気の昇温(caloric resistance)を忘れない --- 1D 化で最も多いミス
  * 拡がり抵抗を忘れない --- 小さい熱源ほど効く
  * 感度分析とスイープは 1D の独壇場。CFD は校正と最終確認に使う
  * 熱容量を入れれば、そのまま制御・保護ロジック用のモデルになる

  演習は exercises/exercises.md へ。
""")
