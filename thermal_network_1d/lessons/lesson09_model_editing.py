"""lesson09: モデルファイルを編集して回す --- 実務のワークフロー

lesson01-08 では抵抗値を Python で計算していたが、実務では
「熱伝導率・発熱量・温度指定・流体への排熱」を書き換えて何度も回す。
そのためのモデル定義層(thermalnet.model)の使い方を体験する。

学ぶこと
  1. YAML に物理量を書くだけで熱回路が組めること
  2. ファイルを書き換えずに値を上書きして条件を振る方法(--set / apply_override)
  3. パラメータスイープと過渡解析をコマンド 1 本で回すこと
  4. 拡がり抵抗の h を下流と自己整合させる仕組み(h: auto)
  5. モデルの作り間違い(浮きノード・境界なし)が自動検出されること

実行:  python lessons/lesson09_model_editing.py
"""

from _util import header, section, savefig, plt
import copy
import os

import numpy as np

from thermalnet.model import (apply_override, build_network, load_spec,
                              report_steady)

MODELS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "models")

header("lesson09  モデルファイルを編集して回す")

print("""
これまでは抵抗値を Python で直接計算していた。実務では条件を何十通りも振るので、
「物理量を書いたファイル」と「解くコード」を分けたほうが速い。

    models/cpu_heatsink.yaml   <- 熱伝導率・寸法・発熱量・流量を書く
    python -m thermalnet.model models/cpu_heatsink.yaml

コマンドラインからも Python からも同じモデルを扱える。
""")

# ---------------------------------------------------------------- 読み込み
section("(1) モデルファイルを読んで解く")

spec = load_spec(os.path.join(MODELS, "cpu_heatsink.yaml"))
print(f"  モデル名: {spec['name']}")
print(f"  ノード {len(spec['nodes'])} 個 / 要素 {len(spec['elements'])} 個")
net, descriptions = build_network(spec)
sol = net.solve_steady(guess=60.0)
print()
print(report_steady(net, {}, sol))

# ---------------------------------------------------------------- 上書き
section("(2) ファイルを書き換えずに条件を振る")

print("  apply_override(spec, 'nodes.die.heat=125') のように書くだけ。")
print("  コマンドラインなら --set nodes.die.heat=125 と同じ。\n")

print("     条件                                   T_die [degC]   変化")
base = sol["die"]
cases = [
    ("基準(95 W, TIM2 = 2.0e-5)", []),
    ("発熱を 125 W に", ["nodes.die.heat=125"]),
    ("TIM2 を高性能品に (1.0e-5)", ["elements.TIM2.resistivity=1.0e-5"]),
    ("ベースを銅に (k=390)", ["elements.sink base conduction.k=390",
                              "elements.sink base spreading.k=390"]),
    ("風量 1.5 倍", ["elements.air caloric.mass_flow=0.0056"]),
    ("フィンの h を 60 に", ["elements.fin convection.h=60"]),
    ("吸気を 40 degC に", ["nodes.air_in.fixed=40"]),
]
results = {}
for label, overrides in cases:
    local = copy.deepcopy(spec)
    for assignment in overrides:
        apply_override(local, assignment)
    net_i, _ = build_network(local)
    t_die = net_i.solve_steady(guess=60.0)["die"]
    results[label] = t_die
    print(f"    {label:<38s} {t_die:9.2f}   {t_die - base:+7.2f} K")

print("""
  こうして「どの設計変更が何 K 効くか」を並べるのが熱設計の第一歩。
  CFD で同じ表を作ると数日かかるが、1D なら 1 秒で出る。""")

# ---------------------------------------------------------------- スイープ
section("(3) パラメータスイープ")

flows = np.linspace(0.0015, 0.008, 12)
tj = []
for mdot in flows:
    local = copy.deepcopy(spec)
    apply_override(local, f"elements.air caloric.mass_flow={mdot}")
    net_i, _ = build_network(local)
    tj.append(net_i.solve_steady(guess=60.0)["die"])
print("     質量流量 [g/s]   T_die [degC]")
for mdot, t in zip(flows[::2], tj[::2]):
    print(f"      {1000 * mdot:12.2f}   {t:11.2f}")
print(f"""
  流量を {1000 * flows[0]:.1f} -> {1000 * flows[-1]:.1f} g/s に増やしても
  T_die は {tj[0]:.1f} -> {tj[-1]:.1f} degC までしか下がらない。
  空気の昇温(caloric)は流量で効くが、フィン対流抵抗は変わらないため。
  「ファンを強くしても頭打ちになる」現象が回路の上で読み取れる。
  コマンドラインなら:
      python -m thermalnet.model models/cpu_heatsink.yaml \\
          --sweep "elements.air caloric.mass_flow=0.0015:0.008:12" --watch die""")

# ---------------------------------------------------------------- 自己整合
section("(4) 拡がり抵抗の h を自動で整合させる (h: auto)")

print("""  拡がり抵抗の相関式は「板の裏面から熱が抜ける速さ h」を必要とする。
  この h は板より下流の全抵抗で決まるので、本来は手で反復して合わせる:

      h_eq = 1 / (R_downstream x A_plate)

  モデルファイルに h: auto と書くと、この反復を自動でやってくれる。""")

spec_pm = load_spec(os.path.join(MODELS, "power_module_liquid.yaml"))
net_pm, desc_pm = build_network(spec_pm)
sol_pm = net_pm.solve_steady(guess=70.0)
for label in ("substrate spreading", "base spreading"):
    print(f"    {label:<22s} {desc_pm[label]}")
print(f"\n  結果: T_chip = {sol_pm['chip']:.2f} degC, "
      f"R(chip->water) = {(sol_pm['chip'] - 45.0) / 300.0:.4f} K/W")

# h を決め打ちしたときとの比較
for h_guess in (1000.0, 3000.0, 30000.0):
    local = copy.deepcopy(spec_pm)
    apply_override(local, f"elements.substrate spreading.h={h_guess}")
    apply_override(local, f"elements.base spreading.h={h_guess}")
    net_g, _ = build_network(local)
    t_chip = net_g.solve_steady(guess=70.0)["chip"]
    print(f"    h を {h_guess:6.0f} と決め打ちすると T_chip = {t_chip:6.2f} degC"
          f"  ({t_chip - sol_pm['chip']:+6.2f} K のずれ)")
print("  h の当てずっぽうがそのまま数十 K の誤差になる。auto を使うのが安全。")

# ---------------------------------------------------------------- 検証
section("(5) モデルの作り間違いは自動検出される")

for label, broken in [
        ("温度指定ノードが無い",
         {"nodes": {"a": {"heat": 10.0}, "b": {}},
          "elements": [{"type": "resistance", "from": "a", "to": "b", "value": 1.0}]}),
        ("浮きノードがある",
         {"nodes": {"a": {"heat": 10.0}, "amb": {"fixed": 25.0}, "orphan": {}},
          "elements": [{"type": "resistance", "from": "a", "to": "amb",
                        "value": 1.0}]})]:
    try:
        build_network(broken)
        print(f"    {label:<22s} -> 検出できず(バグ)")
    except ValueError as exc:
        print(f"    {label:<22s} -> {str(exc).splitlines()[0]}")
print("""
  熱回路は「絵を描けてしまう」ぶん、つなぎ忘れに気づきにくい。
  解く前に構造を検査する癖をつけると事故が減る。""")

# ---------------------------------------------------------------- 図
fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))

labels_en = ["baseline", "125 W", "better TIM2", "copper base",
             "1.5x airflow", "h = 60", "inlet 40C"]
values = [results[label] for label, _ in cases]
colors = ["#4C72B0"] + ["#C44E52" if v > base else "#55A868" for v in values[1:]]
axes[0].barh(labels_en, values, color=colors)
axes[0].axvline(base, color="k", ls="--", lw=1)
axes[0].set_xlim(min(values) - 5, max(values) + 5)
axes[0].set_xlabel("junction temperature [degC]")
axes[0].set_title("Design variants (dashed = baseline)")
axes[0].grid(alpha=0.3, axis="x")

axes[1].plot(1000 * flows, tj, "o-", lw=2, color="#4C72B0")
axes[1].set_xlabel("air mass flow [g/s]")
axes[1].set_ylabel("junction temperature [degC]")
axes[1].set_title("Airflow sweep from the model file")
axes[1].grid(alpha=0.3)

savefig(fig, "lesson09_model_editing.png")

print("""
まとめ
  * 物理量は YAML に、解く手順はコードに。分けると条件出しが速い。
  * --set / apply_override で「ファイルを汚さずに」条件を振れる。
  * h: auto のように、整合が必要なパラメータは自動化しておくと事故が減る。
  * 書式の全リファレンスは docs/model_format.md。

  次は自分の装置をモデル化してみること。手順は:
    1. 熱の経路を紙に描く(どこからどこへ流れるか)
    2. 経路の切れ目にノードを置き、区間ごとに type を選ぶ
    3. 温度指定境界(外気・冷却水)を必ず 1 つ置く
    4. 解いて、抵抗の内訳の大きい順に対策を考える
""")
