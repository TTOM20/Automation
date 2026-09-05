# モデルファイル書式リファレンス

`thermalnet.model` は YAML(または JSON)に書いた**物理量**から熱回路を組み立てます。
抵抗値を自分で計算する必要はありません。熱伝導率・寸法・発熱量・温度指定・流量を書けば、
対応する熱抵抗が自動で作られます。

```bash
python -m thermalnet.model models/cpu_heatsink.yaml
python -m thermalnet.model models/cpu_heatsink.yaml --set nodes.die.heat=125
python -m thermalnet.model models/cpu_heatsink.yaml --sweep nodes.die.heat=50:150:11 --watch die
python -m thermalnet.model models/power_module_liquid.yaml --transient --plot out.png
```

---

## ファイル全体の構造

```yaml
name: モデル名                # 任意
materials: {...}             # 任意。組み込み材料の追加・上書き
nodes: {...}                 # 必須。温度が決まる点
elements: [...]              # 必須。ノード間の熱抵抗
solve: {...}                 # 任意。定常解のオプション
transient: {...}             # 任意。過渡解析のオプション
```

---

## nodes --- 温度が決まる点

```yaml
nodes:
  die:
    heat: 95                 # [W] 発熱(電流源に相当)
    material: silicon        # 熱容量を材料から作る場合
    volume: 1.0e-7           #   volume [m^3] か mass [kg] のどちらか
  sink:
    capacitance: 300         # [J/K] 直接指定してもよい
  air_in:
    fixed: 30                # [degC] 温度指定境界(電圧源に相当)
  midpoint: {}               # 何も持たない中継ノード
```

| キー | 単位 | 意味 |
|---|---|---|
| `heat` | W | 発熱量。定常解で使われる |
| `heat_schedule` | — | `[[時刻, W], ...]` の階段状の発熱。過渡用。指定すると `heat` より優先 |
| `fixed` | degC | 温度指定境界。**最低 1 つ必要** |
| `capacitance` | J/K | 熱容量。過渡解析でのみ効く |
| `material` + `volume` | — | `C = rho cp V` を自動計算 |
| `material` + `mass` | — | `C = m cp` を自動計算 |

> **温度指定ノードが 1 つも無いとエラー**になります。熱の逃げ先(外気・冷却水)を必ず置いてください。
> 同じく、温度指定ノードにつながらない**浮きノード**もエラーで検出されます。

---

## elements --- ノード間の熱抵抗

共通のキー:

| キー | 意味 |
|---|---|
| `type` | 要素の種類(下表) |
| `from`, `to` | つなぐノード名(存在しなければ自動で作られる) |
| `name` | 表示用のラベル。`--set` で参照するときの名前にもなる |
| `disabled: true` | 一時的にこの要素を無効化(効果の切り分けに便利) |

面積は `area` [m²] のほか、`size`(正方形の一辺)、`width` + `length`、`diameter` でも書けます。
熱伝導率は `k` を直接書くか、`material: copper` のように材料名で指定します。

### type 一覧

| type | 必要な物理量 | 式 |
|---|---|---|
| `resistance` | `value` [K/W] | 直接指定 |
| `conduction` | `k` or `material`, `thickness`, 面積 | `L/(kA)` |
| `cylinder` | `k`, `r_in`, `r_out`, `length` | `ln(r2/r1)/(2πkL)` |
| `sphere` | `k`, `r_in`, `r_out` | `(1/r1-1/r2)/(4πk)` |
| `interface` | `resistivity` [m²K/W] + 面積、または `k`+`thickness`+面積 | TIM のカタログ値向け |
| `contact` | `conductance` [W/(m²K)], 面積 | `1/(h_c A)` |
| `convection` | `h` または `correlation`, 面積 | `1/(hA)` |
| `radiation` | `emissivity`, 面積, `view_factor`(既定 1) | `1/(h_r A)`、温度依存 |
| `fin_array` | `h`/`correlation`, `k`, `fin_thickness`, `fin_height`, `fin_width`, `n_fins`, `base_area` | フィン効率込み |
| `spreading` | `k`, `source_*`, `plate_*`, `thickness`, `h` | 拡がり抵抗(Lee et al.) |
| `flow` | `mass_flow` / `volume_flow` / `lpm` / `cfm`, `fluid` or `cp` | `1/(ṁcp)` |

### convection の相関式

`h` の代わりに `correlation` を書くと、**熱伝達率が温度の関数になり自動で非線形反復**されます。

| correlation | 必要なキー | 内容 |
|---|---|---|
| `natural_vertical` | `height` | 垂直平板の自然対流(Churchill–Chu) |
| `natural_horizontal_up` | `char_length` または 面積+`perimeter` | 上向き高温水平面 |
| `forced_plate` | `velocity`, `length` | 平板の強制対流 |
| `duct` | `velocity`, `hydraulic_diameter` | ダクト内強制対流 |

```yaml
- name: conv side
  type: convection
  from: case            # 面(表面温度)側を from に書く
  to: ambient           # 流体側を to に書く
  correlation: natural_vertical
  height: 0.05
  area: 0.035
```

### spreading の `h` は下流と整合させる

拡がり抵抗の相関式は「板の裏面から熱が抜ける速さ」を `h` として要求します。
この `h` は**その板より下流にある全抵抗**で決まります:

```
h_eq = 1 / (R_downstream × A_plate)
```

適当な値を入れると拡がり抵抗が数倍ずれます。**`h: auto` と書けば自動で反復して整合**させます:

```yaml
- name: substrate spreading
  type: spreading
  from: substrate_top
  to: substrate_bottom
  material: aln
  source_area: 1.0e-3      # source_size / source_radius / source_diameter でも可
  plate_area: 4.0e-3
  thickness: 0.00064
  h: auto                  # <- 下流の合成抵抗から自動決定(結果は要素内訳に表示)
```

### flow --- 流体への排熱(1D 化で最も忘れられる要素)

流路を通る流体は熱をもらって昇温します。これを入れないと冷却性能を過大評価します。

```yaml
- name: air caloric
  type: flow
  from: air_mean         # 流体の代表温度ノード
  to: air_in             # 入口(温度指定)
  fluid: air             # cp を材料表から引く。cp: 1007 と直接書いてもよい
  mass_flow: 0.00373     # [kg/s]。lpm / cfm / volume_flow でも書ける
  mean: true             # true なら 1/(2 ṁcp)、false なら 1/(ṁcp)
```

- `mean: true` … `from` のノードが**平均**流体温度を表す場合(既定)
- `mean: false` … `from` のノードが**出口**温度を表す場合

---

## 直列と並列の書き分け(最頻出のミス)

**同じ 2 ノード間に 2 本書くと並列合成されます。** 直列にしたいなら中間ノードを挟んでください。

```yaml
# NG: 拡がりと伝導が並列になってしまう
- {name: spreading,  type: spreading,  from: ihs_top, to: ihs_bottom, ...}
- {name: conduction, type: conduction, from: ihs_top, to: ihs_bottom, ...}

# OK: 中間ノードで直列にする
- {name: spreading,  type: spreading,  from: ihs_top, to: ihs_mid,    ...}
- {name: conduction, type: conduction, from: ihs_mid, to: ihs_bottom, ...}
```

逆に、対流と輻射のように**本当に並列**な経路は同じ 2 ノード間に並べて書きます。

---

## report --- レポート生成のオプション

```yaml
report:
  limit: 100         # 設計上限温度 [degC]。合否判定に使う
  target: die        # 評価対象ノード(省略時は最高温度のノードを自動選択)
```

```bash
python -m thermalnet.report models/cpu_heatsink.yaml --limit 100 -o report.html
```

CLI の `--limit` / `--target` / `--title` はこのブロックより優先されます。
発熱源のないモデル(壁など)では、投入熱量の代わりに境界を通過する熱量を基準に
合成熱抵抗を計算します(U 値の評価に使えます)。

## solve / transient

```yaml
solve:
  guess: 60          # 反復の初期推定温度 [degC]
  relaxation: 0.7    # 非線形反復の緩和係数(発散するときは下げる)
  source: die        # 合成熱抵抗を表示する起点(省略時は発熱ノードから自動)

transient:
  t_end: 300         # [s]
  dt: 0.05           # [s]
  initial: 30        # [degC]
  theta: 1.0         # 1.0 後退オイラー / 0.5 クランク・ニコルソン / 0.0 陽解法
```

---

## CLI

| オプション | 意味 |
|---|---|
| `--set PATH=VALUE` | 値を上書き(複数指定可)。例 `--set nodes.die.heat=120` |
| `--sweep PATH=LO:HI:N` | パラメータを振って表にする。`PATH=1,2,5` と列挙も可 |
| `--watch NODE` | スイープ・過渡で追跡するノード(複数指定可) |
| `--transient` | 過渡解析。`--t-end` `--dt` `--initial` で上書き可 |
| `--plot FILE.png` | 結果を図に保存 |
| `--quiet` | 要素内訳の表示を省略 |

`--set` のパスは `nodes.<ノード名>.<キー>` または `elements.<要素名>.<キー>` です。
要素名に空白がある場合は引用符で囲みます:

```bash
python -m thermalnet.model models/cpu_heatsink.yaml \
    --set "elements.fin convection.h=60" \
    --set "elements.air caloric.mass_flow=0.006"
```

---

## materials --- 材料を足す/上書きする

組み込みの材料表は `thermalnet/materials.py`(代表値。実設計ではデータシート値を使うこと)。
モデルファイル内で追加・上書きもできます:

```yaml
materials:
  my_tim:   {k: 5.0,  rho: 2500, cp: 1000}
  my_alloy: {k: 120,  rho: 2800, cp: 880}
```

---

## Python から使う

```python
from thermalnet.model import load_spec, build_network, apply_override, report_steady

spec = load_spec("models/cpu_heatsink.yaml")
apply_override(spec, "nodes.die.heat=125")
net, desc = build_network(spec)
sol = net.solve_steady(guess=60)
print(report_steady(net, desc, sol))
print(sol["die"], sol.flow_of("fin convection"))
```
