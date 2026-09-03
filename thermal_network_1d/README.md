# 1D 熱回路網(サーマルネットワーク)ハンズオン教材

伝熱 CFD 解析を **1D の熱回路網に落とし込む** ための、手を動かして学ぶ教材です。
Python で書かれた小さなソルバ(約 300 行)と、8 本のレッスン、6 問の演習、
検証テストで構成されています。

> **狙い**
> CFD は「答えを出す道具」ですが、設計を回すには遅すぎます。
> 一方 1D 熱回路網は一瞬で解けますが、モデルの立て方を誤ると平気で 20 % 外します。
> この教材では、**CFD 相当の多次元解を自分で解いて正解を作り**、
> それを 1D に還元する過程を全部自分の手で確かめます。

## セットアップ

```bash
pip install -r requirements.txt      # numpy, scipy, matplotlib
python run_all.py                    # 全レッスン + テストを実行(図は figures/ へ)
```

Python 3.10 以上。個別に走らせる場合:

```bash
cd lessons && python lesson01_thermal_resistance.py
```

## 学習ロードマップ

| # | レッスン | 学ぶこと | キーワード |
|---|---|---|---|
| 01 | `lesson01_thermal_resistance.py` | 熱抵抗の定義とオームの法則アナロジー | 伝導/対流/輻射、Bi 数 |
| 02 | `lesson02_series_parallel.py` | 直列・並列合成、多層壁、熱橋 | U 値、支配抵抗、並列枝 |
| 03 | `lesson03_nodal_method.py` | 節点法の行列を自分で組む | コンダクタンス行列、保存則 |
| 04 | `lesson04_discretization.py` | **CFD の離散化 = 熱回路**であること | 有限体積法、格子収束、半セル抵抗 |
| 05 | `lesson05_transient.py` | 過渡応答、時定数、陰解法 | RC、τ=RC、θ 法、安定限界 |
| 06 | `lesson06_nonlinear.py` | 温度依存 h の反復解法 | Churchill–Chu、輻射、緩和係数 |
| 07 | `lesson07_cfd_to_1d.py` | **CFD から 1D へ還元する型** | 拡がり抵抗、抵抗の抽出、適用範囲 |
| 08 | `lesson08_project_heatsink.py` | 総合演習:空冷 CPU のフルモデル | フィン効率、流体昇温、感度分析、ファン曲線 |

各レッスンは単体で完結しており、実行すると解説・数値・検算がターミナルに出て、
`figures/` に図が保存されます。**まずそのまま実行し、次に数値を書き換えて遊ぶ**のが
一番早い学び方です。

## この教材の核心(lesson04 と 07)

**lesson04**: 1D 伝導方程式を有限体積法で離散化すると、それは文字どおり
「抵抗のはしご回路」になります。CFD ソルバがやっていることと、
熱回路網ソルバがやっていることは同じです。違いは次元とセル数だけ。

**lesson07**: では何が違うのか。3D では熱が横方向にも広がるため、
小さな熱源では素朴な 1D 式 `R = t/(kA) + 1/(hA)` が **15 % ほど温度を過小評価** します。
この差の正体が **拡がり抵抗**。教材では軸対称 FVM ソルバ(= ミニ CFD)で正解を作り、
差分を抵抗として抽出し、1D 回路に戻して一致することを確認し、
さらに「その抵抗はどの条件まで流用してよいか」を数値で確かめます。

これが 1D 化の型です:

```
経路を切る → CFD から代表温度を面積加重平均で拾う
          → R = ΔT/Q で抵抗を作る → 回路を再構成して一致を確認
          → 適用範囲を明記する
```

## ディレクトリ構成

```
thermal_network_1d/
├── thermalnet/              学習用ミニライブラリ
│   ├── elements.py          熱抵抗の部品カタログ(伝導/対流/輻射/フィン/拡がり/流体)
│   ├── network.py           節点法ソルバ(定常・過渡・非線形)
│   ├── fluids.py            空気物性と対流相関式
│   └── fvm.py               軸対称有限体積ソルバ(= 答え合わせ用のミニ CFD)
├── lessons/                 レッスン 01-08(実行するだけで解説が出る)
├── exercises/
│   ├── exercises.md         演習 6 問(答え合わせ用の数値つき)
│   └── solutions/           解答スクリプト
├── tests/test_thermalnet.py 解析解・保存則との検証テスト(19 件)
├── docs/cheatsheet.md       公式・無次元数・よくあるミス一覧
└── figures/                 実行すると図が出力される
```

## ライブラリの使い方(最小例)

```python
from thermalnet import ThermalNetwork

net = ThermalNetwork("cpu")
net.add_node("die", heat=95.0, capacitance=0.11)   # 発熱 [W], 熱容量 [J/K]
net.add_node("sink", capacitance=300.0)
net.add_node("air", fixed=30.0)                    # 温度固定境界 [degC]
net.add_resistor("die", "sink", 0.15, "TIM")       # 熱抵抗 [K/W]
net.add_resistor("sink", "air", 0.35, "heatsink")

sol = net.solve_steady()
print(sol["die"], sol.flow_of("TIM"))
print(sol.report())

tr = net.solve_transient(t_end=600.0, dt=0.1, initial=30.0)   # 過渡
```

温度依存の抵抗(自然対流・輻射)は関数として渡すと自動で反復されます:

```python
net.add_resistor("case", "air",
                 lambda ts, ta: r_conv(h_natural_vertical_plate(ts, ta, 0.05), area))
```

## 検証について

`tests/test_thermalnet.py` は 19 件のテストで、以下を確認しています:

- 平板/円筒伝導、輻射線形化、フィン効率、拡がり抵抗の解析解との一致
- 直列・並列合成、分流比、エネルギー保存
- 1D はしご回路の**収束次数 2** と表面熱流束の厳密性
- 集中定数の過渡解と指数関数解の一致、C–N が後退オイラーより高精度なこと
- 陽解法が `Δt > 2RC` で発散すること
- 軸対称 FVM の格子収束と、`R_total = R_1d + R_spread` の分解
- 拡がり抵抗の相関式(Lee et al.)が FVM 解と 10 % 以内で一致すること

```bash
python tests/test_thermalnet.py     # pytest なしでも動く
```

## 次に読むもの

`docs/cheatsheet.md` に公式・無次元数・**1D 化でやりがちなミス一覧**をまとめてあります。
実務で 1D モデルを作るときは、そのミス一覧をチェックリストとして使ってください。
