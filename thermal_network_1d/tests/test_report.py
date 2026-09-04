"""レポート生成(thermalnet.report / charts)の検証テスト。

実行:
    python tests/test_report.py
    pytest tests/test_report.py
"""

import glob
import io
import os
import re
import sys
import tempfile
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from thermalnet import charts
from thermalnet.model import build_network, load_spec
from thermalnet.report import (analyze, assumptions, main, main_heat_path,
                               render_html)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_DIR = os.path.join(ROOT, "models")

SIMPLE = {
    "name": "test model",
    "nodes": {"src": {"heat": 10.0, "capacitance": 5.0},
              "mid": {"capacitance": 50.0},
              "amb": {"fixed": 20.0}},
    "elements": [
        {"name": "upper", "type": "resistance", "from": "src", "to": "mid",
         "value": 2.0},
        {"name": "lower", "type": "resistance", "from": "mid", "to": "amb",
         "value": 3.0},
        {"name": "leak", "type": "resistance", "from": "src", "to": "amb",
         "value": 50.0},
    ],
}


# ------------------------------------------------------------------ 解析
def test_source_free_model_uses_throughput_heat():
    """発熱源が無い(温度差で駆動される)モデルでも解析できること。"""
    wall = load_spec(os.path.join(MODEL_DIR, "multilayer_wall.yaml"))
    a = analyze(wall)
    assert a["heat_label"] == "通過熱量"
    assert a["reference_heat"] > 0
    assert a["path_start"] == "room" and a["path"][-1] == "outdoor"
    assert a["target"] not in a["sinks"], "境界ノードを評価対象にしてはいけない"
    # U 値 = 1/(R A)。面積 10 m^2 で 0.3-0.5 W/(m^2 K) のオーダーになるはず
    u_value = 1.0 / (a["r_total"] * 10.0)
    assert 0.2 < u_value < 0.8, u_value


def test_analyze_basic_quantities():
    a = analyze(SIMPLE, limit=100.0)
    # 直列 2 + 3 と漏れ 50 の並列 -> 合成 4.717...
    assert abs(a["r_total"] - (a["temps"]["src"] - 20.0) / 10.0) < 1e-12
    assert a["target"] == "src"
    assert abs(a["residual"]) < 1e-9
    assert a["limit"] == 100.0


def test_main_heat_path_follows_largest_flow():
    """主熱流経路は熱流の大きい枝(漏れではなく本流)をたどること。"""
    net, _ = build_network(dict(SIMPLE))
    sol = net.solve_steady(guess=20.0)
    path, labels = main_heat_path(net, dict(sol.temperatures), "src", ["amb"])
    assert path == ["src", "mid", "amb"]
    assert labels == ["upper", "lower"]


def test_branches_are_classified_on_path_or_not():
    a = analyze(SIMPLE)
    on_path = {r["label"] for r in a["rows"] if r["on_path"]}
    assert on_path == {"upper", "lower"}
    assert not next(r for r in a["rows"] if r["label"] == "leak")["on_path"]


def test_sensitivity_matches_linear_expectation():
    """純直列の線形回路なら、区間を -20 % した効果は厳密に -0.2 R Q。"""
    series_only = {
        "nodes": {"src": {"heat": 10.0}, "mid": {}, "amb": {"fixed": 20.0}},
        "elements": [
            {"name": "upper", "type": "resistance", "from": "src", "to": "mid",
             "value": 2.0},
            {"name": "lower", "type": "resistance", "from": "mid", "to": "amb",
             "value": 3.0}],
    }
    a = analyze(series_only)
    upper = next(r for r in a["rows"] if r["label"] == "upper")
    delta = dict(a["sensitivity"])["upper を -20 %"]
    assert abs(delta - (-0.2 * upper["R"] * upper["Q"])) < 1e-9, delta
    # 発熱 +20 % は温度上昇が 20 % 増えること
    cond = dict(a["conditions"])["src の発熱 +20 %"]
    assert abs(cond - 0.2 * (a["temps"]["src"] - 20.0)) < 1e-9


def test_parallel_branch_softens_the_first_order_estimate():
    """並列に逃げ道があると、-0.2 R Q の一次近似より効果は小さくなる。

    抵抗を下げると分流が変わるため。内訳表から暗算した効果を
    そのまま設計判断に使うと過大評価になる、という実務的な注意点。
    """
    a = analyze(SIMPLE)
    upper = next(r for r in a["rows"] if r["label"] == "upper")
    naive = -0.2 * upper["R"] * upper["Q"]
    actual = dict(a["sensitivity"])["upper を -20 %"]
    assert naive < actual < 0, (naive, actual)


def test_condition_sensitivity_includes_boundary_shift():
    """境界温度を +5 K すると、線形回路では全ノードが +5 K 動く。"""
    a = analyze(SIMPLE)
    assert abs(dict(a["conditions"])["amb の温度 +5 K"] - 5.0) < 1e-9


def test_transient_is_included_when_capacitance_and_options_exist():
    spec = dict(SIMPLE, transient={"t_end": 200.0, "dt": 0.5, "initial": 20.0})
    a = analyze(spec)
    assert a["transient"] is not None
    assert 1 <= len(a["transient"]["series"]) <= 3
    # 熱容量が無ければ過渡は出さない
    no_cap = {"name": "x",
              "nodes": {"a": {"heat": 1.0}, "b": {"fixed": 0.0}},
              "elements": [{"type": "resistance", "from": "a", "to": "b",
                            "value": 1.0}],
              "transient": {"t_end": 10.0, "dt": 0.1}}
    assert analyze(no_cap)["transient"] is None


def test_assumptions_reflect_the_model_contents():
    """使った要素の種類に応じて、前提と適用限界の項目が変わること。"""
    plain = [name for name, _ in assumptions(analyze(SIMPLE))]
    assert "1 次元化" in plain
    assert "拡がり抵抗" not in plain

    spec = load_spec(os.path.join(MODEL_DIR, "cpu_heatsink.yaml"))
    rich = [name for name, _ in assumptions(analyze(spec))]
    for expected in ("拡がり抵抗", "流体の昇温", "界面熱抵抗"):
        assert expected in rich, f"{expected} が前提に出ていない"


def test_verdict_reacts_to_the_limit():
    html_ok = render_html(analyze(SIMPLE, limit=200.0), "t")
    html_ng = render_html(analyze(SIMPLE, limit=30.0), "t")
    assert "verdict-good" in html_ok
    assert "verdict-critical" in html_ng


# ------------------------------------------------------------------ 出力
def test_render_html_is_well_formed_enough():
    html = render_html(analyze(SIMPLE, limit=100.0), "テストレポート")
    assert html.count("<svg") == html.count("</svg>")
    assert html.count("<table") == html.count("</table>")
    assert html.count("<section") == html.count("</section>")
    assert "<title>テストレポート</title>" in html
    # アーティファクトのラッパが付けるので、自前の doctype/html/body は出さない
    for forbidden in ("<!doctype", "<html", "<body"):
        assert forbidden not in html.lower()


def test_charts_stay_inside_their_viewbox():
    """SVG の描画要素が viewBox の外にはみ出していないこと(ラベル切れの検出)。"""
    svgs = [
        charts.bar_chart([("very long element name here", 0.275),
                          ("short", 0.01)], "K/W"),
        charts.diverging_bar([("a を -20 %", -5.2), ("b の発熱 +20 %", 12.2)], "K"),
        charts.temperature_ladder(["die", "mid", "air"], [90.0, 70.0, 30.0],
                                  ["cond", "conv"]),
        charts.line_chart([0.0, 1.0, 2.0], [("die", [30.0, 60.0, 70.0])]),
    ]
    for svg in svgs:
        vb = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
        assert vb, "viewBox がない"
        width, height = float(vb.group(1)), float(vb.group(2))
        for x_attr in re.findall(r'\b(?:x|cx|x1|x2)="(-?[\d.]+)"', svg):
            assert -1.0 <= float(x_attr) <= width + 1.0, f"x={x_attr} が範囲外"
        for y_attr in re.findall(r'\b(?:y|cy|y1|y2)="(-?[\d.]+)"', svg):
            assert -1.0 <= float(y_attr) <= height + 1.0, f"y={y_attr} が範囲外"


def test_charts_reference_theme_tokens_only():
    """図の色は CSS クラス経由(= テーマ追従)で、直書きの色が無いこと。"""
    svg = charts.bar_chart([("a", 1.0)], "K/W") + charts.line_chart(
        [0.0, 1.0], [("a", [1.0, 2.0])])
    assert "fill=\"#" not in svg and "stroke=\"#" not in svg


def test_every_chart_mark_has_a_tooltip():
    svg = charts.bar_chart([("a", 1.0), ("b", 2.0)], "K/W")
    assert svg.count("data-tip=") == svg.count("<title>") == 2


def test_cli_generates_reports_for_every_model():
    with tempfile.TemporaryDirectory() as tmp:
        for path in sorted(glob.glob(os.path.join(MODEL_DIR, "*.yaml"))):
            out = os.path.join(tmp, os.path.basename(path) + ".html")
            with redirect_stdout(io.StringIO()):
                assert main([path, "-o", out, "--limit", "125"]) == 0
            body = open(out, encoding="utf-8").read()
            assert len(body) > 10000
            assert "検算" in body and "前提と適用限界" in body
            assert "エネルギー保存" in body


# ------------------------------------------------------------------ runner
def _main():
    tests = [(name, obj) for name, obj in sorted(globals().items())
             if name.startswith("test_") and callable(obj)]
    failures = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception as exc:                       # noqa: BLE001
            failures += 1
            print(f"  FAIL  {name}: {exc}")
    print(f"\n{len(tests) - failures} / {len(tests)} 件が成功")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_main())
