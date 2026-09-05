"""レポート用のインライン SVG チャート生成。

外部ライブラリに依存せず、テーマ(明/暗)に追従する SVG を文字列で組み立てる。
色は CSS カスタムプロパティ経由で参照するので、明暗の切り替えは CSS 側だけで済む。

配色は検証済みパレット(カテゴリ 3 色 + 発散 2 色 + ステータス 4 色)に従う。
系列が 3 本を超える場合は「その他」にまとめるか分割すること。
"""

from __future__ import annotations

import html
from typing import Iterable, List, Sequence, Tuple

# --- レイアウト定数 ---------------------------------------------------
BAR_H = 22           # 棒 1 本の高さ [px]
BAR_GAP = 10         # 棒の間隔(2px 以上の面ギャップを確保)
LABEL_W = 200        # 左のラベル幅
VALUE_W = 108        # 右の数値幅
PAD = 12


def _esc(text) -> str:
    return html.escape(str(text), quote=True)


def _fmt(value: float, digits: int = 3) -> str:
    if value == 0:
        return "0"
    if abs(value) >= 1000 or abs(value) < 1e-3:
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def _ticks(lo: float, hi: float, count: int = 5) -> List[float]:
    """きりのよい目盛りを返す。"""
    if hi <= lo:
        return [lo]
    raw = (hi - lo) / count
    magnitude = 10 ** int(__import__("math").floor(__import__("math").log10(raw)))
    for mult in (1, 2, 2.5, 5, 10):
        step = magnitude * mult
        if raw <= step:
            break
    start = step * int(lo / step) - (step if lo < 0 else 0)
    values, current = [], start
    while current <= hi + step * 0.5:
        if current >= lo - step * 0.5:
            values.append(round(current, 10))
        current += step
    return values


# =====================================================================
# 横棒グラフ(単一系列 = 大きさの比較)
# =====================================================================
def bar_chart(items: Sequence[Tuple[str, float]], unit: str,
              highlight: Iterable[str] = (), value_digits: int = 4) -> str:
    """items = [(ラベル, 値), ...] の横棒グラフ。値は正の量を想定。

    highlight に入れたラベルは主熱流経路として濃い色で描く。
    """
    if not items:
        return ""
    highlight = set(highlight)
    max_value = max(v for _, v in items) or 1.0
    plot_w = 420
    height = PAD * 2 + len(items) * (BAR_H + BAR_GAP) - BAR_GAP + 42
    width = LABEL_W + plot_w + VALUE_W

    parts = [f'<svg class="chart" viewBox="0 0 {width} {height}" '
             f'role="img" aria-label="熱抵抗の内訳">']
    # 目盛り
    for tick in _ticks(0.0, max_value):
        x = LABEL_W + plot_w * tick / max_value
        parts.append(f'<line x1="{x:.1f}" y1="{PAD}" x2="{x:.1f}" '
                     f'y2="{height - 42:.1f}" class="grid" />')
        parts.append(f'<text x="{x:.1f}" y="{height - 26}" class="tick" '
                     f'text-anchor="middle">{_fmt(tick, 3)}</text>')
    parts.append(f'<text x="{LABEL_W + plot_w / 2:.0f}" y="{height - 10}" '
                 f'class="axis-title" text-anchor="middle">{_esc(unit)}</text>')

    for i, (label, value) in enumerate(items):
        y = PAD + i * (BAR_H + BAR_GAP)
        w = max(plot_w * value / max_value, 1.5)
        cls = "bar bar-strong" if label in highlight else "bar"
        tip = f"{label}: {_fmt(value, value_digits)} {unit}"
        parts.append(
            f'<rect x="{LABEL_W}" y="{y}" width="{w:.1f}" height="{BAR_H}" '
            f'rx="4" class="{cls}" data-tip="{_esc(tip)}"><title>{_esc(tip)}'
            f'</title></rect>')
        parts.append(f'<text x="{LABEL_W - 8}" y="{y + BAR_H * 0.72:.1f}" '
                     f'class="bar-label" text-anchor="end">{_esc(label)}</text>')
        parts.append(f'<text x="{LABEL_W + plot_w + 8}" '
                     f'y="{y + BAR_H * 0.72:.1f}" class="bar-value">'
                     f'{_fmt(value, value_digits)}</text>')
    parts.append('</svg>')
    return "\n".join(parts)


# =====================================================================
# 発散棒グラフ(符号のある変化量 = 感度分析)
# =====================================================================
def diverging_bar(items: Sequence[Tuple[str, float]], unit: str) -> str:
    """items = [(ラベル, 変化量), ...]。負(改善)は青、正(悪化)は赤。"""
    if not items:
        return ""
    span = max(abs(v) for _, v in items) or 1.0
    plot_w = 420
    zero_x = LABEL_W + plot_w / 2
    height = PAD * 2 + len(items) * (BAR_H + BAR_GAP) - BAR_GAP + 42
    width = LABEL_W + plot_w + VALUE_W

    parts = [f'<svg class="chart" viewBox="0 0 {width} {height}" '
             f'role="img" aria-label="感度分析">']
    for tick in _ticks(-span, span, 4):
        x = zero_x + (plot_w / 2) * tick / span
        parts.append(f'<line x1="{x:.1f}" y1="{PAD}" x2="{x:.1f}" '
                     f'y2="{height - 42:.1f}" class="grid" />')
        parts.append(f'<text x="{x:.1f}" y="{height - 26}" class="tick" '
                     f'text-anchor="middle">{tick:+.1f}</text>')
    parts.append(f'<line x1="{zero_x:.1f}" y1="{PAD}" x2="{zero_x:.1f}" '
                 f'y2="{height - 42:.1f}" class="zero-line" />')
    parts.append(f'<text x="{zero_x:.0f}" y="{height - 10}" class="axis-title" '
                 f'text-anchor="middle">改善 &#8592; {_esc(unit)} &#8594; 悪化</text>')

    for i, (label, value) in enumerate(items):
        y = PAD + i * (BAR_H + BAR_GAP)
        w = max(abs(value) / span * (plot_w / 2), 1.5)
        x = zero_x if value >= 0 else zero_x - w
        cls = "bar-worse" if value > 0 else "bar-better"
        sign = "悪化" if value > 0 else "改善"
        tip = f"{label}: {value:+.2f} {unit}({sign})"
        parts.append(
            f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{BAR_H}" '
            f'rx="4" class="{cls}" data-tip="{_esc(tip)}"><title>{_esc(tip)}'
            f'</title></rect>')
        parts.append(f'<text x="{LABEL_W - 8}" y="{y + BAR_H * 0.72:.1f}" '
                     f'class="bar-label" text-anchor="end">{_esc(label)}</text>')
        parts.append(f'<text x="{LABEL_W + plot_w + 8}" '
                     f'y="{y + BAR_H * 0.72:.1f}" class="bar-value">'
                     f'{value:+.2f}</text>')
    parts.append('</svg>')
    return "\n".join(parts)


# =====================================================================
# 温度ラダー(主熱流経路に沿った温度の階段)
# =====================================================================
def temperature_ladder(path: Sequence[str], temps: Sequence[float],
                       labels: Sequence[str]) -> str:
    """主熱流経路に沿った温度降下を階段で描く。

    区間名は長くなりがちで図の上では重なるため、図には通し番号だけを置き、
    番号と区間名の対応は本文の表(表 2)で引く。
    """
    if len(path) < 2:
        return ""
    left, right, top, bottom = 62, 20, 20, 86
    plot_w = max(560, 74 * (len(path) - 1))
    width = left + plot_w + right
    height = 320
    plot_h = height - top - bottom
    t_max, t_min = max(temps), min(temps)
    span = (t_max - t_min) or 1.0
    t_hi, t_lo = t_max + span * 0.10, t_min - span * 0.06

    def y_of(t):
        return top + plot_h * (t_hi - t) / (t_hi - t_lo)

    step_w = plot_w / (len(path) - 1)
    parts = [f'<svg class="chart" viewBox="0 0 {width} {height}" '
             f'role="img" aria-label="主熱流経路に沿った温度">']
    for tick in _ticks(t_lo, t_hi, 5):
        y = y_of(tick)
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" '
                     f'y2="{y:.1f}" class="grid" />')
        parts.append(f'<text x="{left - 8}" y="{y + 4:.1f}" class="tick" '
                     f'text-anchor="end">{tick:.0f}</text>')
    parts.append(f'<text x="18" y="{top + plot_h / 2:.0f}" class="axis-title" '
                 f'transform="rotate(-90 18 {top + plot_h / 2:.0f})" '
                 f'text-anchor="middle">温度 [degC]</text>')

    points = [(left + i * step_w, y_of(t)) for i, t in enumerate(temps)]
    # 階段線: 水平(踏面)= その区間の抵抗、垂直(蹴上げ)= そこで生じる温度降下
    d = [f"M {points[0][0]:.1f} {points[0][1]:.1f}"]
    for i in range(1, len(points)):
        d.append(f"L {points[i][0]:.1f} {points[i - 1][1]:.1f}")
        d.append(f"L {points[i][0]:.1f} {points[i][1]:.1f}")
    parts.append(f'<path d="{" ".join(d)}" class="line" />')

    # 区間番号を踏面の上に丸で置く(表 2 の # と対応)
    for i, label in enumerate(labels):
        if i + 1 >= len(points):
            break
        x = (points[i][0] + points[i + 1][0]) / 2
        y = points[i][1]
        drop = temps[i] - temps[i + 1]
        tip = f"{i + 1}. {label}: {drop:.2f} K の降下"
        parts.append(f'<circle cx="{x:.1f}" cy="{y + 15:.1f}" r="9" '
                     f'class="span-marker" data-tip="{_esc(tip)}">'
                     f'<title>{_esc(tip)}</title></circle>')
        parts.append(f'<text x="{x:.1f}" y="{y + 19:.1f}" class="span-no" '
                     f'text-anchor="middle">{i + 1}</text>')

    for i, ((x, y), name, t) in enumerate(zip(points, path, temps)):
        tip = f"{name}: {t:.2f} degC"
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" '
                     f'class="dot" data-tip="{_esc(tip)}"><title>{_esc(tip)}'
                     f'</title></circle>')
        parts.append(f'<text x="{x:.1f}" y="{y - 11:.1f}" class="point-label" '
                     f'text-anchor="middle">{t:.1f}</text>')
        label_y = height - bottom + 16
        parts.append(f'<text x="{x:.1f}" y="{label_y:.1f}" class="tick" '
                     f'text-anchor="end" transform="rotate(-40 {x:.1f} '
                     f'{label_y:.1f})">{_esc(name)}</text>')
    parts.append('</svg>')
    return "\n".join(parts)


# =====================================================================
# 折れ線(過渡応答。系列は最大 3 本)
# =====================================================================
def line_chart(time: Sequence[float], series: Sequence[Tuple[str, Sequence[float]]],
               y_label: str = "温度 [degC]",
               shade: Tuple[float, float] | None = None,
               shade_label: str = "") -> str:
    if not series:
        return ""
    width, height = 760, 300
    left, right, top, bottom = 62, 110, 18, 46
    plot_w = width - left - right
    plot_h = height - top - bottom
    t0, t1 = float(time[0]), float(time[-1])
    all_values = [v for _, values in series for v in values]
    v_min, v_max = min(all_values), max(all_values)
    span = (v_max - v_min) or 1.0
    v_hi, v_lo = v_max + span * 0.10, v_min - span * 0.10

    def x_of(t):
        return left + plot_w * (t - t0) / ((t1 - t0) or 1.0)

    def y_of(v):
        return top + plot_h * (v_hi - v) / (v_hi - v_lo)

    parts = [f'<svg class="chart" viewBox="0 0 {width} {height}" '
             f'role="img" aria-label="過渡応答">']
    if shade:
        x_a, x_b = x_of(shade[0]), x_of(shade[1])
        parts.append(f'<rect x="{x_a:.1f}" y="{top}" width="{max(x_b - x_a, 1):.1f}" '
                     f'height="{plot_h:.1f}" class="shade" />')
        parts.append(f'<text x="{(x_a + x_b) / 2:.1f}" y="{top + 14:.0f}" '
                     f'class="tick" text-anchor="middle">{_esc(shade_label)}</text>')
    for tick in _ticks(v_lo, v_hi, 5):
        y = y_of(tick)
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" '
                     f'y2="{y:.1f}" class="grid" />')
        parts.append(f'<text x="{left - 8}" y="{y + 4:.1f}" class="tick" '
                     f'text-anchor="end">{tick:.0f}</text>')
    for tick in _ticks(t0, t1, 6):
        x = x_of(tick)
        parts.append(f'<text x="{x:.1f}" y="{height - 22:.0f}" class="tick" '
                     f'text-anchor="middle">{tick:.0f}</text>')
    parts.append(f'<text x="{left + plot_w / 2:.0f}" y="{height - 6:.0f}" '
                 f'class="axis-title" text-anchor="middle">時間 [s]</text>')
    parts.append(f'<text x="18" y="{top + plot_h / 2:.0f}" class="axis-title" '
                 f'transform="rotate(-90 18 {top + plot_h / 2:.0f})" '
                 f'text-anchor="middle">{_esc(y_label)}</text>')

    for idx, (name, values) in enumerate(series[:3]):
        pts = " ".join(f"{x_of(t):.1f},{y_of(v):.1f}" for t, v in zip(time, values))
        parts.append(f'<polyline points="{pts}" class="line s{idx + 1}" />')
        # 直接ラベル(凡例と二重化して色だけに頼らない)
        parts.append(f'<text x="{left + plot_w + 8:.0f}" '
                     f'y="{y_of(values[-1]) + 4:.1f}" class="series-label s{idx + 1}">'
                     f'{_esc(name)}</text>')
        peak = max(range(len(values)), key=lambda i: values[i])
        tip = f"{name} 最高 {values[peak]:.2f} degC (t = {time[peak]:.1f} s)"
        parts.append(f'<circle cx="{x_of(time[peak]):.1f}" '
                     f'cy="{y_of(values[peak]):.1f}" r="4.5" class="dot s{idx + 1}" '
                     f'data-tip="{_esc(tip)}"><title>{_esc(tip)}</title></circle>')
    parts.append('</svg>')
    return "\n".join(parts)
