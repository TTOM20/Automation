"""レッスン共通の小道具(パス設定・見出し・図の保存)。"""

from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

FIGDIR = os.path.join(ROOT, "figures")

import matplotlib          # noqa: E402
matplotlib.use("Agg")      # ヘッドレス環境でも図を保存できるようにする
import matplotlib.pyplot as plt  # noqa: E402


def header(title: str) -> None:
    bar = "=" * 78
    print(f"\n{bar}\n {title}\n{bar}")


def section(title: str) -> None:
    print(f"\n--- {title} " + "-" * max(0, 70 - len(title)))


def savefig(fig, filename: str) -> str:
    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, filename)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"\n[図を保存] {os.path.relpath(path, ROOT)}")
    return path
