"""歩容 (gait) の生成。

各脚は周期 `period` で位相 0→1 を繰り返し、位相が `duty` 未満の間は
接地 (stance)、それ以降は遊脚 (swing) となる。脚ごとの位相差で
トロット・ウォークなどの歩容を表現する。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class GaitPattern:
    name: str
    duty: float  # 1周期のうち接地している割合
    offsets: dict[str, float] = field(default_factory=dict)  # 脚ごとの位相差


GAITS: dict[str, GaitPattern] = {
    "stand": GaitPattern("stand", 1.0, {"FL": 0.0, "FR": 0.0, "RL": 0.0, "RR": 0.0}),
    # 対角の脚を同時に動かす
    "trot": GaitPattern("trot", 0.5, {"FL": 0.0, "RR": 0.0, "FR": 0.5, "RL": 0.5}),
    # 常に3脚接地する静的に安定な歩容 (RL→FL→RR→FR の順に遊脚)
    "walk": GaitPattern("walk", 0.75, {"RL": 0.0, "FL": 0.25, "RR": 0.5, "FR": 0.75}),
    # 同じ側の脚を同時に動かす
    "pace": GaitPattern("pace", 0.5, {"FL": 0.0, "RL": 0.0, "FR": 0.5, "RR": 0.5}),
    # 前脚と後脚を交互に動かす
    "bound": GaitPattern("bound", 0.5, {"FL": 0.0, "FR": 0.0, "RL": 0.5, "RR": 0.5}),
}


@dataclass
class VelocityCommand:
    """ボディの目標速度 (ボディ座標系)。"""

    vx: float = 0.0  # 前進速度 [m/s]
    vy: float = 0.0  # 横移動速度 [m/s]
    yaw_rate: float = 0.0  # 旋回角速度 [rad/s]

    def foot_velocity(self, neutral: np.ndarray) -> np.ndarray:
        """中立足先位置 neutral におけるボディの水平速度。"""
        return np.array(
            [self.vx - self.yaw_rate * neutral[1], self.vy + self.yaw_rate * neutral[0], 0.0]
        )


class GaitGenerator:
    """位相に基づいて各脚の足先軌道を生成する。"""

    def __init__(
        self,
        pattern: GaitPattern | str = "trot",
        period: float = 0.5,
        step_height: float = 0.04,
    ) -> None:
        self.pattern = GAITS[pattern] if isinstance(pattern, str) else pattern
        if not 0.0 < self.pattern.duty <= 1.0:
            raise ValueError("duty は (0, 1] の範囲で指定してください")
        if period <= 0:
            raise ValueError("period は正の値で指定してください")
        self.period = period
        self.step_height = step_height
        self.time = 0.0

    def set_pattern(self, pattern: GaitPattern | str) -> None:
        self.pattern = GAITS[pattern] if isinstance(pattern, str) else pattern

    def advance(self, dt: float) -> None:
        self.time += dt

    def phase(self, leg: str) -> float:
        """脚の現在位相 [0, 1)。"""
        return (self.time / self.period + self.pattern.offsets[leg]) % 1.0

    def in_stance(self, leg: str) -> bool:
        return self.pattern.duty >= 1.0 or self.phase(leg) < self.pattern.duty

    def foot_position(
        self, leg: str, neutral: np.ndarray, cmd: VelocityCommand
    ) -> np.ndarray:
        """接地面基準フレームでの足先目標位置を返す。"""
        duty = self.pattern.duty
        stance_time = duty * self.period
        stride = cmd.foot_velocity(neutral) * stance_time
        if duty >= 1.0:
            return neutral.copy()

        phase = self.phase(leg)
        if phase < duty:
            # 接地中: 足先はボディに対して後方へ一定速度で移動
            s = phase / duty
            return neutral + stride * (0.5 - s)

        # 遊脚中: サイクロイド軌道で前方へ振り出す (始点・終点で速度0)
        s = (phase - duty) / (1.0 - duty)
        progress = s - math.sin(2 * math.pi * s) / (2 * math.pi)
        lift = self.step_height * (1 - math.cos(2 * math.pi * s)) / 2
        return neutral + stride * (progress - 0.5) + np.array([0.0, 0.0, lift])
