"""関節角からサーボのパルス幅への変換。

実機 (例: PCA9685 + ホビー用サーボ) を動かす際に、
キャリブレーション値を保持して角度 [rad] をパルス幅 [µs] に変換する。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol

from .config import LEG_NAMES


@dataclass
class ServoCalibration:
    channel: int
    center_us: float = 1500.0  # 関節角 0 rad のときのパルス幅
    us_per_rad: float = 2000.0 / math.pi  # 500-2500µs で 180度
    direction: int = 1  # 取り付け向きによる回転方向の反転
    offset_rad: float = 0.0  # 取り付け誤差の補正
    min_us: float = 500.0
    max_us: float = 2500.0

    def to_pulse(self, angle: float) -> float:
        us = self.center_us + self.direction * (angle + self.offset_rad) * self.us_per_rad
        return min(self.max_us, max(self.min_us, us))


class ServoDriver(Protocol):
    def set_pulse(self, channel: int, pulse_us: float) -> None: ...


class PrintDriver:
    """ハードウェアなしで動作確認するためのドライバ。"""

    def __init__(self) -> None:
        self.last: dict[int, float] = {}

    def set_pulse(self, channel: int, pulse_us: float) -> None:
        self.last[channel] = pulse_us


def default_calibration() -> dict[str, list[ServoCalibration]]:
    """脚ごとに [外転, 股, 膝] の3チャンネルを順に割り当てる。

    右脚は左脚と鏡像に取り付ける前提で回転方向を反転する。
    """
    cal = {}
    for i, leg in enumerate(LEG_NAMES):
        sign = 1 if leg[1] == "L" else -1
        cal[leg] = [ServoCalibration(channel=i * 3 + j, direction=sign) for j in range(3)]
    return cal


def apply_joint_angles(
    driver: ServoDriver,
    joint_angles: dict[str, "list[float]"],
    calibration: dict[str, list[ServoCalibration]],
) -> None:
    for leg, angles in joint_angles.items():
        for angle, cal in zip(angles, calibration[leg]):
            driver.set_pulse(cal.channel, cal.to_pulse(float(angle)))
