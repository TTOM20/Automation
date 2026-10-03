"""ロボットの寸法・脚配置などの設定。

座標系 (ボディ座標系):
    x: 前方, y: 左方, z: 上方  (単位: メートル, ラジアン)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

LEG_NAMES: tuple[str, ...] = ("FL", "FR", "RL", "RR")


@dataclass(frozen=True)
class LegGeometry:
    """3自由度脚のリンク長。

    hip_offset: 股関節(外転軸)から大腿の付け根までの横方向オフセット
    upper: 大腿リンク長
    lower: 下腿リンク長
    """

    hip_offset: float = 0.05
    upper: float = 0.12
    lower: float = 0.12


@dataclass(frozen=True)
class RobotConfig:
    """四足歩行ロボット全体の設定。"""

    body_length: float = 0.24
    body_width: float = 0.10
    leg: LegGeometry = field(default_factory=LegGeometry)
    stand_height: float = 0.17
    # 膝の向き: +1 で膝が後方に突き出す (犬型), -1 で前方に突き出す
    knee_direction: int = 1

    def side(self, leg: str) -> int:
        """左脚なら +1, 右脚なら -1。"""
        return 1 if leg[1] == "L" else -1

    def hip_position(self, leg: str) -> np.ndarray:
        """ボディ中心から見た各脚の股関節位置。"""
        fx = 1 if leg[0] == "F" else -1
        return np.array(
            [fx * self.body_length / 2, self.side(leg) * self.body_width / 2, 0.0]
        )

    def neutral_foot_position(self, leg: str) -> np.ndarray:
        """立位姿勢での足先位置 (ボディ座標系)。"""
        hip = self.hip_position(leg)
        return hip + np.array(
            [0.0, self.side(leg) * self.leg.hip_offset, -self.stand_height]
        )
