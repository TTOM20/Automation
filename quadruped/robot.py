"""四足歩行ロボットのコントローラ。"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .body import BodyPose
from .config import LEG_NAMES, RobotConfig
from .gait import GaitGenerator, VelocityCommand
from .kinematics import inverse_kinematics, leg_joint_positions


@dataclass
class RobotState:
    time: float
    joint_angles: dict[str, np.ndarray]  # 脚名 -> [θ1, θ2, θ3]
    foot_positions: dict[str, np.ndarray]  # 接地面基準フレームでの足先位置
    stance: dict[str, bool]


class Quadruped:
    """歩容生成・姿勢制御・逆運動学をまとめたコントローラ。

    例:
        robot = Quadruped()
        robot.set_gait("trot")
        state = robot.update(0.01, VelocityCommand(vx=0.2))
    """

    def __init__(
        self,
        config: RobotConfig | None = None,
        gait: GaitGenerator | None = None,
    ) -> None:
        self.config = config or RobotConfig()
        self.gait = gait or GaitGenerator("stand")
        self.pose = BodyPose()
        self.command = VelocityCommand()
        self.state = self._compute()

    def set_gait(self, name: str) -> None:
        self.gait.set_pattern(name)

    def set_command(self, cmd: VelocityCommand) -> None:
        self.command = cmd

    def set_pose(self, pose: BodyPose) -> None:
        self.pose = pose

    def update(self, dt: float, cmd: VelocityCommand | None = None) -> RobotState:
        """時間を dt 進めて、新しい関節角を計算する。"""
        if cmd is not None:
            self.command = cmd
        self.gait.advance(dt)
        self.state = self._compute()
        return self.state

    def _compute(self) -> RobotState:
        cfg = self.config
        angles: dict[str, np.ndarray] = {}
        feet: dict[str, np.ndarray] = {}
        stance: dict[str, bool] = {}
        for leg in LEG_NAMES:
            neutral = cfg.neutral_foot_position(leg)
            foot_ref = self.gait.foot_position(leg, neutral, self.command)
            foot_body = self.pose.to_body_frame(foot_ref)
            foot_hip = foot_body - cfg.hip_position(leg)
            angles[leg] = inverse_kinematics(
                foot_hip, cfg.leg, cfg.side(leg), cfg.knee_direction
            )
            feet[leg] = foot_ref
            stance[leg] = self.gait.in_stance(leg)
        return RobotState(self.gait.time, angles, feet, stance)

    def joint_positions(self) -> dict[str, np.ndarray]:
        """描画用: 接地面基準フレームでの各脚の関節位置 (4x3)。"""
        cfg = self.config
        out = {}
        for leg, a in self.state.joint_angles.items():
            pts = leg_joint_positions(a, cfg.leg, cfg.side(leg)) + cfg.hip_position(leg)
            out[leg] = np.array([self.pose.to_reference_frame(p) for p in pts])
        return out

    def body_corners(self) -> np.ndarray:
        """描画用: 接地面基準フレームでのボディ四隅 (FL, FR, RR, RL)。"""
        cfg = self.config
        order = ("FL", "FR", "RR", "RL")
        return np.array([self.pose.to_reference_frame(cfg.hip_position(l)) for l in order])
