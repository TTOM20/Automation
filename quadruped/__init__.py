"""四足歩行ロボットの運動学・歩容生成ライブラリ。"""

from .body import BodyPose
from .config import LEG_NAMES, LegGeometry, RobotConfig
from .gait import GAITS, GaitGenerator, GaitPattern, VelocityCommand
from .kinematics import IKError, forward_kinematics, inverse_kinematics
from .robot import Quadruped, RobotState

__all__ = [
    "BodyPose",
    "GAITS",
    "GaitGenerator",
    "GaitPattern",
    "IKError",
    "LEG_NAMES",
    "LegGeometry",
    "Quadruped",
    "RobotConfig",
    "RobotState",
    "VelocityCommand",
    "forward_kinematics",
    "inverse_kinematics",
]
