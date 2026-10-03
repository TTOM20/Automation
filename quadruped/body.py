"""ボディ姿勢 (並進・回転) の変換。"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


def rotation_matrix(roll: float, pitch: float, yaw: float) -> np.ndarray:
    """ZYX オイラー角 (yaw→pitch→roll) から回転行列を作る。"""
    cr, sr = math.cos(roll), math.sin(roll)
    cp, sp = math.cos(pitch), math.sin(pitch)
    cy, sy = math.cos(yaw), math.sin(yaw)
    rx = np.array([[1, 0, 0], [0, cr, -sr], [0, sr, cr]])
    ry = np.array([[cp, 0, sp], [0, 1, 0], [-sp, 0, cp]])
    rz = np.array([[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]])
    return rz @ ry @ rx


@dataclass
class BodyPose:
    """接地面基準フレームに対するボディの姿勢オフセット。"""

    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    roll: float = 0.0
    pitch: float = 0.0
    yaw: float = 0.0

    def to_body_frame(self, point: np.ndarray) -> np.ndarray:
        """基準フレーム上の点をボディ座標系に変換する。"""
        r = rotation_matrix(self.roll, self.pitch, self.yaw)
        return r.T @ (np.asarray(point, dtype=float) - np.array([self.x, self.y, self.z]))

    def to_reference_frame(self, point: np.ndarray) -> np.ndarray:
        """ボディ座標系の点を基準フレームに変換する。"""
        r = rotation_matrix(self.roll, self.pitch, self.yaw)
        return r @ np.asarray(point, dtype=float) + np.array([self.x, self.y, self.z])
