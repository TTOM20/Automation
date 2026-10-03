"""3自由度脚の順運動学・逆運動学。

関節角:
    abduction (θ1): x軸まわりの股関節外転
    hip       (θ2): 外転後の y軸まわりの股関節ピッチ (正で足先が前方へ)
    knee      (θ3): 膝関節 (θ2 と同じ回転方向)

全関節 0 rad で脚はまっすぐ真下に伸び、足先は股関節から
(0, side * hip_offset, -(upper + lower)) の位置にある。
"""

from __future__ import annotations

import math

import numpy as np

from .config import LegGeometry


class IKError(ValueError):
    """目標位置に足先が届かない場合の例外。"""


def _rot_x(p: np.ndarray, angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([p[0], c * p[1] - s * p[2], s * p[1] + c * p[2]])


def leg_joint_positions(
    angles: np.ndarray, geom: LegGeometry, side: int
) -> np.ndarray:
    """股関節基準の [股関節, 大腿付け根, 膝, 足先] の位置 (4x3) を返す。"""
    a1, a2, a3 = angles
    d = side * geom.hip_offset
    root = np.array([0.0, d, 0.0])
    knee = root + np.array(
        [geom.upper * math.sin(a2), 0.0, -geom.upper * math.cos(a2)]
    )
    foot = knee + np.array(
        [geom.lower * math.sin(a2 + a3), 0.0, -geom.lower * math.cos(a2 + a3)]
    )
    pts = [np.zeros(3), root, knee, foot]
    return np.array([_rot_x(p, a1) for p in pts])


def forward_kinematics(angles: np.ndarray, geom: LegGeometry, side: int) -> np.ndarray:
    """関節角から股関節基準の足先位置を求める。"""
    return leg_joint_positions(angles, geom, side)[-1]


def inverse_kinematics(
    foot: np.ndarray, geom: LegGeometry, side: int, knee_direction: int = 1
) -> np.ndarray:
    """股関節基準の足先位置から関節角 [θ1, θ2, θ3] を求める。

    届かない位置が指定された場合は IKError を送出する。
    """
    x, y, z = (float(v) for v in foot)
    d = side * geom.hip_offset

    # 前額面 (y-z): 外転角を求める
    r2 = y * y + z * z - d * d
    if r2 < 0:
        raise IKError(f"足先が股関節に近すぎます: {foot}")
    zp = -math.sqrt(r2)  # 外転後の座標系での足先高さ (下向き)
    a1 = math.atan2(z, y) - math.atan2(zp, d)
    a1 = math.atan2(math.sin(a1), math.cos(a1))

    # 矢状面 (x-z'): 2リンクの逆運動学
    l1, l2 = geom.upper, geom.lower
    dist2 = x * x + zp * zp
    cos_knee = (dist2 - l1 * l1 - l2 * l2) / (2 * l1 * l2)
    if not -1.0 - 1e-9 <= cos_knee <= 1.0 + 1e-9:
        raise IKError(f"足先が可動範囲外です: {foot}")
    cos_knee = min(1.0, max(-1.0, cos_knee))
    a3 = knee_direction * math.acos(cos_knee)
    a2 = math.atan2(x, -zp) - math.atan2(l2 * math.sin(a3), l1 + l2 * math.cos(a3))
    return np.array([a1, a2, a3])
