import math

import numpy as np
import pytest

from quadruped import (
    GAITS,
    LEG_NAMES,
    BodyPose,
    GaitGenerator,
    IKError,
    LegGeometry,
    Quadruped,
    RobotConfig,
    VelocityCommand,
    forward_kinematics,
    inverse_kinematics,
)
from quadruped.__main__ import main
from quadruped.servo import PrintDriver, apply_joint_angles, default_calibration

GEOM = LegGeometry()


@pytest.mark.parametrize("side", [1, -1])
@pytest.mark.parametrize("knee", [1, -1])
@pytest.mark.parametrize(
    "target",
    [(0.0, 0.0, -0.17), (0.05, 0.02, -0.15), (-0.04, -0.03, -0.2), (0.03, 0.0, -0.1)],
)
def test_ik_fk_roundtrip(side, knee, target):
    foot = np.array(target) + np.array([0, side * GEOM.hip_offset, 0])
    angles = inverse_kinematics(foot, GEOM, side, knee)
    np.testing.assert_allclose(forward_kinematics(angles, GEOM, side), foot, atol=1e-9)
    assert np.sign(angles[2]) in (0, knee)


def test_straight_leg_is_zero_angles():
    foot = np.array([0.0, GEOM.hip_offset, -(GEOM.upper + GEOM.lower)])
    np.testing.assert_allclose(inverse_kinematics(foot, GEOM, 1), [0, 0, 0], atol=1e-6)


def test_unreachable_raises():
    with pytest.raises(IKError):
        inverse_kinematics(np.array([0.0, 0.05, -0.5]), GEOM, 1)
    with pytest.raises(IKError):
        inverse_kinematics(np.array([0.0, 0.0, 0.0]), GEOM, 1)


def test_body_pose_roundtrip():
    pose = BodyPose(0.01, -0.02, 0.03, 0.1, -0.2, 0.3)
    p = np.array([0.1, -0.05, -0.15])
    np.testing.assert_allclose(pose.to_reference_frame(pose.to_body_frame(p)), p)


@pytest.mark.parametrize("name", sorted(GAITS))
def test_gaits_produce_valid_angles(name):
    robot = Quadruped(gait=GaitGenerator(name))
    cmd = VelocityCommand(0.2, 0.05, 0.3)
    for _ in range(100):
        state = robot.update(0.01, cmd)
        for leg in LEG_NAMES:
            assert np.all(np.isfinite(state.joint_angles[leg]))


def test_stance_ratio_matches_duty():
    for name, pattern in GAITS.items():
        gen = GaitGenerator(name, period=1.0)
        n, count = 400, 0
        for _ in range(n):
            gen.advance(1.0 / n)
            count += sum(gen.in_stance(leg) for leg in LEG_NAMES)
        assert count / (n * 4) == pytest.approx(pattern.duty, abs=0.01)


def test_trot_diagonal_pairs_in_sync():
    gen = GaitGenerator("trot")
    for _ in range(50):
        gen.advance(0.013)
        assert gen.in_stance("FL") == gen.in_stance("RR")
        assert gen.in_stance("FR") == gen.in_stance("RL")
        assert gen.in_stance("FL") != gen.in_stance("FR")


def test_walk_keeps_three_feet_down():
    gen = GaitGenerator("walk")
    for _ in range(200):
        gen.advance(0.005)
        assert sum(gen.in_stance(leg) for leg in LEG_NAMES) >= 3


def test_stance_foot_moves_backward_and_stays_on_ground():
    cfg = RobotConfig()
    gen = GaitGenerator("trot", period=0.5)
    cmd = VelocityCommand(vx=0.2)
    neutral = cfg.neutral_foot_position("FL")
    xs = []
    while gen.in_stance("FL"):
        p = gen.foot_position("FL", neutral, cmd)
        assert p[2] == pytest.approx(neutral[2])
        xs.append(p[0])
        gen.advance(0.01)
    assert all(b < a for a, b in zip(xs, xs[1:]))
    # 接地中の足先速度はボディ速度の逆向き
    assert (xs[0] - xs[-1]) / (0.01 * (len(xs) - 1)) == pytest.approx(0.2, rel=1e-6)


def test_swing_lifts_foot_and_trajectory_is_continuous():
    cfg = RobotConfig()
    gen = GaitGenerator("trot", period=0.5, step_height=0.04)
    cmd = VelocityCommand(vx=0.2)
    neutral = cfg.neutral_foot_position("FL")
    prev = gen.foot_position("FL", neutral, cmd)
    max_z = prev[2]
    for _ in range(500):
        gen.advance(0.001)
        p = gen.foot_position("FL", neutral, cmd)
        assert np.linalg.norm(p - prev) < 0.002
        max_z = max(max_z, p[2])
        prev = p
    assert max_z - neutral[2] == pytest.approx(0.04, abs=1e-4)


def test_body_pitch_changes_leg_angles():
    robot = Quadruped()
    base = robot.state.joint_angles["FL"].copy()
    robot.set_pose(BodyPose(pitch=math.radians(10)))
    robot.update(0.0)
    assert not np.allclose(base, robot.state.joint_angles["FL"])


def test_servo_mapping():
    robot = Quadruped()
    driver = PrintDriver()
    apply_joint_angles(driver, robot.state.joint_angles, default_calibration())
    assert sorted(driver.last) == list(range(12))
    assert all(500 <= us <= 2500 for us in driver.last.values())


def test_cli_table(capsys):
    main(["--gait", "walk", "--duration", "0.1"])
    out = capsys.readouterr().out
    assert "FL" in out and "接地中" in out


def test_cli_save_gif(tmp_path):
    out = tmp_path / "trot.gif"
    main(["--duration", "0.1", "--save", str(out)])
    assert out.stat().st_size > 0
