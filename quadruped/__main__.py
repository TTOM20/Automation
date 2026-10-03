"""コマンドラインからの実行: python -m quadruped --help"""

from __future__ import annotations

import argparse
import math

import numpy as np

from .body import BodyPose
from .config import LEG_NAMES
from .gait import GAITS, GaitGenerator, VelocityCommand
from .robot import Quadruped


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="四足歩行ロボットのシミュレーション")
    p.add_argument("--gait", choices=sorted(GAITS), default="trot")
    p.add_argument("--period", type=float, default=0.5, help="歩行周期 [s]")
    p.add_argument("--step-height", type=float, default=0.04, help="足上げ高さ [m]")
    p.add_argument("--vx", type=float, default=0.15, help="前進速度 [m/s]")
    p.add_argument("--vy", type=float, default=0.0, help="横移動速度 [m/s]")
    p.add_argument("--yaw-rate", type=float, default=0.0, help="旋回角速度 [rad/s]")
    p.add_argument("--roll", type=float, default=0.0, help="ボディのロール角 [deg]")
    p.add_argument("--pitch", type=float, default=0.0, help="ボディのピッチ角 [deg]")
    p.add_argument("--duration", type=float, default=2.0, help="シミュレーション時間 [s]")
    p.add_argument("--dt", type=float, default=0.02, help="制御周期 [s]")
    p.add_argument("--animate", action="store_true", help="matplotlib でアニメーション表示")
    p.add_argument("--save", metavar="FILE", help="アニメーションを GIF に保存")
    return p.parse_args(argv)


def build_robot(args: argparse.Namespace) -> tuple[Quadruped, VelocityCommand]:
    robot = Quadruped(gait=GaitGenerator(args.gait, args.period, args.step_height))
    robot.set_pose(BodyPose(roll=math.radians(args.roll), pitch=math.radians(args.pitch)))
    return robot, VelocityCommand(args.vx, args.vy, args.yaw_rate)


def print_table(robot: Quadruped, cmd: VelocityCommand, duration: float, dt: float) -> None:
    header = "  t[s] " + " ".join(f"{leg:>22}" for leg in LEG_NAMES)
    print(header)
    print("       " + " ".join(f"{'abd    hip   knee [deg]':>22}" for _ in LEG_NAMES))
    for _ in range(int(round(duration / dt))):
        state = robot.update(dt, cmd)
        cols = []
        for leg in LEG_NAMES:
            a = np.degrees(state.joint_angles[leg])
            mark = "*" if state.stance[leg] else " "
            cols.append(f"{mark}{a[0]:7.1f}{a[1]:7.1f}{a[2]:7.1f}".rjust(22))
        print(f"{state.time:6.2f} " + " ".join(cols))
    print("(* = 接地中)")


def animate(robot: Quadruped, cmd: VelocityCommand, duration: float, dt: float, save: str | None) -> None:
    import matplotlib

    if save:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(projection="3d")
    lim = 0.25
    body_line, = ax.plot([], [], [], "k-", lw=3)
    leg_lines = {leg: ax.plot([], [], [], "o-", lw=2, label=leg)[0] for leg in LEG_NAMES}
    title = ax.set_title("")
    ax.legend(loc="upper right")

    def setup_axes() -> None:
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(0, 2 * lim)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")

    setup_axes()
    ground = robot.config.stand_height

    def draw(_frame: int):
        robot.update(dt, cmd)
        corners = robot.body_corners()
        corners = np.vstack([corners, corners[:1]])
        body_line.set_data_3d(corners[:, 0], corners[:, 1], corners[:, 2] + ground)
        for leg, pts in robot.joint_positions().items():
            leg_lines[leg].set_data_3d(pts[:, 0], pts[:, 1], pts[:, 2] + ground)
            leg_lines[leg].set_linestyle("-" if robot.state.stance[leg] else "--")
        title.set_text(f"{robot.gait.pattern.name}  t={robot.state.time:.2f}s")
        return [body_line, title, *leg_lines.values()]

    frames = int(round(duration / dt))
    anim = FuncAnimation(fig, draw, frames=frames, interval=dt * 1000, blit=False)
    if save:
        anim.save(save, writer=PillowWriter(fps=int(round(1 / dt))))
        print(f"保存しました: {save}")
    else:
        plt.show()


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    robot, cmd = build_robot(args)
    if args.animate or args.save:
        animate(robot, cmd, args.duration, args.dt, args.save)
    else:
        print_table(robot, cmd, args.duration, args.dt)


if __name__ == "__main__":
    main()
