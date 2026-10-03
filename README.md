# Automation
Develop code to automate the process using Python.

## quadruped — 四足歩行ロボット

3自由度 × 4脚の四足歩行ロボット用の運動学・歩容生成ライブラリです。

| モジュール | 内容 |
| --- | --- |
| `quadruped/config.py` | 寸法・脚配置 (`RobotConfig`, `LegGeometry`) |
| `quadruped/kinematics.py` | 脚の順運動学・逆運動学 |
| `quadruped/gait.py` | 歩容生成 (stand / trot / walk / pace / bound) |
| `quadruped/body.py` | ボディ姿勢 (ロール・ピッチ・ヨー・並進) |
| `quadruped/robot.py` | 上記をまとめたコントローラ `Quadruped` |
| `quadruped/servo.py` | 関節角 → サーボパルス幅 [µs] の変換 |

座標系は x: 前方, y: 左方, z: 上方 (m, rad)。脚名は `FL`, `FR`, `RL`, `RR`。

### セットアップ

```bash
pip install -r requirements.txt
```

### 使い方

```python
from quadruped import Quadruped, GaitGenerator, VelocityCommand

robot = Quadruped(gait=GaitGenerator("trot", period=0.5, step_height=0.04))
cmd = VelocityCommand(vx=0.2, vy=0.0, yaw_rate=0.3)
for _ in range(100):
    state = robot.update(0.01, cmd)
    print(state.joint_angles["FL"])  # [外転, 股, 膝] [rad]
```

### コマンドライン

```bash
python -m quadruped --gait walk --vx 0.1            # 関節角の表を出力
python -m quadruped --gait trot --animate           # 3D アニメーション表示
python -m quadruped --gait trot --save trot.gif     # GIF に保存
python -m quadruped --help
```

### テスト

```bash
python -m pytest
```
