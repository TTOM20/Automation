"""全レッスン・全解答・テストを順に実行する(教材の動作確認用)。

    python run_all.py            # 全部
    python run_all.py lessons    # レッスンだけ
    python run_all.py tests      # テストだけ
    python run_all.py models     # モデルファイルだけ
    python run_all.py exercises  # 演習の解答だけ
"""

import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))


def run(script, cwd):
    rel = os.path.relpath(os.path.join(cwd, script), ROOT)
    print(f"\n{'=' * 78}\n>>> {rel}\n{'=' * 78}")
    result = subprocess.run([sys.executable, script], cwd=cwd)
    return result.returncode


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    failures = []

    if what in ("all", "tests"):
        test_dir = os.path.join(ROOT, "tests")
        for name in sorted(os.listdir(test_dir)):
            if name.startswith("test_") and name.endswith(".py"):
                if run(name, test_dir):
                    failures.append(name)

    if what in ("all", "lessons"):
        lesson_dir = os.path.join(ROOT, "lessons")
        for name in sorted(os.listdir(lesson_dir)):
            if name.startswith("lesson") and name.endswith(".py"):
                if run(name, lesson_dir):
                    failures.append(name)

    if what in ("all", "exercises"):
        sol_dir = os.path.join(ROOT, "exercises", "solutions")
        for name in sorted(os.listdir(sol_dir)):
            if name.startswith("solution") and name.endswith(".py"):
                if run(name, sol_dir):
                    failures.append(name)

    if what in ("all", "models"):
        model_dir = os.path.join(ROOT, "models")
        for name in sorted(os.listdir(model_dir)):
            if name.endswith((".yaml", ".yml", ".json")):
                rel = os.path.join("models", name)
                print(f"\n{'=' * 78}\n>>> python -m thermalnet.model {rel}\n"
                      f"{'=' * 78}")
                if subprocess.run([sys.executable, "-m", "thermalnet.model", rel],
                                  cwd=ROOT).returncode:
                    failures.append(rel)

    print(f"\n{'=' * 78}")
    if failures:
        print(f"失敗: {', '.join(failures)}")
        return 1
    print("すべて正常に実行されました。figures/ に図が出力されています。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
