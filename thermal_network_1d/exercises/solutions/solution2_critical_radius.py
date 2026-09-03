"""Q2 解答: 保温材の臨界半径(円筒座標の落とし穴)"""

import _path  # noqa: F401
import numpy as np

from thermalnet import r_cond_cylinder, r_conv

K_INS, H, LENGTH = 0.15, 12.0, 1.0
R_CABLE = 0.003          # ケーブル外半径 [m]

print("Q2: 保温材の臨界半径\n")


def total_resistance(r_out: float) -> tuple:
    """被覆外半径 r_out のときの (伝導, 対流, 合計) 熱抵抗 [K/W]。"""
    r_cond = (0.0 if r_out <= R_CABLE
              else r_cond_cylinder(R_CABLE, r_out, K_INS, LENGTH))
    r_out_conv = r_conv(H, 2 * np.pi * r_out * LENGTH)
    return r_cond, r_out_conv, r_cond + r_out_conv


radii = np.linspace(R_CABLE, 0.030, 400)
totals = np.array([total_resistance(r)[2] for r in radii])
i_min = int(np.argmin(totals))
r_cr_numeric = radii[i_min]
r_cr_theory = K_INS / H

print("   被覆外半径 [mm]   R_cond    R_conv    R_total [K/W]")
for r in (0.003, 0.006, 0.0125, 0.020, 0.030):
    rc, rv, rt = total_resistance(r)
    print(f"    {1000 * r:10.1f}   {rc:8.4f}  {rv:8.4f}  {rt:10.4f}")

print(f"""
  数値的に求めた臨界半径 = {1000 * r_cr_numeric:.2f} mm
  理論値 r_cr = k/h      = {1000 * r_cr_theory:.2f} mm

  裸のケーブル(半径 {1000 * R_CABLE:.0f} mm)の熱抵抗 {total_resistance(R_CABLE)[2]:.4f} K/W
  臨界半径まで巻いたとき             {total_resistance(r_cr_theory)[2]:.4f} K/W
  -> 保温材を巻いたのに熱抵抗が
     {100 * (1 - total_resistance(r_cr_theory)[2] / total_resistance(R_CABLE)[2]):.0f} % 下がった = 放熱が増えた

  理由: 被覆を厚くすると
    * 伝導抵抗 ln(r/r1)/(2 pi k L) は増える(対数なのでゆっくり)
    * 対流抵抗 1/(h 2 pi r L) は減る(表面積が増えるので 1/r で急に)
  細い管では後者が勝つ。r > k/h になって初めて保温効果が出る。
  電線の被覆はむしろ放熱に有利、という有名な例。

  1D 熱回路で考えると「面積が場所によって変わる」という円筒/球の特徴が
  素直に効いてくることが分かる。平板の直感をそのまま持ち込まないこと。
""")
