"""ミニ CFD(伝導のみ)--- 1D 回路の答え合わせに使う軸対称有限体積ソルバ。

「CFD 解析を 1D 化する」練習をするには、比較対象となる多次元の解が要る。
ここでは商用 CFD の代わりに、円板スプレッダの軸対称定常熱伝導を
有限体積法で解く 100 行のソルバを用意する。

モデル(lesson07 / lesson08 で使用):

      z=0   +---------------------------+   上面: r<a に一様熱流束 Q/(pi a^2)
            |####|                      |        r>a は断熱
            |    |     スプレッダ k     |
      z=t   +---------------------------+   下面: 熱伝達率 h で T_inf へ
            r=0                       r=b     側面 r=b: 断熱(対称セル境界)

有限体積法の離散式そのものが熱回路網である:
    sum_faces G_face (T_nb - T_P) + S_P = 0,   G_face = k A_face / d
これは thermalnet.ThermalNetwork の節点方程式と完全に同じ形をしている。
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


def axisymmetric_spreader(a: float, b: float, t: float, k: float, h: float,
                          q_total: float = 1.0, nr: int = 80, nz: int = 40,
                          t_inf: float = 0.0) -> dict:
    """円板スプレッダの温度場を解き、拡がり抵抗を抽出する。

    Parameters
    ----------
    a, b, t : 熱源半径 / 板半径 / 板厚 [m]
    k       : 板の熱伝導率 [W/(m K)]
    h       : 下面の等価熱伝達率 [W/(m^2 K)]
    q_total : 全発熱量 [W]
    nr, nz  : 半径方向 / 板厚方向のセル数
    t_inf   : 下面側の周囲温度 [degC]

    Returns
    -------
    dict : T(2D 配列 [nr, nz]), r, z, R_total, R_1d, R_spread ほか

    抵抗の定義(業界標準の分解):
        R_total  = (熱源平均温度 - T_inf) / Q
        R_1d     = t / (k pi b^2) + 1 / (h pi b^2)     <- 素朴な 1D 回路
        R_spread = R_total - R_1d                       <- 1D では抜け落ちる分
    """
    dr = b / nr
    dz = t / nz
    r_face = np.linspace(0.0, b, nr + 1)                 # セル境界 [m]
    r_cell = 0.5 * (r_face[:-1] + r_face[1:])            # セル中心 [m]
    z_cell = (np.arange(nz) + 0.5) * dz
    # 各リングセルの水平断面積 pi (r_out^2 - r_in^2)
    a_cell = np.pi * (r_face[1:] ** 2 - r_face[:-1] ** 2)

    def node(i, j):
        return i * nz + j

    n = nr * nz
    rows, cols, vals = [], [], []
    rhs = np.zeros(n)

    def add(i_node, j_node, value):
        rows.append(i_node)
        cols.append(j_node)
        vals.append(value)

    # 熱源: セルと半径 a の円の重なり面積に比例配分(格子依存を減らす)
    r_clip = np.clip(r_face, 0.0, a)
    a_src = np.pi * (r_clip[1:] ** 2 - r_clip[:-1] ** 2)
    flux_area_total = np.pi * a ** 2
    src = q_total * a_src / flux_area_total

    for i in range(nr):
        for j in range(nz):
            p = node(i, j)
            diag = 0.0
            # 半径方向(内側・外側)。r=0 と r=b は断熱なので寄与なし
            if i > 0:
                g = k * (2.0 * np.pi * r_face[i] * dz) / dr
                diag += g
                add(p, node(i - 1, j), -g)
            if i < nr - 1:
                g = k * (2.0 * np.pi * r_face[i + 1] * dz) / dr
                diag += g
                add(p, node(i + 1, j), -g)
            # 板厚方向
            if j > 0:
                g = k * a_cell[i] / dz
                diag += g
                add(p, node(i, j - 1), -g)
            else:
                rhs[p] += src[i]           # 上面: 熱流束流入
            if j < nz - 1:
                g = k * a_cell[i] / dz
                diag += g
                add(p, node(i, j + 1), -g)
            else:
                # 下面: 対流境界。半セル分の伝導と対流の直列合成
                g = 1.0 / (dz / (2.0 * k * a_cell[i]) + 1.0 / (h * a_cell[i]))
                diag += g
                rhs[p] += g * t_inf
            add(p, p, diag)

    mat = csr_matrix((vals, (rows, cols)), shape=(n, n))
    temp = spsolve(mat.tocsc(), rhs).reshape(nr, nz)

    # 熱源直下(上面 j=0, r<a)の面積加重平均温度・最高温度
    mask = a_src > 0
    t_src_mean = float(np.sum(temp[mask, 0] * a_src[mask]) / np.sum(a_src[mask]))
    t_src_max = float(temp[0, 0])
    area_b = np.pi * b ** 2
    r_total = (t_src_mean - t_inf) / q_total
    r_1d = t / (k * area_b) + 1.0 / (h * area_b)
    return {
        "T": temp, "r": r_cell, "z": z_cell, "a_cell": a_cell,
        "T_source_mean": t_src_mean, "T_source_max": t_src_max,
        "T_base_mean": float(np.sum(temp[:, -1] * a_cell) / area_b),
        "R_total": r_total, "R_1d": r_1d, "R_spread": r_total - r_1d,
        "R_total_max": (t_src_max - t_inf) / q_total,
        "grid": (nr, nz),
    }
