"""空気物性と対流熱伝達率の相関式。

CFD を 1D 化するというのは、要するに
    「CFD が数値的に解いていた境界層を、相関式 h = Nu k / L で置き換える」
ということ。ここでは代表的な相関式を最小限そろえる。
出典は Incropera & DeWitt, "Fundamentals of Heat and Mass Transfer"。
"""

from __future__ import annotations

import math

import numpy as np

from .elements import T_ABS

G_ACC = 9.80665  # 重力加速度 [m/s^2]

# 1 atm 乾燥空気の物性(Incropera 表 A.4)。温度は絶対温度 [K]。
_T_TABLE = np.array([250.0, 300.0, 350.0, 400.0, 450.0, 500.0])
_RHO = np.array([1.3947, 1.1614, 0.9950, 0.8711, 0.7740, 0.6964])   # [kg/m^3]
_CP = np.array([1006.0, 1007.0, 1009.0, 1014.0, 1021.0, 1030.0])    # [J/(kg K)]
_K = np.array([0.02227, 0.02624, 0.03003, 0.03365, 0.03707, 0.04038])  # [W/(m K)]
_NU = np.array([11.44e-6, 15.89e-6, 20.92e-6, 26.41e-6, 32.39e-6, 38.79e-6])  # [m^2/s]
_PR = np.array([0.720, 0.707, 0.700, 0.690, 0.686, 0.684])


def air_properties(t_celsius: float) -> dict:
    """1 atm 空気の物性を線形補間で返す(表の範囲外は端点で頭打ち)。"""
    tk = float(np.clip(t_celsius + T_ABS, _T_TABLE[0], _T_TABLE[-1]))
    rho = float(np.interp(tk, _T_TABLE, _RHO))
    cp = float(np.interp(tk, _T_TABLE, _CP))
    k = float(np.interp(tk, _T_TABLE, _K))
    nu = float(np.interp(tk, _T_TABLE, _NU))
    pr = float(np.interp(tk, _T_TABLE, _PR))
    return {"T_K": tk, "rho": rho, "cp": cp, "k": k, "nu": nu,
            "alpha": nu / pr, "Pr": pr, "beta": 1.0 / tk, "mu": rho * nu}


def rayleigh(t_surf_c: float, t_inf_c: float, length: float) -> float:
    """レイリー数 Ra = g beta dT L^3 / (nu alpha)(膜温度で物性評価)。"""
    p = air_properties(0.5 * (t_surf_c + t_inf_c))
    dt = abs(t_surf_c - t_inf_c)
    return G_ACC * p["beta"] * dt * length ** 3 / (p["nu"] * p["alpha"])


def h_natural_vertical_plate(t_surf_c: float, t_inf_c: float, height: float) -> float:
    """垂直平板の自然対流(Churchill & Chu, 全 Ra 域)。

        Nu = {0.825 + 0.387 Ra^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/27)}^2
    """
    if abs(t_surf_c - t_inf_c) < 1e-9:
        return 1e-6  # ゼロ割回避(実質断熱)
    p = air_properties(0.5 * (t_surf_c + t_inf_c))
    ra = rayleigh(t_surf_c, t_inf_c, height)
    denom = (1.0 + (0.492 / p["Pr"]) ** (9.0 / 16.0)) ** (8.0 / 27.0)
    nu_num = (0.825 + 0.387 * ra ** (1.0 / 6.0) / denom) ** 2
    return nu_num * p["k"] / height


def h_natural_horizontal_plate_up(t_surf_c: float, t_inf_c: float,
                                  char_length: float) -> float:
    """上向き高温水平面の自然対流。特性長 Lc = A / P。

        Nu = 0.54 Ra^(1/4)   (1e4 <= Ra <= 1e7)
        Nu = 0.15 Ra^(1/3)   (1e7 <  Ra <= 1e11)
    """
    if abs(t_surf_c - t_inf_c) < 1e-9:
        return 1e-6
    p = air_properties(0.5 * (t_surf_c + t_inf_c))
    ra = rayleigh(t_surf_c, t_inf_c, char_length)
    nu = 0.54 * ra ** 0.25 if ra <= 1e7 else 0.15 * ra ** (1.0 / 3.0)
    return nu * p["k"] / char_length


def h_forced_flat_plate(velocity: float, length: float, t_film_c: float) -> float:
    """平板上の強制対流(平均 Nu)。

        層流 (Re < 5e5) : Nu = 0.664 Re^(1/2) Pr^(1/3)
        混合           : Nu = (0.037 Re^(4/5) - 871) Pr^(1/3)
    """
    p = air_properties(t_film_c)
    re = velocity * length / p["nu"]
    if re < 5e5:
        nu = 0.664 * math.sqrt(re) * p["Pr"] ** (1.0 / 3.0)
    else:
        nu = (0.037 * re ** 0.8 - 871.0) * p["Pr"] ** (1.0 / 3.0)
    return nu * p["k"] / length


def h_duct(velocity: float, hydraulic_diameter: float, t_film_c: float,
           uniform_flux: bool = True) -> float:
    """ダクト内(ヒートシンクのフィン間流路)の強制対流。

        層流 (Re < 2300) : 十分発達 Nu = 4.36 (等熱流束) / 3.66 (等温)
        乱流             : Dittus-Boelter  Nu = 0.023 Re^(4/5) Pr^0.4
    """
    p = air_properties(t_film_c)
    re = velocity * hydraulic_diameter / p["nu"]
    if re < 2300.0:
        nu = 4.36 if uniform_flux else 3.66
    else:
        nu = 0.023 * re ** 0.8 * p["Pr"] ** 0.4
    return nu * p["k"] / hydraulic_diameter


def reynolds(velocity: float, length: float, t_film_c: float) -> float:
    return velocity * length / air_properties(t_film_c)["nu"]
