"""熱抵抗エレメント集 --- CFD を 1D 熱回路に落とすときの「部品カタログ」。

すべての関数は SI 単位、戻り値の熱抵抗は [K/W] を返す。
温度は特記なき限り [degC] で受け取り、内部で必要に応じ [K] に直す。

熱回路のアナロジー:
    電気              熱
    ------------------------------------
    電圧 V [V]        温度 T [K]
    電流 I [A]        熱流 Q [W]
    抵抗 R [Ohm]      熱抵抗 R [K/W]
    静電容量 C [F]    熱容量 C [J/K]
    オームの法則      Q = (T1 - T2) / R
"""

from __future__ import annotations

import math

SIGMA = 5.670374419e-8  # ステファン・ボルツマン定数 [W/(m^2 K^4)]
T_ABS = 273.15


# ---------------------------------------------------------------- 熱伝導
def r_cond_plane(thickness: float, k: float, area: float) -> float:
    """平板の熱伝導抵抗 R = L / (k A)。

    thickness [m], k [W/(m K)], area [m^2]
    """
    _positive(thickness=thickness, k=k, area=area)
    return thickness / (k * area)


def r_cond_cylinder(r_in: float, r_out: float, k: float, length: float) -> float:
    """円筒殻(配管の管壁・断熱材)の熱伝導抵抗 R = ln(r2/r1) / (2 pi k L)。"""
    _positive(r_in=r_in, r_out=r_out, k=k, length=length)
    if r_out <= r_in:
        raise ValueError("r_out は r_in より大きい必要があります")
    return math.log(r_out / r_in) / (2.0 * math.pi * k * length)


def r_cond_sphere(r_in: float, r_out: float, k: float) -> float:
    """球殻の熱伝導抵抗 R = (1/r1 - 1/r2) / (4 pi k)。"""
    _positive(r_in=r_in, r_out=r_out, k=k)
    if r_out <= r_in:
        raise ValueError("r_out は r_in より大きい必要があります")
    return (1.0 / r_in - 1.0 / r_out) / (4.0 * math.pi * k)


def r_contact(contact_conductance: float, area: float) -> float:
    """接触熱抵抗 R = 1 / (h_c A)。

    contact_conductance [W/(m^2 K)] : 接触熱コンダクタンス
    TIM のカタログ値が熱抵抗率 [m^2 K/W] で与えられる場合は
    r_interface(resistivity, area) を使う。
    """
    _positive(contact_conductance=contact_conductance, area=area)
    return 1.0 / (contact_conductance * area)


def r_interface(resistivity: float, area: float) -> float:
    """界面熱抵抗率 [m^2 K/W] と面積から R [K/W] を作る(TIM のカタログ値向け)。"""
    _positive(resistivity=resistivity, area=area)
    return resistivity / area


# ---------------------------------------------------------------- 対流
def r_conv(h: float, area: float) -> float:
    """対流熱抵抗 R = 1 / (h A)。h [W/(m^2 K)], area [m^2]"""
    _positive(h=h, area=area)
    return 1.0 / (h * area)


def r_flow(mass_flow: float, cp: float) -> float:
    """流体の熱容量流量による抵抗(caloric resistance) R = 1 / (mdot cp)。

    強制空冷を 1D 化するとき必須。流路を通過する空気は Q を受け取って
    dT = Q / (mdot cp) だけ昇温する。この昇温を「入口温度と出口温度の間の
    抵抗」として回路に入れないと、放熱器の性能を過大評価してしまう。
    """
    _positive(mass_flow=mass_flow, cp=cp)
    return 1.0 / (mass_flow * cp)


# ---------------------------------------------------------------- 輻射
def h_radiation(emissivity: float, t_surf_c: float, t_surr_c: float,
                view_factor: float = 1.0) -> float:
    """輻射を線形化した等価熱伝達率 h_r [W/(m^2 K)]。

        q" = eps F sigma (Ts^4 - Tsur^4) = h_r (Ts - Tsur)
        h_r = eps F sigma (Ts + Tsur)(Ts^2 + Tsur^2)     (絶対温度)

    温度に依存するので、非線形反復の中で毎回更新して使う。
    """
    ts = t_surf_c + T_ABS
    tr = t_surr_c + T_ABS
    return emissivity * view_factor * SIGMA * (ts + tr) * (ts * ts + tr * tr)


def r_radiation(emissivity: float, area: float, t_surf_c: float, t_surr_c: float,
                view_factor: float = 1.0) -> float:
    """輻射熱抵抗 R = 1 / (h_r A)。"""
    return r_conv(h_radiation(emissivity, t_surf_c, t_surr_c, view_factor), area)


# ---------------------------------------------------------------- フィン
def fin_efficiency(h: float, k: float, thickness: float, length: float) -> float:
    """矩形ストレートフィンの効率 eta = tanh(m Lc) / (m Lc)。

    先端断熱の補正長さ Lc = L + t/2 を使う。
    """
    _positive(h=h, k=k, thickness=thickness, length=length)
    lc = length + thickness / 2.0
    m = math.sqrt(2.0 * h / (k * thickness))  # 単位幅あたり P/Ac = 2/t
    return math.tanh(m * lc) / (m * lc)


def r_fin_array(h: float, k: float, fin_thickness: float, fin_length: float,
                fin_width: float, n_fins: float, base_area: float) -> float:
    """フィンアレイ全体の対流抵抗。

        R = 1 / [ h (A_base_exposed + N eta A_fin) ]

    A_fin = 2 (L + t/2) W  (両面), A_base_exposed = base_area - N t W
    """
    eta = fin_efficiency(h, k, fin_thickness, fin_length)
    a_fin = 2.0 * (fin_length + fin_thickness / 2.0) * fin_width
    a_base = base_area - n_fins * fin_thickness * fin_width
    if a_base < 0:
        raise ValueError("フィン根元面積がベース面積を超えています")
    return 1.0 / (h * (a_base + n_fins * eta * a_fin))


# ---------------------------------------------------------------- 拡がり抵抗
def r_spreading_halfspace(radius: float, k: float, isothermal: bool = True) -> float:
    """半無限体に接する円形熱源の拡がり抵抗。

    等温源      R = 1 / (4 k a)
    等熱流束源  R = 8 / (3 pi^2 k a)
    3D CFD を 1D 化するとき、「面積で割った 1D 伝導抵抗」に加えて
    必ずこの項が要る(小さい熱源ほど支配的)。
    """
    _positive(radius=radius, k=k)
    if isothermal:
        return 1.0 / (4.0 * k * radius)
    return 8.0 / (3.0 * math.pi ** 2 * k * radius)


def r_spreading_disk(a: float, b: float, t: float, k: float, h: float) -> float:
    """有限厚みの円板スプレッダの拡がり抵抗(Lee et al., 1995 の近似式)。

    a : 熱源半径 [m], b : 円板半径 [m], t : 板厚 [m]
    k : 板の熱伝導率 [W/(m K)], h : 裏面の等価熱伝達率 [W/(m^2 K)]

    Reference: S. Lee, S. Song, V. Au, K. Moran,
      "Constriction/Spreading Resistance Model for Electronics Packaging",
      ASME/JSME Thermal Engineering Conference, 1995.

    注意: これは相関式であり万能ではない。lesson07 で軸対称の数値解と
    突き合わせて誤差を自分の目で確認すること。
    """
    _positive(a=a, b=b, t=t, k=k, h=h)
    if a >= b:
        raise ValueError("熱源半径 a は板半径 b より小さい必要があります")
    eps = a / b
    tau = t / b
    bi = h * b / k
    lam = math.pi + 1.0 / (eps * math.sqrt(math.pi))
    phi = ((math.tanh(lam * tau) + lam / bi)
           / (1.0 + (lam / bi) * math.tanh(lam * tau)))
    psi = 0.5 * (1.0 - eps) ** 1.5 * phi
    return psi / (k * a * math.sqrt(math.pi))


# ---------------------------------------------------------------- 熱容量
def capacitance(rho: float, cp: float, volume: float) -> float:
    """熱容量 C = rho cp V [J/K]。過渡解析の「コンデンサ」。"""
    _positive(rho=rho, cp=cp, volume=volume)
    return rho * cp * volume


def biot_number(h: float, k: float, length_scale: float) -> float:
    """ビオ数 Bi = h Lc / k。

    Bi < 0.1 なら物体内部の温度分布を無視して 1 ノード(集中定数)で
    扱ってよい、というのが 1D 化の最初の判定基準。
    """
    _positive(h=h, k=k, length_scale=length_scale)
    return h * length_scale / k


def _positive(**kwargs) -> None:
    for name, value in kwargs.items():
        if value <= 0:
            raise ValueError(f"{name} は正の値である必要があります (given {value})")
