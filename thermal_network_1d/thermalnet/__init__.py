"""thermalnet --- 1D 熱回路網の学習用ミニライブラリ。"""

from .network import Branch, Node, Solution, ThermalNetwork, TransientSolution
from . import elements
from .elements import (
    SIGMA,
    T_ABS,
    biot_number,
    capacitance,
    fin_efficiency,
    h_radiation,
    r_cond_cylinder,
    r_cond_plane,
    r_cond_sphere,
    r_contact,
    r_conv,
    r_fin_array,
    r_flow,
    r_interface,
    r_radiation,
    r_spreading_disk,
    r_spreading_halfspace,
)

__all__ = [
    "ThermalNetwork", "Node", "Branch", "Solution", "TransientSolution",
    "elements", "SIGMA", "T_ABS",
    "r_cond_plane", "r_cond_cylinder", "r_cond_sphere", "r_contact", "r_interface",
    "r_conv", "r_flow", "r_radiation", "h_radiation",
    "fin_efficiency", "r_fin_array",
    "r_spreading_halfspace", "r_spreading_disk",
    "capacitance", "biot_number",
]
