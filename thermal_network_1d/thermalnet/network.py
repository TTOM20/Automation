"""学習用の熱回路網ソルバ(節点法)。

CFD の伝熱解析コードも、突き詰めれば
    「各コントロールボリューム = ノード」
    「面を通る伝熱 = ノード間のコンダクタンス」
    「発熱 = ノードへの電流源」
    「rho cp V = ノードの熱容量」
という同じ構造をしている。ここでは 100 行程度でその骨格を実装し、
1D 熱回路網と CFD 離散化が地続きであることを体感する。

支配方程式(節点 i について):
    C_i dT_i/dt = Q_i + sum_j (T_j - T_i) / R_ij
定常なら左辺 = 0。行列で書くと  A T = b,  A = ラプラシアン型コンダクタンス行列。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Sequence

import numpy as np

Number = float
RSpec = float | Callable[[float, float], float]
QSpec = float | Callable[[float], float]


@dataclass
class Node:
    name: str
    capacitance: float = 0.0      # [J/K]
    fixed: float | None = None    # 温度固定境界 [degC]
    heat: QSpec = 0.0             # 発熱 [W] (定数 or t の関数)


@dataclass
class Branch:
    a: str
    b: str
    resistance: RSpec             # [K/W] (定数 or (Ta, Tb) の関数)
    label: str = ""

    def value(self, ta: float, tb: float) -> float:
        r = self.resistance(ta, tb) if callable(self.resistance) else self.resistance
        if r <= 0:
            raise ValueError(f"熱抵抗は正である必要があります: {self.a}-{self.b} R={r}")
        return r


@dataclass
class Solution:
    """定常解。solution['node'] で温度 [degC] を取り出せる。"""
    temperatures: Dict[str, float]
    network: "ThermalNetwork"
    iterations: int = 1
    residual: float = 0.0

    def __getitem__(self, name: str) -> float:
        return self.temperatures[name]

    def flow(self, a: str, b: str, label: str | None = None) -> float:
        """枝 a->b を流れる熱流 [W](正なら a から b へ)。

        同じ 2 ノード間に並列枝が複数あるときは label で指定する。
        """
        return self.network.branch_flow(a, b, self.temperatures, label)

    def flow_of(self, label: str) -> float:
        """ラベルで枝を指定して熱流 [W] を得る。"""
        for br in self.network.branches:
            if br.label == label:
                return self.network.branch_flow(br.a, br.b, self.temperatures, label)
        raise KeyError(f"ラベル {label} の枝がありません")

    def report(self, title: str = "定常解") -> str:
        lines = [title, "-" * len(title) * 2, "  ノード温度 [degC]"]
        for name, t in self.temperatures.items():
            tag = " (固定)" if self.network.nodes[name].fixed is not None else ""
            lines.append(f"    {name:<22s} {t:9.3f}{tag}")
        lines.append("  枝の熱流 [W] と熱抵抗 [K/W]")
        for br in self.network.branches:
            ta, tb = self.temperatures[br.a], self.temperatures[br.b]
            r = br.value(ta, tb)
            q = (ta - tb) / r
            name = br.label or f"{br.a}->{br.b}"
            lines.append(f"    {name:<22s} Q={q:8.3f}  R={r:8.4f}  dT={ta - tb:7.3f}")
        return "\n".join(lines)


@dataclass
class TransientSolution:
    """過渡解。t [s] と各ノードの温度履歴を保持する。"""
    time: np.ndarray
    temperatures: Dict[str, np.ndarray]
    names: List[str] = field(default_factory=list)

    def __getitem__(self, name: str) -> np.ndarray:
        return self.temperatures[name]

    def final(self) -> Dict[str, float]:
        return {k: float(v[-1]) for k, v in self.temperatures.items()}


class ThermalNetwork:
    """熱回路網。ノードと抵抗を登録して解くだけの素朴な実装。"""

    def __init__(self, name: str = "network") -> None:
        self.name = name
        self.nodes: Dict[str, Node] = {}
        self.branches: List[Branch] = []

    # ------------------------------------------------------------ 構築
    def add_node(self, name: str, capacitance: float = 0.0,
                 fixed: float | None = None, heat: QSpec = 0.0) -> "ThermalNetwork":
        if name in self.nodes:
            raise ValueError(f"ノード {name} は既にあります")
        self.nodes[name] = Node(name, capacitance, fixed, heat)
        return self

    def add_resistor(self, a: str, b: str, resistance: RSpec,
                     label: str = "") -> "ThermalNetwork":
        for n in (a, b):
            if n not in self.nodes:
                self.add_node(n)
        self.branches.append(Branch(a, b, resistance, label))
        return self

    def add_series(self, path: Sequence[str], resistances: Sequence[RSpec],
                   labels: Sequence[str] | None = None) -> "ThermalNetwork":
        """path = [n0, n1, n2, ...] を直列につなぐ糖衣。"""
        if len(path) != len(resistances) + 1:
            raise ValueError("path の長さは resistances より 1 だけ多い必要があります")
        for i, r in enumerate(resistances):
            lab = labels[i] if labels else ""
            self.add_resistor(path[i], path[i + 1], r, lab)
        return self

    def set_heat(self, name: str, heat: QSpec) -> "ThermalNetwork":
        self.nodes[name].heat = heat
        return self

    # ------------------------------------------------------------ 行列組み立て
    def _index(self) -> Dict[str, int]:
        return {name: i for i, name in enumerate(self.nodes)}

    def assemble(self, temps: np.ndarray, time: float = 0.0):
        """コンダクタンス行列 A [W/K] と熱源ベクトル b [W] を作る。

        A[i,i] = sum_j G_ij,  A[i,j] = -G_ij,  b[i] = Q_i
        温度固定ノードの行は後段(_apply_fixed)で単位行に置き換える。
        """
        idx = self._index()
        n = len(idx)
        a = np.zeros((n, n))
        b = np.zeros(n)
        for br in self.branches:
            i, j = idx[br.a], idx[br.b]
            g = 1.0 / br.value(temps[i], temps[j])
            a[i, i] += g
            a[j, j] += g
            a[i, j] -= g
            a[j, i] -= g
        for name, node in self.nodes.items():
            q = node.heat(time) if callable(node.heat) else node.heat
            b[idx[name]] += q
        return a, b

    def _fixed_mask(self) -> np.ndarray:
        return np.array([node.fixed is not None for node in self.nodes.values()])

    def _initial_vector(self, guess: float | Dict[str, float]) -> np.ndarray:
        idx = self._index()
        t = np.full(len(idx), guess if np.isscalar(guess) else 0.0, dtype=float)
        if isinstance(guess, dict):
            base = float(np.mean(list(guess.values()))) if guess else 20.0
            t[:] = base
            for k, v in guess.items():
                t[idx[k]] = v
        for name, node in self.nodes.items():
            if node.fixed is not None:
                t[idx[name]] = node.fixed
        return t

    def _apply_fixed(self, a: np.ndarray, b: np.ndarray) -> None:
        for i, node in enumerate(self.nodes.values()):
            if node.fixed is not None:
                a[i, :] = 0.0
                a[i, i] = 1.0
                b[i] = node.fixed

    def _is_nonlinear(self) -> bool:
        return any(callable(br.resistance) for br in self.branches)

    # ------------------------------------------------------------ 定常解
    def solve_steady(self, guess: float | Dict[str, float] = 20.0,
                     tol: float = 1e-8, max_iter: int = 200,
                     relaxation: float = 0.7) -> Solution:
        """定常解を解く。非線形抵抗があれば不動点反復(低速だが分かりやすい)。"""
        t = self._initial_vector(guess)
        iterations, residual = 0, 0.0
        for iterations in range(1, max_iter + 1):
            a, b = self.assemble(t)
            self._apply_fixed(a, b)
            t_new = np.linalg.solve(a, b)
            residual = float(np.max(np.abs(t_new - t)))
            t = t + relaxation * (t_new - t) if self._is_nonlinear() else t_new
            if residual < tol or not self._is_nonlinear():
                break
        temps = {name: float(t[i]) for i, name in enumerate(self.nodes)}
        return Solution(temps, self, iterations, residual)

    # ------------------------------------------------------------ 過渡解
    def solve_transient(self, t_end: float, dt: float,
                        initial: float | Dict[str, float] = 20.0,
                        theta: float = 1.0, store_every: int = 1,
                        picard_iter: int = 3) -> TransientSolution:
        """theta 法で時間積分する。

        theta = 1.0 : 後退オイラー(陰解法、無条件安定)  <- 既定
        theta = 0.5 : クランク・ニコルソン(2 次精度)
        theta = 0.0 : 前進オイラー(陽解法、dt < 2 R C で不安定化)

        熱容量 0 のノードは代数拘束なので、theta によらず陰的に扱う。
        """
        if not 0.0 <= theta <= 1.0:
            raise ValueError("theta は 0..1")
        idx = self._index()
        n = len(idx)
        caps = np.array([node.capacitance for node in self.nodes.values()], dtype=float)
        fixed = self._fixed_mask()
        algebraic = (caps <= 0.0) & (~fixed)
        if theta == 0.0 and algebraic.any():
            raise ValueError("陽解法では熱容量 0 のノードを扱えません(C を与えてください)")

        t = self._initial_vector(initial)
        nsteps = int(round(t_end / dt))
        times = [0.0]
        history = [t.copy()]

        a_old, b_old = self.assemble(t, 0.0)
        for step in range(1, nsteps + 1):
            time_new = step * dt
            t_new = t.copy()
            loops = picard_iter if self._is_nonlinear() else 1
            for _ in range(loops):
                a_new, b_new = self.assemble(t_new, time_new)
                lhs = np.diag(caps / dt) + theta * a_new
                rhs = (caps / dt) * t - (1.0 - theta) * (a_old @ t - b_old) + theta * b_new
                # 熱容量ゼロのノードは常に陰的な代数式 A T = b
                for i in np.where(algebraic)[0]:
                    lhs[i, :] = a_new[i, :]
                    rhs[i] = b_new[i]
                self._apply_fixed(lhs, rhs)
                t_new = np.linalg.solve(lhs, rhs)
            t = t_new
            a_old, b_old = self.assemble(t, time_new)
            if step % store_every == 0 or step == nsteps:
                times.append(time_new)
                history.append(t.copy())

        arr = np.array(history)
        temps = {name: arr[:, i] for i, name in enumerate(self.nodes)}
        return TransientSolution(np.array(times), temps, list(self.nodes))

    # ------------------------------------------------------------ 後処理
    def branch_flow(self, a: str, b: str, temps: Dict[str, float],
                    label: str | None = None) -> float:
        """枝 a->b の熱流 [W]。並列枝があるときは label で区別する。"""
        matches = [br for br in self.branches
                   if {br.a, br.b} == {a, b} and (label is None or br.label == label)]
        if not matches:
            raise KeyError(f"枝 {a}-{b} (label={label}) は存在しません")
        if len(matches) > 1:
            labels = [br.label or "(無名)" for br in matches]
            raise KeyError(f"枝 {a}-{b} は並列に {len(matches)} 本あります "
                           f"{labels}。label= を指定してください")
        br = matches[0]
        ta, tb = temps[br.a], temps[br.b]
        return (temps[a] - temps[b]) / br.value(ta, tb)

    def total_resistance(self, source: str, sink: str, power: float = 1.0) -> float:
        """source に power [W] を入れ sink を 0 degC 固定にしたときの合成熱抵抗。

        枝が線形な回路にのみ意味がある(非線形なら動作点ごとに変わる)。
        """
        saved = {name: (node.fixed, node.heat) for name, node in self.nodes.items()}
        try:
            for node in self.nodes.values():
                node.heat = 0.0
                if node.fixed is not None:
                    node.fixed = None
            self.nodes[sink].fixed = 0.0
            self.nodes[source].heat = power
            sol = self.solve_steady(guess=0.0)
            return sol[source] / power
        finally:
            for name, (fixed, heat) in saved.items():
                self.nodes[name].fixed = fixed
                self.nodes[name].heat = heat

    def describe(self) -> str:
        lines = [f"ThermalNetwork '{self.name}': "
                 f"{len(self.nodes)} ノード / {len(self.branches)} 枝"]
        for name, node in self.nodes.items():
            bits = []
            if node.fixed is not None:
                bits.append(f"T固定={node.fixed} degC")
            if node.capacitance:
                bits.append(f"C={node.capacitance:.4g} J/K")
            if node.heat:
                bits.append("Q=可変" if callable(node.heat) else f"Q={node.heat:.4g} W")
            lines.append(f"  node {name:<20s} {' '.join(bits)}")
        for br in self.branches:
            r = "f(T)" if callable(br.resistance) else f"{br.resistance:.5g}"
            lines.append(f"  branch {br.a:<14s}-{br.b:<14s} R={r} {br.label}")
        return "\n".join(lines)
