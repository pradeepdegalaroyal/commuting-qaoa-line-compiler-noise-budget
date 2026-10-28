"""Compilers that map a p-level QAOA MaxCut circuit onto a line of qubits (native gates cx, rz, sx, x).

All compilers return a :class:`Compiled` object: the circuit on physical qubits plus the final
logical-to-physical map ``pos`` (the QAOA mixer is permutation invariant, so no restoring swaps are used).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import Optimize1qGatesDecomposition

from .problem import Edge

BASIS = ["cx", "rz", "sx", "x"]


@dataclass
class Compiled:
    """Compiled circuit and bookkeeping."""

    name: str
    circuit: QuantumCircuit
    pos: list[int]  # pos[logical vertex] = physical qubit holding it at the end
    cx_count: int
    depth: int


@dataclass(frozen=True)
class Policy:
    """Switches of the greedy commutation-aware router (the ablation ladder)."""

    free_order: bool = True  # execute whichever pending ZZ term is adjacent; order is free (terms commute)
    fuse: bool = True  # merge a SWAP with a pending ZZ term on the same pair into 3 cx
    restarts: int = 1  # randomised restarts (random tie-breaks and initial layout), keep the cheapest


POLICIES = {
    "ordered_route": Policy(free_order=False, fuse=False, restarts=1),
    "commute_order": Policy(free_order=True, fuse=False, restarts=1),
    "commute_fuse": Policy(free_order=True, fuse=True, restarts=1),
    "commute_fuse_multistart": Policy(free_order=True, fuse=True, restarts=16),
}


def _finalise(qc: QuantumCircuit, name: str, pos: list[int]) -> Compiled:
    """Translate to the native basis and merge single-qubit gates; two-qubit gates are never cancelled."""
    tqc = transpile(qc, basis_gates=BASIS, optimization_level=0)
    tqc = PassManager([Optimize1qGatesDecomposition(basis=["rz", "sx", "x"])]).run(tqc)
    return Compiled(name, tqc, pos, tqc.count_ops().get("cx", 0), tqc.depth())


def _zz(qc: QuantumCircuit, a: int, b: int, gamma: float) -> None:
    """exp(+i gamma/2 Z_a Z_b), the phase-separator term of the cut operator (2 cx)."""
    qc.cx(a, b)
    qc.rz(-gamma, b)
    qc.cx(a, b)


def _zz_swap(qc: QuantumCircuit, a: int, b: int, gamma: float) -> None:
    """ZZ term followed by SWAP(a, b) in 3 cx (the two middle cx cancel)."""
    qc.cx(a, b)
    qc.rz(-gamma, b)
    qc.cx(b, a)
    qc.cx(a, b)


def _swap(qc: QuantumCircuit, a: int, b: int) -> None:
    qc.cx(a, b)
    qc.cx(b, a)
    qc.cx(a, b)


def _mixer(qc: QuantumCircuit, n: int, beta: float) -> None:
    for q in range(n):
        qc.rx(2 * beta, q)


class _Layout:
    """Logical <-> physical map on the line with adjacent-swap updates."""

    def __init__(self, n: int, order: list[int]) -> None:
        self.at = list(order)  # at[physical] = logical
        self.pos = [0] * n
        for ph, lg in enumerate(self.at):
            self.pos[lg] = ph

    def swap(self, k: int) -> None:
        la, lb = self.at[k], self.at[k + 1]
        self.at[k], self.at[k + 1] = lb, la
        self.pos[la], self.pos[lb] = k + 1, k

    def dist(self, e: Edge) -> int:
        return abs(self.pos[e[0]] - self.pos[e[1]])

    def key(self, k: int) -> Edge:
        a, b = self.at[k], self.at[k + 1]
        return (min(a, b), max(a, b))


def _route_layer_ordered(qc: QuantumCircuit, lay: _Layout, edges: list[Edge], g: float) -> None:
    """Fixed term order; for each term move the left endpoint next to the right one, then apply ZZ."""
    for e in edges:
        while lay.dist(e) > 1:
            k = min(lay.pos[e[0]], lay.pos[e[1]])
            _swap(qc, k, k + 1)
            lay.swap(k)
        _zz(qc, lay.pos[e[0]], lay.pos[e[1]], g)


def _route_layer_free(qc: QuantumCircuit, lay: _Layout, edges: list[Edge], g: float, policy: Policy,
                      rng: np.random.Generator) -> None:
    """Commutation-aware greedy routing: adjacent terms run first, swaps minimise the total remaining
    distance, and (if ``policy.fuse``) a swap that carries a pending term is fused into 3 cx."""
    n = len(lay.at)
    pending = list(edges)
    for _ in range(10 * n * n + 100):
        adj = [e for e in pending if lay.dist(e) == 1]
        if len(adj) == len(pending):
            for e in adj:
                _zz(qc, lay.pos[e[0]], lay.pos[e[1]], g)
            return
        phi0 = sum(lay.dist(e) for e in pending)
        best: list[int] = []
        best_score = -np.inf
        for k in range(n - 1):
            carries = lay.key(k) in pending
            if carries and not policy.fuse:
                continue  # without fusion such a pair is simply executed as an ordinary ZZ term
            lay.swap(k)
            gain = phi0 - sum(lay.dist(e) for e in pending)
            lay.swap(k)
            score = gain + (0.5 if carries else 0.0)
            if score > best_score + 1e-9:
                best, best_score = [k], score
            elif abs(score - best_score) < 1e-9:
                best.append(k)
        k = int(rng.choice(best)) if best else -1
        carried = lay.key(k) if k >= 0 and policy.fuse and lay.key(k) in pending else None
        if k >= 0 and best_score <= 0 and carried is None:
            k = -1  # no swap improves anything
        for e in adj:
            if e != carried:
                _zz(qc, lay.pos[e[0]], lay.pos[e[1]], g)
        pending = [e for e in pending if e not in adj or e == carried]
        if k < 0:  # step along the path of the closest pending term
            e = min(pending, key=lay.dist)
            k = min(lay.pos[e[0]], lay.pos[e[1]])
            _swap(qc, k, k + 1)
        elif carried is not None:
            _zz_swap(qc, k, k + 1, g)
            pending.remove(carried)
        else:
            _swap(qc, k, k + 1)
        lay.swap(k)
        if not pending:
            return
    raise RuntimeError("router did not terminate")


def _route_circuit(n: int, edges: list[Edge], params: np.ndarray, policy: Policy, rng: np.random.Generator,
                   random_layout: bool) -> tuple[QuantumCircuit, list[int]]:
    p = len(params) // 2
    order = [int(x) for x in rng.permutation(n)] if random_layout else list(range(n))
    lay = _Layout(n, order)
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for layer in range(p):
        g = float(params[layer])
        if policy.free_order:
            _route_layer_free(qc, lay, edges, g, policy, rng)
        else:
            _route_layer_ordered(qc, lay, edges, g)
        _mixer(qc, n, float(params[p + layer]))
    return qc, lay.pos


def compile_commute(n: int, edges: list[Edge], params: np.ndarray, policy_name: str, seed: int = 0) -> Compiled:
    """Greedy commutation-aware router on a line (see :data:`POLICIES` for the ablation ladder)."""
    policy = POLICIES[policy_name]
    rng = np.random.default_rng(seed)
    best: Compiled | None = None
    for r in range(policy.restarts):
        qc, pos = _route_circuit(n, edges, params, policy, rng, random_layout=r > 0)
        comp = _finalise(qc, policy_name, pos)
        if best is None or (comp.cx_count, comp.depth) < (best.cx_count, best.depth):
            best = comp
    assert best is not None
    return best


def compile_swap_network(n: int, edges: list[Edge], params: np.ndarray) -> Compiled:
    """Odd-even transposition swap network: every vertex pair meets once per layer (fixed cost n(n-1)/2 gates)."""
    p = len(params) // 2
    eset = set(edges)
    lay = _Layout(n, list(range(n)))
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for layer in range(p):
        g = float(params[layer])
        for rnd in range(n):
            last = layer == p - 1 and rnd == n - 1  # final swaps would only relabel qubits: skip them
            for k in range(rnd % 2, n - 1, 2):
                has = lay.key(k) in eset
                if last:
                    if has:
                        _zz(qc, k, k + 1, g)
                    continue
                if has:
                    _zz_swap(qc, k, k + 1, g)
                else:
                    _swap(qc, k, k + 1)
                lay.swap(k)
        _mixer(qc, n, float(params[p + layer]))
    return _finalise(qc, "swap_network", lay.pos)


def compile_qiskit(n: int, edges: list[Edge], params: np.ndarray, level: int, seed: int = 7,
                   restarts: int = 1) -> Compiled:
    """Standard toolchain baseline: Qiskit transpile with SABRE layout/routing on a line.

    With ``restarts > 1`` the transpiler seed is varied and the circuit with the fewest cx (then depth) is kept,
    which gives the baseline the same compile-time search budget as the multi-start router."""
    p = len(params) // 2
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for layer in range(p):
        for u, v in edges:
            qc.rzz(-float(params[layer]), u, v)
        _mixer(qc, n, float(params[p + layer]))
    best: Compiled | None = None
    for r in range(restarts):
        tqc = transpile(qc, coupling_map=CouplingMap.from_line(n), basis_gates=BASIS,
                        optimization_level=level, seed_transpiler=seed + r)
        pos = [int(x) for x in tqc.layout.final_index_layout(filter_ancillas=True)]
        comp = Compiled(f"qiskit_L{level}", tqc, pos, tqc.count_ops().get("cx", 0), tqc.depth())
        if best is None or (comp.cx_count, comp.depth) < (best.cx_count, best.depth):
            best = comp
    assert best is not None
    return best


COMPILERS = ["qiskit_L1", "qiskit_L3", "qiskit_L3_best16", "swap_network", *POLICIES]


def compile_all(name: str, n: int, edges: list[Edge], params: np.ndarray, seed: int = 0) -> Compiled:
    """Dispatch by compiler name (see :data:`COMPILERS`)."""
    if name in POLICIES:
        return compile_commute(n, edges, params, name, seed)
    if name == "swap_network":
        return compile_swap_network(n, edges, params)
    if name == "qiskit_L3_best16":
        return compile_qiskit(n, edges, params, 3, restarts=16)
    if name.startswith("qiskit_L"):
        return compile_qiskit(n, edges, params, int(name[-1]))
    raise ValueError(name)
