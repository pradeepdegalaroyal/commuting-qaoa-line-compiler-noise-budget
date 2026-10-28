"""MaxCut instances, exact cut values, QAOA state-vector simulation and parameter optimisation."""
from __future__ import annotations

import networkx as nx
import numpy as np
from scipy.optimize import minimize

Edge = tuple[int, int]


def random_graph(n: int, density: float, seed: int) -> list[Edge]:
    """Connected Erdos-Renyi graph on ``n`` vertices with edge probability ``density``."""
    rng = np.random.default_rng(seed)
    while True:
        g = nx.gnp_random_graph(n, density, seed=int(rng.integers(2**31)))
        if nx.is_connected(g):
            return sorted((min(u, v), max(u, v)) for u, v in g.edges())


def cut_values(n: int, edges: list[Edge]) -> np.ndarray:
    """Cut value of every bitstring; bit ``i`` of the index is the side of vertex ``i``."""
    idx = np.arange(2**n)
    cut = np.zeros(2**n)
    for u, v in edges:
        cut += ((idx >> u) ^ (idx >> v)) & 1
    return cut


def normalised_gain(expected_cut: float, cut: np.ndarray, n_edges: int) -> float:
    """Map the expected cut to 0 (uniform random guessing, m/2) ... 1 (optimal cut)."""
    return (expected_cut - n_edges / 2) / (cut.max() - n_edges / 2)


def qaoa_state(n: int, cut: np.ndarray, gammas: np.ndarray, betas: np.ndarray) -> np.ndarray:
    """QAOA state exp(-i b X) exp(-i g C) ... |+>^n with C the cut operator (up to global phase)."""
    psi = np.full(2**n, 2 ** (-n / 2), dtype=complex)
    for g, b in zip(gammas, betas):
        psi = psi * np.exp(-1j * g * cut)
        c, s = np.cos(b), -1j * np.sin(b)
        psi = psi.reshape((2,) * n)
        for ax in range(n):
            a0 = np.take(psi, 0, axis=ax)
            a1 = np.take(psi, 1, axis=ax)
            psi = np.stack([c * a0 + s * a1, s * a0 + c * a1], axis=ax)
        psi = psi.reshape(-1)
    return psi


def expected_cut(n: int, cut: np.ndarray, params: np.ndarray) -> float:
    p = len(params) // 2
    psi = qaoa_state(n, cut, params[:p], params[p:])
    return float(np.real(np.vdot(psi, cut * psi)))


def optimise_qaoa(n: int, cut: np.ndarray, p: int, seed: int, prev: np.ndarray | None = None) -> np.ndarray:
    """Noiseless parameter optimisation; p>1 starts from the INTERP initialisation of the p-1 optimum
    (Zhou et al., 2020) plus a few random restarts. Returns [gammas..., betas...]."""
    rng = np.random.default_rng(seed)
    starts = []
    if p == 1:
        for g in np.linspace(0.2, 1.2, 4):
            for b in np.linspace(0.2, 0.8, 3):
                starts.append(np.array([g, b]))
    else:
        assert prev is not None
        q = p - 1
        gam, bet = prev[:q], prev[q:]

        def interp(x: np.ndarray) -> np.ndarray:
            return np.array(
                [((i) / q) * (x[i - 1] if i >= 1 else 0.0) + ((q - i) / q) * (x[i] if i < q else 0.0)
                 for i in range(p)]
            )

        starts.append(np.concatenate([interp(gam), interp(bet)]))
        for _ in range(3):
            starts.append(np.concatenate([rng.uniform(0.1, 1.2, p), rng.uniform(0.1, 0.8, p)]))
    best, best_val = None, -np.inf
    for x0 in starts:
        res = minimize(lambda x: -expected_cut(n, cut, x), x0, method="BFGS", options={"gtol": 1e-6})
        if -res.fun > best_val:
            best, best_val = res.x, -res.fun
    return best


def local_search_gain(n: int, cut: np.ndarray, n_edges: int, seed: int, starts: int = 200) -> float:
    """Classical baseline: steepest-ascent 1-flip local search from uniformly random starts.
    Returns the mean normalised gain (same scale as :func:`normalised_gain`) over the starts."""
    rng = np.random.default_rng(seed)
    flips = 1 << np.arange(n)
    vals = []
    for x in rng.integers(0, 2**n, size=starts):
        while True:
            nb = x ^ flips
            j = int(np.argmax(cut[nb]))
            if cut[nb[j]] <= cut[x]:
                break
            x = int(nb[j])
        vals.append(cut[x])
    return normalised_gain(float(np.mean(vals)), cut, n_edges)
