"""Experiment drivers: main compiler comparison (n = 8) and qubit-count scaling."""
from __future__ import annotations

import time

import pandas as pd

from .compilers import COMPILERS, compile_all
from .evaluate import noisy_expected_cut
from .problem import (cut_values, expected_cut, local_search_gain, normalised_gain, optimise_qaoa,
                      random_graph)

FULL = {
    "n": 8,
    "densities": [0.3, 0.5, 0.8],
    "n_graphs": 5,
    "lams": [0.0, 0.0005, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02],
    "compilers": COMPILERS,
    "scaling_n": [6, 8, 10],
    "scaling_graphs": 5,
    "scaling_lams": [0.005, 0.01],
    "scaling_compilers": ["qiskit_L3_best16", "swap_network", "commute_fuse_multistart"],
}
QUICK = {
    "n": 6,
    "densities": [0.5, 0.8],
    "n_graphs": 2,
    "lams": [0.0, 0.005, 0.01],
    "compilers": COMPILERS,
    "scaling_n": [6, 7],
    "scaling_graphs": 1,
    "scaling_lams": [0.01],
    "scaling_compilers": ["qiskit_L3_best16", "swap_network", "commute_fuse_multistart"],
}


def _instance(n: int, density: float, seed: int) -> dict:
    edges = random_graph(n, density, seed)
    cut = cut_values(n, edges)
    p1 = optimise_qaoa(n, cut, 1, seed)
    p2 = optimise_qaoa(n, cut, 2, seed, p1)
    return {"edges": edges, "cut": cut, "params": {1: p1, 2: p2}}


def _rows(n: int, density: float, seed: int, compilers: list[str], lams: list[float], ps=(1, 2)) -> list[dict]:
    inst = _instance(n, density, seed)
    edges, cut = inst["edges"], inst["cut"]
    m = len(edges)
    base = {"n": n, "density": density, "seed": seed, "m": m, "maxcut": float(cut.max()),
            "classical_local_search_gain": local_search_gain(n, cut, m, seed)}
    out = []
    for p in ps:
        params = inst["params"][p]
        ideal = normalised_gain(expected_cut(n, cut, params), cut, m)
        for name in compilers:
            comp = compile_all(name, n, edges, params, seed)
            for lam in lams:
                val = noisy_expected_cut(comp, cut, lam)
                out.append({**base, "p": p, "compiler": name, "cx": comp.cx_count, "depth": comp.depth,
                            "lam2": lam, "gain": normalised_gain(val, cut, m), "ideal_gain": ideal})
    return out


def run_main(cfg: dict, log=print) -> pd.DataFrame:
    rows = []
    t0 = time.time()
    for d in cfg["densities"]:
        for s in range(cfg["n_graphs"]):
            rows += _rows(cfg["n"], d, s, cfg["compilers"], cfg["lams"])
            log(f"  main: n={cfg['n']} density={d} seed={s} done ({time.time() - t0:.0f}s)")
    return pd.DataFrame(rows)


def run_scaling(cfg: dict, log=print) -> pd.DataFrame:
    rows = []
    t0 = time.time()
    for n in cfg["scaling_n"]:
        for s in range(cfg["scaling_graphs"]):
            rows += _rows(n, 0.5, s, cfg["scaling_compilers"], cfg["scaling_lams"], ps=(2,))
            log(f"  scaling: n={n} seed={s} done ({time.time() - t0:.0f}s)")
    return pd.DataFrame(rows)
