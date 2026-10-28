"""Exact noisy evaluation of compiled circuits (Aer density matrix, depolarising gate noise)."""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

from .compilers import Compiled


@lru_cache(maxsize=None)
def noise_model(lam2: float, ratio_1q: float = 0.1) -> NoiseModel | None:
    """Depolarising channel after every cx (parameter lam2) and every sx/x (lam2 * ratio_1q); rz is virtual."""
    if lam2 <= 0:
        return None
    nm = NoiseModel(basis_gates=["cx", "rz", "sx", "x"])
    nm.add_all_qubit_quantum_error(depolarizing_error(lam2, 2), ["cx"])
    nm.add_all_qubit_quantum_error(depolarizing_error(lam2 * ratio_1q, 1), ["sx", "x"])
    return nm


def logical_cut_of_physical(n: int, cut: np.ndarray, pos: list[int]) -> np.ndarray:
    """Cut value of each *physical* bitstring given the final logical->physical map."""
    idx = np.arange(2**n)
    x = np.zeros_like(idx)
    for i in range(n):
        x |= ((idx >> pos[i]) & 1) << i
    return cut[x]


def noisy_expected_cut(comp: Compiled, cut: np.ndarray, lam2: float) -> float:
    """Exact expected cut of the compiled circuit under depolarising noise of strength ``lam2``."""
    n = comp.circuit.num_qubits
    qc = comp.circuit.copy()
    qc.save_probabilities()
    sim = AerSimulator(method="density_matrix", noise_model=noise_model(lam2), max_parallel_threads=6)
    probs = np.asarray(sim.run(qc).result().data()["probabilities"])
    return float(probs @ logical_cut_of_physical(n, cut, comp.pos))
