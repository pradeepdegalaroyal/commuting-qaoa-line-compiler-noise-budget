import numpy as np
import pytest

from qaoacomp.compilers import COMPILERS, compile_all
from qaoacomp.evaluate import noisy_expected_cut
from qaoacomp.problem import cut_values, expected_cut, optimise_qaoa, random_graph

N = 6


@pytest.fixture(scope="module")
def inst():
    edges = random_graph(N, 0.6, seed=3)
    cut = cut_values(N, edges)
    p1 = optimise_qaoa(N, cut, 1, 0)
    p2 = optimise_qaoa(N, cut, 2, 0, p1)
    return edges, cut, p2


def test_cut_values_triangle():
    cut = cut_values(3, [(0, 1), (1, 2), (0, 2)])
    assert cut.max() == 2 and cut[0] == 0


def test_optimiser_beats_random(inst):
    edges, cut, params = inst
    assert expected_cut(N, cut, params) > len(edges) / 2 + 0.3


@pytest.mark.parametrize("name", COMPILERS)
def test_noiseless_compilation_matches_ideal(inst, name):
    edges, cut, params = inst
    comp = compile_all(name, N, edges, params)
    assert noisy_expected_cut(comp, cut, 0.0) == pytest.approx(expected_cut(N, cut, params), abs=1e-7)


def test_only_native_gates_and_line_connectivity(inst):
    edges, _, params = inst
    for name in COMPILERS:
        comp = compile_all(name, N, edges, params)
        assert set(comp.circuit.count_ops()) <= {"cx", "rz", "sx", "x"}
        for ins in comp.circuit.data:
            if ins.operation.name == "cx":
                a, b = (comp.circuit.find_bit(q).index for q in ins.qubits)
                assert abs(a - b) == 1


def test_multistart_never_worse_than_single_start(inst):
    edges, _, params = inst
    single = compile_all("commute_fuse", N, edges, params).cx_count
    multi = compile_all("commute_fuse_multistart", N, edges, params).cx_count
    assert multi <= single


def test_fused_zz_swap_saves_two_cx():
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator

    from qaoacomp.compilers import _swap, _zz, _zz_swap

    a, b = QuantumCircuit(2), QuantumCircuit(2)
    _zz(a, 0, 1, 0.7)
    _swap(a, 0, 1)
    _zz_swap(b, 0, 1, 0.7)
    assert Operator(a).equiv(Operator(b)) and b.count_ops()["cx"] == 3 and a.count_ops()["cx"] == 5


def test_noise_reduces_expected_cut(inst):
    edges, cut, params = inst
    comp = compile_all("commute_fuse", N, edges, params)
    assert noisy_expected_cut(comp, cut, 0.02) < noisy_expected_cut(comp, cut, 0.0)
