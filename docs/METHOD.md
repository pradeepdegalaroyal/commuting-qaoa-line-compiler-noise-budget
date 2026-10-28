# Method

## Problem and circuits

MaxCut on a connected Erdős–Rényi graph with n vertices and edge probability d (unit weights). The Quantum Approximate Optimization Algorithm (QAOA)
of depth p prepares |+>^n and applies p times the cost layer exp(-i gamma C), with C the cut operator, followed by the mixer exp(-i beta X) on every qubit.
Parameters are optimised once, **noiselessly**, by the Broyden–Fletcher–Goldfarb–Shanno (BFGS) method on a state-vector simulator (p = 1: 12 starts; p = 2: the interpolation-based
initialisation of Zhou et al., 2020, plus 3 random restarts) and then frozen. Every compiler receives the same parameters, so any difference in
the noisy result comes from compilation alone.

Device model: n qubits on a line (nearest-neighbour controlled-NOT, "cx"), native gates cx, rz, sx, x.

## Compilers (the thing under study)

| Name | Role | What it does |
|---|---|---|
| `qiskit_L1`, `qiskit_L3` | standard toolchain baseline | Qiskit `transpile` with SABRE (SWAP-based BidiREctional heuristic search) layout and routing, optimisation level 1 or 3 |
| `qiskit_L3_best16` | **primary baseline** | level 3 with 16 transpiler seeds, cheapest circuit kept (same compile-time search budget as the proposed router) |
| `swap_network` | standard quantum-literature baseline for dense cost layers | odd-even transposition network: every qubit pair meets once per layer; a pair with an edge uses a fused ZZ-SWAP (3 cx; "ZZ" is the two-qubit Pauli-Z interaction of one cut-operator term), others a bare SWAP |
| `ordered_route` | ablation step 0 | cost terms in a fixed order; the left qubit is swapped next to the right one before each term |
| `commute_order` | + free term order | cost terms commute, so any currently adjacent term is executed first; each swap is chosen to minimise the summed distance of the remaining terms |
| `commute_fuse` | + ZZ-SWAP fusion | a swap whose pair also carries a pending term is merged with it (ZZ then SWAP = 3 cx instead of 5) |
| `commute_fuse_multistart` | **proposed (full router)** | 16 randomised runs (random tie-breaks, random initial layout), cheapest by cx count then depth |

No compiler applies cx cancellation after routing, except what is explicit in its own rules (fusion), so the ablation steps are not blurred
by hidden clean-up. The final qubit permutation is tracked and absorbed into the read-out (the mixer is permutation invariant, so no restoring swaps).
Every compiled circuit is checked in the tests to reproduce the noiseless expectation value to 1e-7.

## Noise model and metrics

Depolarising channel after every cx with parameter lambda2 (swept 0 ... 0.02) and after every sx and x with lambda2 / 10; rz is virtual and noiseless.
Exact density-matrix simulation (Qiskit Aer), so there is no shot noise: the expected cut is computed from the exact output distribution.
Idle decoherence (energy relaxation and dephasing during waiting time) and read-out error are not modelled (see README limitations).

- **Normalised gain** g = (E[cut] − m/2) / (maxcut − m/2), with m the number of edges: 0 is uniform random guessing, 1 the optimum.
- **Compiled cost:** cx count and circuit depth of the native-gate circuit.
- **Break-even error:** the lambda2 at which the p = 2 gain falls to the p = 1 gain (linear interpolation on the grid; instances where it never happens are reported as right-censored).
- **Half-gain error:** the lambda2 at which the gain falls to half of its noiseless value.
- **Classical reference:** steepest-ascent one-flip local search from 200 random starts, same normalisation.

## Design

- Main study: n = 8; edge probability d in {0.3, 0.5, 0.8}; 5 graphs per density (seeds 0-4) = 15 instances; p in {1, 2}; 8 compilers; 10 noise levels.
- Ablation: the four-step router ladder above, evaluated at lambda2 = 0.01, p = 2.
- Noise / robustness sweep: lambda2 from 0 to 0.02, giving gain curves, break-even and half-gain errors.
- Scaling: n in {6, 8, 10}, d = 0.5, 5 graphs, p = 2, lambda2 in {0.005, 0.01}, three compilers.
- Statistics: mean ± sample standard deviation over graphs; paired differences between compilers on the same graph with a Wilcoxon signed-rank test (reported per density and overall).
- Hypotheses H1-H3 and their thresholds are in `RESEARCH_GAP.md`.

## Reproduction

`python experiments/run_all.py` runs everything (one command; `--quick` is a smoke run). All raw rows are written to `results/main_runs.csv` and `results/scaling_runs.csv`.
