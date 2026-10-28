# Research gap: how much of the noise budget of the Quantum Approximate Optimization Algorithm (QAOA) is a compiler decision?

## State of the art (from the supplied literature)

The Noisy Intermediate-Scale Quantum (NISQ) setting is defined by gate noise that caps the size of circuits that run reliably
(Preskill, 2018, doi:10.22331/q-2018-08-06-79). Variational quantum algorithms (VQAs) are the leading answer to this
constraint, but their trainability, accuracy and efficiency remain open (Cerezo et al., 2021, doi:10.1038/s42254-021-00348-9;
Bharti et al., 2022, doi:10.1103/RevModPhys.94.015004). For QAOA on MaxCut, Zhou et al. (2020, doi:10.1103/PhysRevX.10.021067)
give the standard recipe for choosing parameters at depth p (interpolation-based initialisation, "INTERP") and a resource
analysis based on measurement fluctuations. Noise erodes the benefit of
QAOA variants in recent work: the structured mixers of Lei et al. (2027, doi:10.1016/j.trc.2026.106047) lose their advantage under
hardware-inspired noise "in part because [they] introduce additional two-qubit operations", and Qiskit benchmarks find QAOA more
noise-robust than deep Grover oracles (Lavanya et al., 2027, doi:10.1007/978-3-032-30960-0_42). Fidelity decays with the number of
two-qubit gates (Boixo et al., 2018, doi:10.1038/s41567-018-0124-x), so the gate count of the *compiled* circuit is the quantity that
sets the noise exposure.

## Related work

| Paper | Year | Role | What it does | What it does NOT cover |
|---|---|---|---|---|
| Preskill, doi:10.22331/q-2018-08-06-79 | 2018 | landmark | Defines the NISQ era; gate noise limits circuit size | Any concrete compile-level lever for a given algorithm |
| Cerezo et al., doi:10.1038/s42254-021-00348-9 | 2021 | landmark | Reviews VQAs; names trainability, accuracy, efficiency challenges | Quantifies how much compiled depth varies for one VQA |
| Bharti et al., doi:10.1103/RevModPhys.94.015004 | 2022 | landmark | Surveys NISQ algorithms, benchmarking and software tools | Ablation of compiler components on a QAOA circuit |
| Zhou et al., doi:10.1103/PhysRevX.10.021067 | 2020 | landmark | QAOA on MaxCut, INTERP parameter initialisation, measurement-based resource analysis | Gate noise on a connectivity-restricted, compiled circuit |
| Boixo et al., doi:10.1038/s41567-018-0124-x | 2018 | landmark | Cross-entropy benchmarking; fidelity extrapolated from gate error rates | Optimisation workloads such as QAOA |
| Huang et al., doi:10.1038/s41467-021-22539-9 | 2021 | landmark | Shows classical learners can match quantum models; a protocol for fair advantage claims | Combinatorial optimisation; compilation |
| Terhal, doi:10.1103/RevModPhys.87.307 | 2015 | landmark | Reviews stabiliser codes, thresholds and the route to fault tolerance | The cost of compiled variational circuits before error correction |
| Cheng et al., doi:10.1007/s11467-022-1249-z | 2023 | review | Reviews NISQ hardware progress and next-stage challenges | Algorithm-level compile cost on a given topology |
| Lei et al., doi:10.1016/j.trc.2026.106047 | 2027 | recent | Constraint-aware initialisation and hybrid XY-X mixer for vehicle routing; advantage shrinks under noise | How the compiler's routing of the added two-qubit gates changes that result |
| Lavanya et al., doi:10.1007/978-3-032-30960-0_42 | 2027 | recent | Compares Grover variants and QAOA on knapsack with and without noise | Compile strategy or device connectivity |
| Dahi et al., doi:10.1007/978-3-032-36226-1_21 | 2027 | recent | Machine-learning-assisted multi-objective QAOA; quantum-classical interaction overhead is the bottleneck | Gate-level compiled cost of one QAOA layer |
| Ahmed and Mäkitalo, doi:10.1007/978-3-032-39766-9_6 | 2027 | recent | Runtime telemetry to detect shot, readout and depolarising disturbances in QAOA workflows | A compile-time lever that lowers exposure to those disturbances |
| Kashvi and Verma, doi:10.1007/978-3-032-30957-0_10 | 2027 | recent | Single-qubit geometric MaxCut model motivated by QAOA depth and decoherence limits | Reducing depth of the multi-qubit QAOA itself; noisy comparison |
| Alwan et al., doi:10.1007/978-3-032-30497-1_18 | 2027 | adjacent | Trust-region damped quantum natural gradient for redundant parameterisations | Circuit-level gate overhead; it improves the optimiser, not the circuit |

## The gap

The recent QAOA papers above, and the reviews they sit in, attribute the loss of QAOA quality to "extra two-qubit operations" and "circuit
depth and decoherence", but they treat the gate count as a property of the algorithm. For QAOA it is largely a compiler outcome: the
cost-layer terms are mutually commuting, so their order, their routing on a line-connected device, and the fusion of a SWAP with a
neighbouring interaction are all free choices. The default toolchain (Qiskit with SABRE, the SWAP-based BidiREctional heuristic search
layout and routing pass) and the standard swap network do not exploit all of them. No paper in the supplied set measures how much of the
noise budget (the error rate where p = 2 stops beating p = 1, or where solution quality halves) is determined by these compile choices.

**Novelty in one sentence:** a component-by-component ablation of a commutation-aware line router for the QAOA cost layer, measured with exact
noisy density-matrix simulation against Qiskit and the swap network, reporting the shift of the NISQ noise budget rather than only gate counts.

## Falsifiable hypotheses and success metrics

Scale: simulation only, n = 8 qubits (scaling check n = 6, 8, 10), Erdős–Rényi MaxCut graphs, depolarising noise on every two-qubit gate (exact density matrix, no shot noise). Thresholds were fixed before the full run (a one-graph pilot was used only to debug the code).

- **H1 (gate count).** At equal compile-time search (best of 16), the full router needs at least 15 % fewer controlled-NOT (CNOT) gates than Qiskit level 3 on average over the 15 graphs at p = 2.
- **H2 (quality).** At two-qubit error 0.01 the paired mean gain in normalised cut quality (0 = random guess, 1 = optimum) is at least 0.02 higher than Qiskit level 3 (best of 16).
- **H3 (noise budget).** The error rate at which p = 2 stops beating p = 1 moves up by at least 10 % relative to Qiskit level 3 (best of 16).

Each hypothesis is reported as supported or not supported, per density, with paired statistics. A classical one-flip local search is reported as the reference solver; the study makes no claim of quantum advantage.
