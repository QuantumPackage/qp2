# PT2 on native OpenMP: benchmark and reproducibility contract

Benchmark of the in-process OpenMP `pt2_stochastic` (the only PT2 path
in this branch; the legacy distributed collector/slave is removed) on the
benzene fixture (Cholesky-reduced 2e integrals, 114 MOs, ~26 440
determinants, 24 comb teeth).

- Machines: Ampere compute node `cnix6-1` (Intel Xeon 6730P, 48 cores,
  gfortran 15, system OpenBLAS, `-O2 -mavx2`); Apple M1 macmini,
  8 cores (gfortran, Accelerate).
- Fixture: `benzene.ezfio` (the same wavefunction on both machines).
- Driver: launch contract below; two repeats per point.

## Launch contract (bit-reproducibility)

```
OMP_NUM_THREADS=1  QP_NTHREADS_PT2=<walk team size>
```

`rPT2` is then bit-identical for any walk team size and any repeat, for
a given seed (verified: Xeon 48-core, team 1/8/24/48, all repeats;
M1, team 1/3/8/16, all repeats — tables below). The walk stays
parallel: the window loop carries an explicit
`NUM_THREADS(nproc_target)` clause, and `nproc_target` comes from
`QP_NTHREADS_PT2` (memory-budget capped), independent of the OMP
default.

Why the OMP default must be 1: the one-electron AO->MO transform
(two GEMMs with K = ao_num) is dispatched by OpenBLAS across a thread
pool whose size follows the OMP default team, not
`OPENBLAS_NUM_THREADS`. The K reduction of such small-K GEMMs is split
across that pool, so the bits of `mo_one_e_integrals` track the OMP
default (stable at 1, a different stable value at 8, per-run scatter
at 48 under contention). one_e feeds the Fock matrix, Hii and the E0
denominator, so the drift reaches every walk weight and rPT2 at the
~1e-10 level. Every other table build (jj, jj_exchange, the
MO-integrals cache, the Cholesky transform, H|u>, Hii) was measured
bit-stable across thread counts, but the contract makes them serial
regardless of BLAS vendor.

Cost: the integral-table build runs once per run, serially — about
30 s of wall for this fixture, independent of the walk team size.

Guardrail: the pt2 binary warns at start (C constructor in
`src/utils/c_funcs.c`) when launched without `OMP_NUM_THREADS=1`.

## Performance

### Ampere, Xeon 6730P, 48 cores (team = walk threads)

| team | wall (s) | CPU (s) | rPT2 |
|---|---|---|---|
| 1 | 833.1 / 835.3 | 831.9 / 834.0 | -0.30532354477010087 |
| 8 | 161.8 / 161.5 | 1048.5 / 1046.3 | -0.30532354477010087 |
| 24 | 77.3 / 75.8 | 1038.2 / 1008.0 | -0.30532354477010087 |
| 48 | 73.3 / 72.9 | 1815.8 / 1791.4 | -0.30532354477010087 |

Two repeats per row (the two values in each cell); **rPT2 is
bit-identical on every row** — the reproducibility contract holds on
the 48-core node. Speedup 1→48: 11.4×.

Reading the CPU column: the parallel walk's CPU is 2.2× the serial
CPU at 48 threads (1791–1816 s vs 832 s). This per-position CPU
growth at high thread counts is a property of the gather-bound walk
kernel under contention (measured identically on the pre-contract
launches and on the earlier binary); it is what bounds the
efficiency (24% at 48 threads), not the OpenMP machinery. The
8-thread row is a clean re-run (one earlier attempt hit a noisy
neighbor and ran +42% in CPU; discarded).

Against the pre-contract launch on the same node (older binary,
`OMP_NUM_THREADS=T` + `QP_NTHREADS_PT2=T`): 48 threads, 57.3 s wall
→ 72.9 s now: the +16 s is the one-time serial table build. At 24
threads the contract is *faster* than the old launch (75.8 s vs
89.8 s): under the old launch the walk team shared the cores with a
24-thread BLAS pool (the pool tracked the OMP default), which both
slowed the walk and broke the bit-reproducibility.

### Apple M1, 8 cores

| team | wall (s) | CPU (s) | rPT2 |
|---|---|---|---|
| 1 | 703.4 / 705.6 | 702.5 / 705.1 | -0.30532326656917652 |
| 3 | 298.2 / 323.2 | 842.6 / 918.0 | -0.30532326656917652 |
| 8 | 217.6 / 231.7 | 1524.0 / 1615.7 | -0.30532326656917652 |
| 16 *(SMT)* | 268.1 / 242.8 | 1900.3 / 1701.8 | -0.30532326656917652 |

Two repeats per row; rPT2 bit-identical on every row. The absolute
value differs from the Xeon number above — that is expected: the
Cholesky integral threshold is host-dependent, so the walk inputs
differ slightly between machines; cross-machine rPT2 comparison is
meaningless, within-machine constancy is the contract.

T=16 oversubscribes the 8 physical cores via SMT and is *slower*
than T=8 (243–268 s vs 218–232 s), with more CPU: SMT is a net loss
for the gather-bound walk, as on the P/E-hybrid laptop. Use one
thread per physical core.

### Dell laptop (i7-13850HX, 8 P + 12 E cores, 28 logical)

Findings from the 11-09/10 sweeps (earlier binary, pre-contract
launch): use all 28 logical CPUs without binding; E cores and SMT do
not help (20 no-SMT threads beat 16 pure-SMT; finer task granularity
is what matters). The contract adds the one-time serial table build.

## Reproduce

```bash
# Ampere (login node build, compute node run):
export QP_ROOT=~/qp2-pt2; source $QP_ROOT/quantum_package.rc
ninja -j 1 -f config/build.ninja $QP_ROOT/src/fci/pt2   # login node
# compute node, the launch contract:
OMP_NUM_THREADS=1 QP_NTHREADS_PT2=48 $QP_ROOT/src/fci/pt2 benzene.ezfio

# macmini (local build):
cd ~/qp2-omp-pt2
ninja -f config/build.ninja $PWD/src/fci/pt2
OMP_NUM_THREADS=1 QP_NTHREADS_PT2=8 src/fci/pt2 benzene.ezfio
```
