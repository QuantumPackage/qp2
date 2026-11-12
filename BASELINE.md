# Phase 0 baseline (pre-removal, `~/qp3` @ `1bbbf60e`)

Machine: Mac Mini `lcpqsv225` (8 cores / 8 GiB), homebrew GCC, 2026-11-12.
Driver: `baseline.sh` (each run rc=0; raw logs in `baseline/`).

These are the pre-removal references. After Phases 1-2 the same runs must
reproduce these values, except:

- CIPSI-det (and anything that goes through `ZMQ_selection`): the
  selection is ported to `omp_selection` (deterministic top-N under a
  total order instead of arrival-order buffer merging), so boundary
  determinants at the buffer threshold may differ at the tolerance
  level; the converged energy must match to ~1e-12 (see spec risk R1).
- Everything else (SCF, CISD, MRCI, stochastic PT2, the fci/pt2
  deterministic contract) must be bit-identical.

## SCF

| run | system / basis | converged energy (Ha) |
|---|---|---|
| scf-h2 | H2, sto-3g, Dinfh | -1.116759307489386 |
| scf-h2o | H2O, 6-31g*, Cinfv | -76.0091080303 (SCF iter 15) |

H2 MO energies: -0.5785538590 (Sigma_g+), 0.6711434839 (Sigma_u+).

## CIPSI (LiH, sto-3g, n_det_max 2000)

| run | path | final CI energy | final iteration rPT2 | E+rPT2 |
|---|---|---|---|---|
| cipsi-stoch | `qp run fci` (do_pt2 True → pt2_stochastic, already OMP) | -7.882321608808661 | -2.4408812275700665E-05 | -7.8823239997004180 |
| cipsi-det | do_pt2 False → **`ZMQ_selection` (the code being replaced)** | -7.882321608808661 | -2.4408812275700620E-05 | -7.8823239997004180 |

Both converge to the same 38-determinant FCI space (N_det = 38 = FCI on
LiH/sto-3g); the two rPT2 values already differ in the last digit
pre-removal (stochastic vs deterministic selection).

## CISD / MRCI (LiH, sto-3g)

| run | energy (Ha) | N_det |
|---|---|---|
| cisd | -7.8823110003069159 (CISD); -7.8828382998166111 (CISD+Q) | 93 |
| mrci | -7.882321608808661 | 38 (FCI) |

## Deterministic PT2 contract (benzene.ezfio fixture)

`OMP_NUM_THREADS=1 QP_NTHREADS_PT2=1` vs `QP_NTHREADS_PT2=8` on
`src/fci/pt2 ~/OpenCode/QP2/benzene.ezfio`:

| run | rPT2 | E+rPT2 |
|---|---|---|
| pt2-team1 | -0.30532326708442498 (± 5.9317597705745831E-04) | -231.50692258686118 |
| pt2-team8 | -0.30532326708442498 (± 5.9317597705745831E-04) | -231.50692258686118 |

**Bit-identical across thread counts** (the fci/pt2 BENCHMARK contract
holds on this machine).

## Dressing

Module dropped (spec decision 4, 2026-11-12: "nobody uses it"); no
baseline run.

## Phase 1+2 replay (post-removal worktree) — 2026-11-12

Same commands as `baseline.sh` on `~/qp3-nozmq` (`replay.sh`, logs in
`replay/`), plus a strict-FCI anchor run on **both** trees
(`fci-ref.sh`, `fci-ref.log`).

| run | baseline (pre-removal) | replay (post-removal) | status |
|---|---|---|---|
| scf-h2 | -1.116759307489386 | same (table 10-digit: -1.1167593075, MO energies identical) | bit-identical |
| scf-h2o | -76.0091080303 (iter 15) | -76.0091080303 (iter 15) | bit-identical |
| cipsi-stoch | E -7.882321608808661, 38 dets, rPT2 -2.4409e-5 | E -7.882323850069628, 54 dets, rPT2 -1.1951e-5 | **differs — R1, see below** |
| cipsi-det | E -7.882321608808661, 38 dets | E -7.882323850069628, 54 dets (27 with weight) | **differs — R1** |
| mrci | -7.882321608808661 (38 dets) | -7.882323850069628 (27 dets) | **differs — R1** |
| cisd | -7.8823110003069159 / +Q -7.8828382998166111, 93 dets | identical | bit-identical |
| pt2 team1/team8 | -0.30532326708442498, E+rPT2 -231.50692258686118 | identical, team1 == team8 | bit-identical |

### FCI anchor (strict CIPSI, `pt2_max 1e-8`, `n_det_max 300`)

| tree | N_det | final CI energy |
|---|---|---|
| pristine `~/qp3` (ZMQ selection) | 69 | -7.882324392913223 |
| worktree `~/qp3-nozmq` (omp_selection) | 69 | -7.882324392913223 |

**Bit-identical across the removed-code boundary** — the strongest
verification of the selection port: at strict tolerance both
implementations converge to the same 69-determinant space and energy.
Corollum: the baseline table's "N_det = 38 = FCI" was wrong; 38 is where
the *default* `pt2_max` stops the run, not the exact space (69 dets).

### R1 resolution (accepted 2026-11-12)

The old `ZMQ_selection` was **not a fixed point**:

1. Each worker's local selection buffer ratchets its threshold
   (`buf%mini`) across tasks and is *never reset* after a push
   (`run_selection_slave` sets only `buf%cur = 0`, `selection_buffer.irp.f`
   line 100; `select_connected` filters on `w <= buf%mini`). Task k+1 on a
   worker was therefore filtered by task k's N-th weight.
2. Task-to-worker assignment was a race over ZMQ inproc pulls
   (timing-dependent), so the ratcheted filter saw a different task order
   on every run.
3. The collector merged worker batches in arrival order.

Run-to-run probe (`det-check.sh`, default settings, `det-check.log`):

| tree | run 1 | run 2 | run 3 |
|---|---|---|---|
| pristine (ZMQ) | 50 dets, -7.882323476132574 | 38 dets, -7.882321608808661 | 38 dets, -7.882321608808661 |
| worktree (omp) | 54 dets, -7.882323850069628 | 54 dets, -7.882323850069628 | (replay run: 54 dets, same E) |

The baseline's own 38-det value is not even reproducible on the pristine
tree (a re-run gave 50 dets, E -7.882323476132574). The worktree is
deterministic (3/3 bit-identical). Of the four default-space samples, the
new one is closest to the strict anchor.

Consequences, all confirmed:

- The baseline's 38-det space was *one sample* of a run-dependent
  selection, so the spec's R1 expectation ("converged energy matches to
  ~1e-12") rested on a false premise. The honest contract is: each tree
  deterministic within itself, and agreement with the pre-removal release
  only at strict tolerance (proved above: 69 dets, bit-identical).
- The new `omp_selection` resets the per-task threshold (`mini = 0`) and
  merges each task's full local top-N into a shared heap under a named
  critical, so it offers strictly more candidates and is deterministic.
  Its default run (54 dets) lands **closer** to the strict anchor
  (ΔE 5.4e-7 vs the baseline's 2.8e-6) — the old ratcheting dropped
  candidates systematically.

For RELEASE_NOTES: "CIPSI/MRCI selected spaces at the default threshold
differ from the pre-removal release (run-dependent arrival-order merge
replaced by a deterministic selection); on LiH/sto-3g the converged
energy shifts -2.24 uHa toward the strict-FCI anchor, and strict
(`pt2_max 1e-8`) runs of both releases converge to the same 69-determinant
space, E = -7.882324392913223."
