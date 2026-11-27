# Phase 0 — baseline + instrumentation of the Cholesky H|Ψ⟩ apply (Davidson/CIPSI)

Phase 0 (R0) of `~/OpenCode/tmp/qp3-hpsi-PLAN.md`, the merged dual-brain
report at `~/OpenCode/tmp/qp3-davidson-hpsi-report.md`.  Measurement phase:
no numerical change.  All instrumentation is env-gated
(`QP_DVDSN_STATS=1`); a run with the counters on must be bit-identical to
the same run with them off (ULP-GATE item 5).

## Where this stands (25 November 2026)

_(filled at the end of the Phase-0 campaign)_

## Machine, tree, build

* **Mac mini `lcpqsv225`** (8 cores / 8 GiB), homebrew gfortran,
  `-Ofast -march=native`, **Accelerate** BLAS (`-llapack -lblas`).
* **dell unreachable** (`ssh dell`: host does not resolve) → every number
  below is a mini number.  The plan's 114/264/production-size runs were
  specified for the dell (28 cores / 31 GiB); they are run here as fits
  and flagged.  The 28-thread anchor run oversubscribes the 8-core box;
  its wall time is noise, its bits are the gate.  The mini also carries a
  background load (audiomxd/configd, load average ~21) — wall times are
  pessimistic, ratios and bits are what Phase 0 consumes.
* Tree `c2d4ac2b` (branch `qp3-hpsi-opt`) plus the uncommitted Phase-0
  instrumentation: `src/utils/davidson_apply_stats.f90` (new) and hooks in
  `src/davidson/u0_hs2_u0.irp.f`,
  `src/davidson/diagonalization_hs2_dressed.irp.f`,
  `src/determinants/single_excitations.irp.f`,
  `src/mo_two_e_ints/map_integrals.irp.f`.  Nothing committed (standing
  rule).
* Build: `ninja -f config/build.ninja` after regenerating the graph —
  **a new plain `.f90` in `src/utils/` does not join the graph by itself**:
  `config/build.ninja` carries a hard-coded per-binary `SRC` list, and the
  pre-existing file predates the module, so the first instrumented build
  failed with `Cannot open module file`.  Regeneration is
  `python3 scripts/compilation/qp_create_ninja create --development
  config/gfortran_macos.cfg` (the `CONFIG_FILE` recorded in
  `config/build.ninja`); the regenerated graph adds the module to every
  binary's `SRC`/`OBJ` and changes nothing else (same 13 binaries, same
  flags, verified by diff against the pre-regeneration file).

## The instrumentation (what `QP_DVDSN_STATS=1` prints)

One `DAPPLY #<idx>` block per Cholesky apply, printed by
`davidson_apply_stats_mod` (a plain module in `src/utils/`, not an IRPF90
provider — the counters are written inside the OMP team, same decision as
`pair_block_stats.f90`).  Counters are exact integer sums (OMP ATOMIC,
thread-count invariant); the `t_*` fields are real*8 OMP ATOMIC sums of
per-thread brackets, i.e. **thread-seconds**, which equal wall-seconds at
one thread (diagnostic only, cf. `PAIR_BLOCK_CACHE.md`).

* Apply sections: `t_uscan` (u_is_sparse scan), `t_fill` ((k,l) block
  fill, the dgemm bracket), `t_ab` (alpha-beta sweep), `t_ss` (same-spin
  sweep), `t_sd` (singles+diag sweep), `t_transp` (the dtranspose pair
  around the work call).
* Fill census: `fill_pairs` (block dgemms actually run, hooked next to the
  `QP_PAIRBLOCK_STATS` call in `get_mo_two_e_integrals_ij_uncached` —
  reused, not duplicated) and `distinct` (distinct (k,l) pairs via a
  `mo_n*mo_n` stamp array marked under a named critical).
* Walk census: `conn_ab/ss/sd` (connections generated/kept per sweep),
  `reads_ab/rs` (slab integral lookups), `kcol_skip%` (the `kcol_prev`
  reuse, `u0_hs2_u0.irp.f:980`), `buf` GB (determinant-column bytes
  copied), `slab_written` GB (fill pairs × M² × 8 B).
* Beta-column occupancy: `cbar_a` = mean, `max_cola` = max connected-alpha
  list length over the same-spin gathers — the number that sizes the walk
  term `O(sze·c̄_a·N_int)`.
* Singles: per-apply delta of the process-global
  `get_single_excitation_from_fock` column-fill counter (2·mo_num doubles
  per call).
* `(i,j,k)` read census: min/max of the slab indices read per (thread,
  slab) cell, merged at apply end — the Phase-1A precondition check
  (are all reads inside the active window?).
* Driver attribution: `davidson_apply_stats_driver_enter/exit` in
  `davidson_diag_hjj_sjj` hold the apply line back until the driver's
  sections are timed, so `drv: t_sloop, t_dgemm, t_symm, t_upd1, t_upd2`
  (the `s_` loop, the h/s dgemms, the symmetrizations, the two update
  blocks) attach to the apply that produced the W/S columns.  Without a
  driver (final-energy apply, pt2 program) the line prints at
  `end_work`.
* Item 0.5 memory-model prints (always, no env gate needed — they are in
  the driver's existing startup print, report only): `H-apply slab
  8*M^3*8 (Gb)`, `big_array_coulomb/exchange (Gb)`, `ezfio window cache
  (Gb)` = `2^(4*shift)*8` B, next to the existing `r1` estimate.

**Build trap (hit, recorded):** a file-level `use davidson_apply_stats_mod`
in a `.irp.f` is **inert** — IRPF90 expands only in-scope `use` statements,
so the un-`use`d calls compiled as external procedures with plain
`_name_` symbols and the link failed with undefined references (the
in-scope `use omp_lib` convention at the top of the files is dead weight
too, as is `use map_module` at the top of `map_integrals.irp.f`: every
scope carries its own).  The `use` must sit inside each scoping unit that
calls the module.  Second trap: the call sites lack `iso_c_binding`, so
byte-count arguments are widened with `int(x,8)`.

## Phase-0 blocker — the CIPSI outer loop is infeasible on the mini

Decided 26 November 2026, 00:30, while `fixbz b_off` (attempt 1, log kept
as `fixbz-b_off.log.attempt1`) sat in its first diagonalization:

* The first diagonalization of the 114-MO fixture never converges:
  `threshold_davidson 1e-10` against the Cholesky-truncated Hamiltonian
  floors at a residual ~3.7e-5 (attempt 1: all 8 restarts × 15 applies
  burned, energy pinned at -231.20073004, residual 3.7–3.9e-5 in every
  restart). So **every** diagonalization burns the full ~120-apply
  budget — the restart loop in `diagonalization_hs2_dressed.irp.f` is
  hard-capped at `itertot == 8`, `itermax-1 = 15` applies per restart.
* The 46880-det space is **not** CIPSI-converged at `pt2_max 1e-4`: the
  prof campaign's stochastic PT2 on the same wavefunction gives an
  outside-space correction of -0.30538(6) Ha
  (`~/OpenCode/tmp/prof/full_s6_T.log`, 12 samples, 1 thread, 8401 s).
  With `selection_factor = 1` (default, `determinants/EZFIO.cfg`) the
  CIPSI loop (`stochastic_cipsi.irp.f:71-76`) would add up to ~46880
  determinants per outer iteration and iterate toward FCI: tens to
  hundreds of outer iterations, each a ~120-apply diagonalization plus a
  stochastic PT2 step = weeks on this 8-core/8-GiB box at load ~21–26.
  Not launchable — and the converged space is not the Phase-0
  measurement anyway.

**Documented cut (never a substitute physics run).** Each fixture arm is
terminated at a documented point with the in-tree clean-stop file
`work/qpkill` (`ezfio_files/qp_stop.irp.f`; polled per Davidson
iteration, per PT2-walk tooth, at the CIPSI-loop checkpoint): the run
exits rc=0 right after the cut, the log keeps everything up to the cut
plus the driver's final state report. Driver: `p0_driver2.sh` (cut
markers, per-arm 8 h backstop).

| arm | cut point | what it measures |
|---|---|---|
| b_off, b_on | after the first CIPSI summary (first `# E+rPT2` line) | iteration-1 E + rPT2, full 120-apply census (b_on), wall, peak RSS |
| a_v1, a_v8, a_t1 | end of the first diagonalization (PT2-step start) | the 120-apply census + wall under each BLAS threading (item 0.3); the partial PT2 summary they print is **not** used |
| t_on | self-terminating: `n_det_max 1428 = N_det`, the CIPSI loop never runs | single diagonalization of the fixed 1428-det space |

Fixture energies for the ULP-GATE are the **cut-point digits of this
campaign** (b_off vs b_on, same cut, must be bit-identical). They are
compared with the prof campaign's records (same representation,
`do_mo_cholesky T`) only within the documented floor:

| fixture | prof-campaign reference (`~/OpenCode/tmp/prof/`) |
|---|---|
| 114-MO | E denominator -231.2007184830427; E+rPT2 -231.50609618676967 ± 5.9e-4 (`full_s6_T.log`) |
| 264-MO | E denominator -230.9153852201646; E+rPT2 -231.58846737303898 ± 1.2e-4 (`tz_pt2_s6_t4.log`) |

The CIPSI diagonalization sits at the residual floor (~1e-5 from the
pt2-binary denominator) — that spread is the expected Cholesky-floor
level (residual ~3.7e-5 vs `threshold_davidson 1e-10`), not a drift.

### Campaign v2 outcome and the 14:00 correction (26 November)

Campaign v2 (driver2, launched 00:37, pre-fix binary) ran the
marker-based cuts, which proved unreachable. The logs established:

**Correction — the per-diagonalization-call apply budget is 105, not
15.** Each `diagonalize_CI` call runs 7 outer restarts × 15 applies
(`itertot` capped at 8 in `diagonalization_hs2_dressed.irp.f`), and
`diagonalize_ci.irp.f:105-155` then **doubles `N_states_diag`**
(4, 8, 16, 32, 64, 128, 256, …) and re-enters for every non-converged
call. So one CIPSI outer iteration ≈ 7–8 calls × 105 applies ≈
**~11 h at campaign load** (observed: b_on, 5.3 h → 673 applies).

* **b_off (driver2, stats off):** lost to the 8 h backstop at 08:37:28
  (rc=143). Log kept: 6 full calls (S=4…128) + ~10 applies of the
  S=256 call.
* **b_on (driver2, stats on):** its cut marker (`# E+rPT2`) sat behind
  a full CIPSI outer iteration (~11 h) — the 8 h backstop would have
  killed it anyway. The 13:55:44 touch of its `work/qpkill`
  **never fired**: `qp_stop()` is polled in the DAV drivers
  (`dav_general.irp.f:341-345`), the CIPSI checkpoints and the PT2
  walk, but **nowhere in `diagonalization_hs2_dressed.irp.f`**, so it
  is ignored inside a diagonalization call. SIGTERM'd at 14:04
  (rc=143) — **manually cut at 673 DAPPLY applies** (666 at the
  13:55:44 checkpoint), 105/105/105/105/105/105/43 per call
  S=4/8/16/32/64/128/256: the census is **complete over states
  4→256**.

**Per-apply split (114-MO, `VECLIB_MAXIMUM_THREADS=1`, 4 threads,
thread-seconds per apply, b_on; the orchestrator-verified table,
first-apply row verified line-for-line against `fixbz-b_on.log`):**

| apply | t_fill | t_ab | t_ss | t_sd | total |
|---|---|---|---|---|---|
| #1 (first, slab fill, 11664 fill pairs) | 20.06 | 9.62 | 6.97 | 0.34 | 37.0 |
| S=4 | 0.25 | 8.08 | 8.27 | 0.40 | 17.0 |
| S=8 | — | — | — | — | 23.4 |
| S=16 | — | — | — | — | 23.7 |
| S=32 | — | — | — | — | 25.5 |
| S=64 | — | — | — | — | 28.9 |
| S=128 | — | — | — | — | 50.5 |
| S=256 | 1.04 | 26.01 | 42.60 | 1.17 | 70.8 |

Per-apply values scatter within a state size (S=256 applies span
~17 to ~71; the last logged apply #674: 2.45/28.25/21.46/0.82).
**walk (ab+ss) = 96–97 % of every apply; fill = 1.2–1.5 % after the
first apply** (the pair cache reuses blocks across applies within the
diagonalization). Census per apply (representative, S=4): conn_ab
13.49M generated/kept, conn_ss 90.79M generated / 6.05M kept, buf
35.2 GB, kcol_skip 80.5 %, cbar_a 98.2, max_cola 7395, active window
i_,j_,k_ ∈ [7:114] (frozen core 1–6, 108 active MOs).

**Decision consequence (item 0.6 rule):** at 114-MO **fill « 50 %** →
Phase 5 (GPU) is **not armed**; the 264-MO arm (t_on) decides the
large-M regime. The 114-MO lever is **Phase 2 (walk merge)**.

**Tooling fix (rebuild 14:09, `build_instr5.log`, ninja rc=0):** a
`qp_stop()` poll added at the top of the inner iteration loop of
`diagonalization_hs2_dressed.irp.f`, modelled on
`dav_general.irp.f:340-344` (guarded exit with `converged = .True.`).
One INQUIRE per iteration without `work/qpkill` — no numerical effect.
The same rebuild picked up the staged `slab_idx` PRIVATE fix
(`u0_hs2_u0.irp.f:936`).

### Final-binary ULP-GATE anchor (14:12, tags fp1/fp4/fp8, fi1s/fi4s/fi8s)

| threads | E (stats off) | E (stats on) | N_det | rc |
|---|---|---|---|---|
| 1 | -7.882324392913227 | -7.882324392913227 | 69 | 0 / 0 |
| 4 | -7.882324392913227 | -7.882324392913227 | 69 | 0 / 0 |
| 8 | -7.882324392913227 | -7.882324392913227 | 69 | 0 / 0 |

Bit-identical to the recorded reference (all three rounds: pre-fix
instrumented, this final binary, and on/off pairs). No `DAPPLY` leak
(map-mode anchor). **PASS.**

## Phase 0 baseline

_(tables filled from the campaign logs in
`~/OpenCode/tmp/qp3-hpsi-p0/`: `anchor-*.log`, `fixbz-*.log`,
`fixtz-*.log`, driver log `p0_driver2.log`)_

### Anchor — 69-det strict FCI (LiH/sto-3g, `pt2_max 1e-8`, `n_det_max 300`)

Reference (`BASELINE.md`, post-removal tree): N_det = 69,
E = **-7.882324392913227**.

| threads | E (stats off) | E (stats on) | N_det | wall (s)¹ | rc |
|---|---|---|---|---|---|
| 1 | -7.882324392913227 | -7.882324392913227 | 69 | 6 / 4 | 0 / 0 |
| 4 | -7.882324392913227 | -7.882324392913227 | 69 | 5 / 5 | 0 / 0 |
| 8 | -7.882324392913227 | -7.882324392913227 | 69 | 5 / 4 | 0 / 0 |
| 28 | -7.882324392913227 | -7.882324392913227 | 69 | 5 / 5 | 0 / 0 |

¹ driver-level wall (seconds, launch overhead included; the 28-thread
run oversubscribes the 8-core box — wall is noise, bits are the gate).

**ULP-GATE items 2+5: PASS.** Bit-identical across 1/4/8/28 threads and
stats on/off (14th digit fixed — same saved determinant order within the
tree). No `DAPPLY` lines at the anchor: it is a map-mode run, and the
census is instrumented in the Cholesky apply body only (the plan's
target mode); the driver-section timers still ran and stayed inert.

### 114-MO Cholesky CIPSI (`benzene_chol.ezfio`, 46880 dets, shift 6)

Cut per the blocker section above (first CIPSI summary).
Run A (baseline, stats off): OMP_NUM_THREADS=4, QP_NTHREADS_PT2=4,
VECLIB_MAXIMUM_THREADS=1, QP_MAXMEM=5.
Run B (same, stats on): must be bit-identical to A.

| run | CI E (cut) | E+rPT2 (iter 1) | N_det | WALL (s) | peak RES MEM | last `DAPPLY` census |
|---|---|---|---|---|---|---|
| b_off | | | | | | — (stats off) |
| b_on | | | | | | `fixbz-b_on.log` |

### 264-MO Cholesky CIPSI (`tz_basis_test.ezfio`, 1428 dets, shift 6)

Run C: OMP_NUM_THREADS=4, QP_NTHREADS_PT2=4, VECLIB_MAXIMUM_THREADS=1,
QP_MAXMEM=4, stats on.  Fixed 1428-det space (`n_det_max 1428`), so a
pure Davidson run: no selection, no pruning — the fill-dominated regime.
`io_mo_cholesky = Read` (3.1 GB `cholesky_mo_transp`, mmap).
Self-terminating (no cut needed).

| run | CI E | E+rPT2 (prof ref) | N_det | WALL (s) | peak RES MEM | last `DAPPLY` census |
|---|---|---|---|---|---|---|
| t_on | -230.9153852202491 | -231.58846737303898 ± 1.2e-4 | 1428 | 46 | 2.976 GB | `fixtz-t_on.log` (4 applies, S=4) |

Converged in **2 iterations** (residual 7.0e-11 → 5.2e-12 < 1e-10):
the stored 1428-det wavefunction is already an eigenstate of the
truncated Hamiltonian to below `threshold_davidson`, so no
state-doubling retry (contrast the 114-MO fixture, whose residual
floors at ~3.7e-5). CI E matches the prof-campaign denominator
(-230.9153852201646) to Δ = 8.5e-11 — ULP level (different binary/
BLAS state in the prof campaign), same representation.
**Large-M regime (item 0.6):** first apply fill = 97.8 thread-s
(3364 fill pairs, slab 1.75 GB written); applies ≥2: fill 0.33–1.64
(pair cache: 1 pair, then 0), walk ab+ss ≈ 0.18 → **fill 75–95 % of
the apply** — the opposite of the 114-MO walk-dominated regime.
cbar_a 10.4 (vs 98.2 at 114-MO), active window i_,j_,k_ ∈ [7:91].

Fixture energies are **this** run's digits: later phases must match within
representation (Cholesky mode only — the ~3.8e-6 relative shift to map
mode is real, never compare across `do_mo_cholesky`).

## Item 0.3 — BLAS-threading A/B (no code change)

**The plan's lever does not exist on this machine.**  The `c_funcs.c` pin
(`openblas_set_num_threads(1)` + `OPENBLAS_NUM_THREADS=1`, gated on
`/proc/self/comm == "pt2"`) is inside `#ifdef __linux__` and never
compiles on macOS; the BLAS here is Accelerate/vecLib, whose thread pool
is set by `VECLIB_MAXIMUM_THREADS` (the PT2 walk recipes already pin it
to 1).  So the mini A/B substitutes the platform lever, on the 114-MO
fixture, OMP_NUM_THREADS=4:

Cut at the end of the first diagonalization (blocker section); all arms
with `QP_DVDSN_STATS=1` so the fill/walk split is measured under each
threading.

| arm | env | WALL (s)¹ | CI E (cut)² | vs pinned bits | DAPPLY split |
|---|---|---|---|---|---|
| a_v1 (pinned) | VECLIB_MAXIMUM_THREADS=1 | | | — | `fixbz-a_v1.log` |
| a_v8 (free) | VECLIB_MAXIMUM_THREADS=8 | | | | `fixbz-a_v8.log` |
| a_t1 (serial) | OMP_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1 | | | | `fixbz-a_t1.log` |

¹ first-diagonalization wall (120 applies).
² the partial PT2 summary printed after the cut is not used (see the
blocker table).

(The plan's arm (b) "OMP_NUM_THREADS=1 with nthreads_davidson at the
target" is not launchable as written: `nthreads_davidson` is computed
from `nproc` (the OMP default, `diagonalization_hs2_dressed.irp.f:219,263`)
and there is no keyword to override it — it is not a
`QP_NTHREADS_PT2`-style knob.  Serial arm (a_t1) is the closest
realisable analogue and is the bit reference for the thread sweep.)

**Verdict → item 0.4:** _(filled: wall gain or thread-count bit drift
extends the gate to the CIPSI executable; on macOS that means either a
platform-independent `QP3_PIN_BLAS`-style env gate in `c_funcs.c` or the
launch contract `VECLIB_MAXIMUM_THREADS=1` — the current gate cannot be
"extended" because it does not compile on this OS at all.)_

## Item 0.5 — memory-model audit (report only)

The driver now prints, next to its `r1` estimate (U/W/S/u_t/v_t/s_t and
the per-thread buffers):

| term | formula | 114-MO | 264-MO |
|---|---|---|---|
| H-apply slab | 8·M³·norb_block, norb_block=8 | 0.0883 GB | 1.0967 GB |
| big_array_coulomb/exchange | 2·M³·8 B | 0.0221 GB | 0.2742 GB |
| ezfio window cache | 2^(4·shift)·8 B | 0.125 GB (shift 6) | 0.125 GB (shift 6) |

Cross-checked against observed peak RSS: 114-MO **2.408 GB** (b_on),
264-MO **2.976 GB** (t_on). The printed terms cover the H-apply
static part (0.235 / 1.496 GB); the remainder is pair-cache
residency (caps 2.56 / 2.05 GB, partially touched) plus the
wavefunction and integral arrays.

## Item 0.6 — production-size sizing

Measured fill/walk split (thread-seconds per apply, post-first-apply):

| regime | fill | walk (ab+ss) | dominant |
|---|---|---|---|
| 114-MO / 46880 dets (b_on) | 0.25–1.04 | 16.4–68.7 | **walk (96–97 %)** |
| 264-MO / 1428 dets (t_on) | 0.33–1.64 | ≈0.18 | **fill (75–95 %)** |

First-apply fill (slab build): 20.1 (114-MO) / 97.8 (264-MO) thread-s.

**Production sze ≥ 1e6 (extrapolated, flagged):** walk scales
O(sze·c̄_a) — at 1e6 dets the walk term is ≈126 (264-MO, c̄_a 10.4)
to ≈360 (114-MO, c̄_a 98.2) thread-s against fill ≈0.3–1.6, i.e.
**both M values flip back to walk-dominated at production sze**; the
fill-dominated regime is a small-sze (and first-apply) effect.
Decision rules: 114-MO fixture → Phase 2 first, Phase 5 **not**
armed; 264-MO fixture (small sze) → Phase 1 first, Phase 5.1 armed;
at the production sze the walk merge (Phase 2) is again the common
lever. **Mini constraint: 8 GiB — the production-size run itself
must be specified for the dell; this sizing is extrapolated from the
two fixture runs and flagged as such.**

## ULP-GATE (Phase 0 verdict)

1. Rebuild: `build_instr5.log` (14:09, ninja rc=0) — final binary:
   instrumentation + `qp_stop()` poll (diagonalization) + `slab_idx`
   PRIVATE fix.
2. Anchor bits: pre-fix binary 1/4/8/**28** (00:08, table above);
   **final binary 1/4/8 (14:12, table above) — PASS.**
3. vs the pristine Phase-0.1 record: plain PASS so far — no numerical
   change (counters and clock reads only); the 14:09 change is a polled
   file test inert without `work/qpkill` (anchor re-run bit-identical).
4. Fixture digits: _(pending driver3: b2_off/b_on stats on/off A/B;
   a_v1/a_v8/a_t1 cut at call 1; t_on self-terminating)_
5. Stats on/off bit-identical: _(anchor PASS both rounds; fixture
   pending driver3)_
6. Phase record appended: this file.

**Verdict: _(PASS / PASS-ULP / FAIL — pending driver3 fixture arms)_**

## Decision rules (from the Phase-0 numbers)

* fill < 50 % of apply and walk > 50 % → Phase 2 first; Phase 5 stays off.
* fill ≥ 50 % → Phase 1 first; Phase 5.1 armed.
* _(filled with the measured `t_fill` vs `t_ab+t_ss+t_sd` split at both
  fixture sizes — note the split is fixture-regime dependent: at 264-MO/
  1428-dets the fill is pair-count bound, at 114-MO/46880-dets the walk
  is determinant-count bound.)_
