# Phase 0 — baseline + instrumentation of the Cholesky H|Ψ⟩ apply (Davidson/CIPSI)

Phase 0 (R0) of `~/OpenCode/tmp/qp3-hpsi-PLAN.md`, the merged dual-brain
report at `~/OpenCode/tmp/qp3-davidson-hpsi-report.md`.  Measurement phase:
no numerical change.  All instrumentation is env-gated
(`QP_DVDSN_STATS=1`); a run with the counters on must be bit-identical to
the same run with them off (ULP-GATE item 5).

## Where this stands (Phase 0 CLOSED, 26 November 2026 16:45)

Instrumentation built and bit-gated (anchor ULP-GATE **PASS**:
E = -7.882324392913227, 69 dets, 1/4/8/28 threads, stats on/off,
both binaries). Both fixture regimes measured: 114-MO is
walk-dominated (96–97 % of the apply) → Phase 2 (merge slab walks)
is the lever, Phase 5 (GPU) not armed; 264-MO is fill-dominated
(75–95 %) → Phase 1 first, Phase 5.1 armed. Production-size (sze ≥
1e6) sizing is extrapolated from the two fixtures and the production
run itself is deferred to the dell (31 GiB); the BLAS-pin A/B (item
0.4) is deferred to the Phase-1 spec. Final numbers and verdict:
"Campaign v3 outcome (driver3)" below.

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

_(The table records the original design. In practice the b_off/b_on
marker cuts proved unreachable (Campaign v2: the per-call budget is
105 applies and the first `# E+rPT2` sits behind a full ~11 h outer
iteration) and the a_ arms ran as driver3 count cuts at 105
`DAPPLY` blocks; the first-iteration E+rPT2 the clean cuts print is
used and agrees with the prof record. See "Campaign v2 outcome" and
"Campaign v3 outcome (driver3)" for what actually ran.)_

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
_(The ~11 h is the estimate for a chain that doubles all the way to
S=512; the observed b_on chain stopped at S=256 after ~5.4 h — see
the b_on tail in Campaign v3, which refines the per-outer-iteration
cost to ~5.5–6 h. The infeasibility verdict is unchanged and if
anything strengthened.)_

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
  (rc=143) — **manually cut at 674 file DAPPLY blocks** (666 at the
  13:55:44 checkpoint) = 1 initial guess + 672 in-call
  (105/105/105/105/105/105/42 per call S=4/8/16/32/64/128/256) +
  1 PT2-step apply (the last block, printed after the CI state
  report — Campaign v3, b_on tail): the census is **complete over
  states 4→256**.

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
~17 to ~71; note the log's final block #674 — 2.45/28.25/21.46/0.82
— is the post-cut PT2-step apply, not a chain apply: Campaign v3,
b_on tail).
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

### Campaign v3 outcome (driver3, final binary)

Launched 14:15:39 on the final binary (14:09 rebuild) with
count-based cuts (`p0_driver3.sh`: cut when
`grep -cE pattern log` reaches N, then touch `work/qpkill` — now
honoured mid-diagonalization by the Step-2 poll), 3 h per-arm
backstop (t_on 6 h), `rm -f` of `qpkill` at each arm start, fresh
ezfio dirs for the b2_off/a_* arms (old b_off/b_on dirs untouched).
Machine load ~12 at launch, spiking to 23–32 between 14:15 and
15:30. **All five arms rc=0, finished 16:44:21; every count cut
fired and was honoured:**

| arm | cut | cut fired | rc at |
|---|---|---|---|
| t_on | none (self-terminating) | — | 14:16:25 (46 s) |
| b2_off | 105 iteration lines (reached 125) | 14:36:02 | 14:40:48 |
| a_v1 | 105 `DAPPLY` blocks | 15:00:56 | 15:06:27 |
| a_v8 | 105 `DAPPLY` blocks | 15:25:06 | 15:31:38 |
| a_t1 | 105 `DAPPLY` blocks | 16:19:16 | 16:44:21 |

The b2_off marker counts iteration lines (call 1, S=4, yields 83 —
the cut lands ~22 lines into call 2, S=8); the a_* markers count
`DAPPLY` blocks, so their cut lands at the end of call 1 plus 0–2
cut-tail applies. b2_off carries no `DAPPLY` at all (stats off).

**Per-arm S=4 split, thread-seconds per apply** (first 105 DAPPLY
blocks, excluding the first slab-fill apply; the `t_*` sections are
wall-based thread-seconds, so cross-window comparisons are
load-sensitive — see the WALL/CPU table):

| arm | VECLIB | OMP thr | stats | fill | ab | ss | sd | total |
|---|---|---|---|---|---|---|---|---|
| b_on (ref, load ~12) | 1 | 4 | on | 0.25 | 8.08 | 8.27 | 0.40 | 17.0 |
| a_v1 | 1 | 4 | on | 0.37 | 23.94 | 19.41 | 0.97 | 44.69 |
| a_v8 | 8 | 4 | on | 0.36 | 22.62 | 17.90 | 0.86 | 41.73 |
| a_t1 | 1 | 1 | on | 0.12 | 11.12 | 15.37 | 0.75 | 27.36 |

Spot-checked against the logs (recomputed from the `DAPPLY`
blocks): totals 44.74 / 42.20 / 27.12 vs the table's 44.69 / 41.73
/ 27.36 — all within 1.2 %. (The b_on reference row is the driver2
load-~12 window; the driver3 arms ran under the load spike, hence
the ~2.5–3× inflation of the same kernel.)

**Whole-run WALL/CPU** (setup + cut tail + the post-cut PT2 step
included — every clean cut ends fake-converged, so CIPSI runs the
first-iteration PT2 step, prints its summary, and exits rc=0;
b2_off: cut 14:36:02 → rc=0 14:40:48, i.e. that step takes ~5 min
at 4 threads here):

| arm | WALL (s) | CPU (s) | applies | CPU-s/apply |
|---|---|---|---|---|
| b2_off | 1461.7 | 4761.5 | ~140 (cut count 125) | ≈34 |
| a_v1 | 1536.5 | 4653.8 | 105 | 44.3 |
| a_v8 | 1500.3 | 4648.4 | 106 | 43.9 |
| a_t1 | 4352.2 | 2794.2 | 106 | 26.4 |

Spot-checked: b2_off final line 1461.73 / 4761.49, a_v1
1536.54 / 4653.82, a_v8 1500.29 / 4648.35, a_t1 4352.15 / 2794.23.
"applies" = in-call (S=4) DAPPLY blocks; each driver3 log also
carries one pre-header block (the initial-guess energy apply after
pair-cache init) plus 0–2 cut-tail blocks (a_v1: 1 + 105 + 2 = 108
file blocks).

**Findings.**

1. **The H-apply Davidson start vector is random (unseeded).** The
   extra diagonalization states are built as Box–Muller Gaussian
   twists of the previous state plus a +10 diagonal kick
   (`diagonalization_hs2_dressed.irp.f:360-372`), from
   `random_number` which is seeded nowhere on the non-stochastic
   path (`random_seed` appears only in
   `pt2_stoch_routines.irp.f`). Verified directly: two runs of the
   same freshly compiled test program give different first
   `random_number` values on this machine (gfortran's unseeded
   default seed is not run-reproducible here). Consequence:
   **per-iteration energy comparison between runs at fixture level
   is meaningless.** The comparable pair (call 1, iteration 1,
   stats off): -231.2007297346 (b_off, log line 94) vs
   -231.2007303954 (b2_off, log line 94) — 6.6e-9 apart. (The value
   -231.2007302790 cited in the orchestrator's notes is b_on's
   call-5 (S=64) iteration-10 line, log line 6239, not a b_off
   run's iteration-1 energy; the conclusion is unchanged. The two
   runs also differ in binary — pre/post 14:09 rebuild, an inert
   `qp_stop` poll + the `slab_idx` PRIVATE fix — which does not
   affect the mechanism.) Only the anchor (and converged energies)
   are bit-comparable across runs; the bit contract is adjudicated
   by the anchor (ULP-GATE PASS, E = -7.882324392913227).
2. **`VECLIB_MAXIMUM_THREADS` 1→8: CPU-neutral** (4653.8 vs
   4648.4, 0.1 %), ~7 % in wall-thread time (44.69 vs 41.73) over
   adjacent windows. Not a lever for this kernel; the OMP team size
   and the BLAS pool are the threading levers.
3. **Parallel speedup + spin under load.** 4-thread H-apply wall ≈
   4.25 s/apply at load ~12 (b_on S=4: 105 applies in ≈443 s) vs
   1-thread a_t1 ≈ 27–41 s/apply (load-confounded windows) → ≥2×
   speedup. Under the load spike (load 23–32, 14:15–15:30),
   4 threads consumed MORE total CPU than 1 thread (a_v1 44.3 vs
   a_t1 26.4 CPU-s/apply at S=4) — spin-waste from the unpinned BLAS
   pool (the `c_funcs.c:59` pin is pt2-only and `#ifdef
   __linux__`). This is the operational signature of item 0.4's
   contract hole at fixture scale; the clean pin A/B measurement is
   deferred to the Phase-1 spec (the fill dgemm is the object of
   test there; ~1 arm, ~30 min, on a quiet machine or the dell).
4. **Instrumentation cost under the load spike is unquantified**
   (b2_off 29.2 vs a_v1 44.3 CPU-s/apply at S=4, +52 %, but the
   load-spoke timing confounds the two adjacent windows).
   Bit-inertness is proven at the anchor (on/off bit-identical).
   Design note for future high-load rate campaigns: batch the
   per-connection counters per thread and reduce once per apply
   (104 M contended atomic adds per apply at this fixture), or run
   rate A/B arms with stats off.
5. **The qpkill cut mechanism now works mid-diagonalization** (all
   four driver3 arms cut cleanly at their count markers, rc=0).
   Campaign tooling is fixed for future use.

**b_on tail (driver2, for the record):** its first CIPSI outer
iteration completed the diagonalization chain — 7 calls S=4→256,
672 in-call applies over **~5.4 h** (final banner 19480 s; the
chain started right after a ~16 s setup) — and the S=256 call
returned after 3 restarts (42 applies) with the printed residual
still at ~5.7e-5, far from the 1e-10 threshold (the driver carries
a LAPACK `dsygv info>0` subspace-failure path that shortens
restarts, `diagonalization_hs2_dressed.irp.f:525-531`, consistent
with this exit; no S=512 call appears in the log). The CIPSI loop
then printed the CI state energy **-231.2007298588983** (S² 8.99e-5,
A1g purity 0.9961 — the `CI_energy` provider,
`diagonalize_ci.irp.f:27`) and **started the first PT2 step**
("PT2 Energy denominator -231.2007298588983") before the 14:04
SIGTERM (rc=143). Note: b_on's `DAPPLY #674` is that PT2-step
apply, not a diagonalization apply — it prints after the CI state
report (b_on's 674 file blocks = 1 guess + 672 in-call + 1 PT2).
One full CIPSI outer iteration is thus ≈ 5.4 h chain + ~5 min PT2
step + selection ≈ **~5.5–6 h observed** — shorter than the ~11 h
estimate (which assumed the doubling chain runs to S=512; b_on's
S=256 call returned early via the dsygv path, one call short). The
blocker verdict is unchanged: tens to hundreds of such outer
iterations at `pt2_max 1e-4` is days-to-weeks on the mini either
way, and the per-iteration cost is lower, not higher, than
estimated.

**Verdict (Phase 0 closed).**

* 114-MO (b_on, 672 applies, load ~12): **walk (ab+ss) = 96–97 % of
  the apply, fill = 1.2–1.5 %** (post-first-apply; the pair cache
  reuses blocks across applies within the diagonalization). Phase 5
  (GPU) not armed at 114-MO; the lever is Phase 2 (merge slab
  walks).
* 264-MO (t_on, 46 s, converged in 2 iterations, E =
  -230.9153852202491): **fill = 75–95 %** (first-apply fill 97.8
  thread-s, 3364 (k,l) pairs). Phase 5 armed → the Turpan A/B
  target is the one-time 3364-pair fill (host 97.8 thread-s @4
  threads vs A100 + transfer, at the post-Phase-1 compacted shape:
  n_act=43 window, working set ≈ 50 MB).
* Production sze ≥ 1e6: walk flips back to dominant at both M
  (extrapolated, flagged) → both levers; order Phase 1 (compact +
  persist) first, then Phase 2.
* Open items: item 0.4 (BLAS pin A/B) deferred to Phase 1 as above;
  item 0.6 (production-sze fixture) deferred to the dell (31 GiB) —
  the mini (8 GiB) cannot host sze ≥ 1e6.

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

Cuts as actually run (marker cuts of the blocker table proved
unreachable — Campaign v2; count cuts below — Campaign v3). All
arms: OMP_NUM_THREADS=4, QP_NTHREADS_PT2=4,
VECLIB_MAXIMUM_THREADS=1 (script baseline), QP_MAXMEM=5.
CI E (cut) is the fake-converged state report (a mid-Davidson
snapshot); E+rPT2 (iter 1) is the first-outer-iteration PT2
estimate the clean cut lets CIPSI print. The stats on/off
bit-identity is adjudicated at the anchor, not here (Campaign v3,
finding 1: unseeded random start vector — fixture runs are not
bit-comparable across runs).

| run | cut | CI E (cut) | E+rPT2 (iter 1) | N_det | WALL (s) | peak RES MEM | `DAPPLY` |
|---|---|---|---|---|---|---|---|
| b_off (driver2, off) | 8 h backstop, rc=143 (mid S=128 call) | — (no state report; last iter line -231.2007312413) | — | 46880 | 28819 | 2.363 GB | — (off) |
| b_on (driver2, on) | manual SIGTERM 14:04, rc=143 (mid first PT2 step) | -231.2007298588983 | — (PT2 step in progress) | 46880 | 19496 (5.4 h) | 3.450 GB | 672 in-call + 1 guess + 1 PT2 (`fixbz-b_on.log`) |
| b2_off (driver3, off) | count 125, rc=0 | -231.2007309231487 | **-231.50588378561073 ± 9.0e-4** | 46880 | 1461.7 | 2.251 GB | — (off) |

prof reference (same representation): E+rPT2 -231.50609618676967 ±
5.9e-4 → b2_off agrees to Δ = 2.1e-4 (inside both error bars). The
a_v1/a_v8/a_t1 arms (item 0.3 table below) carry the same
first-iteration E+rPT2: -231.50621446873146 ± 1.5e-3 (a_v1),
-231.50684889450741 ± 1.7e-3 (a_t1); a_v8's cut snapshot put a
S² = 0.57 state in slot 1 (E -229.2446, purity 0.67) → excluded.

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

Cut at the end of call 1 (105 `DAPPLY` blocks, driver3 count cuts);
all arms with `QP_DVDSN_STATS=1`, OMP_NUM_THREADS as listed,
VECLIB_MAXIMUM_THREADS as listed (script baseline = 1), QP_MAXMEM=5.
"applies" = in-call (S=4) DAPPLY blocks; the split row is
post-first-apply thread-seconds (Campaign v3 load caveat).

| arm | OMP | VECLIB | WALL (s) | CPU (s) | applies | fill/ab/ss/sd (thread-s/apply) | cut CI E | E+rPT2 (iter 1) |
|---|---|---|---|---|---|---|---|---|
| a_v1 (pinned) | 4 | 1 | 1536.5 | 4653.8 | 105 | 0.37/23.94/19.41/0.97 | -231.2007229681 | -231.50621446873146 ± 1.5e-3 |
| a_v8 (free) | 4 | 8 | 1500.3 | 4648.4 | 106 | 0.36/22.62/17.90/0.86 | -229.2445867898738 (S² 0.57 artifact) | — (excluded) |
| a_t1 (serial) | 1 | 1 | 4352.2 | 2794.2 | 106 | 0.12/11.12/15.37/0.75 | -231.2001690451476 | -231.50684889450741 ± 1.7e-3 |

The original "vs pinned bits" column is void: the bit contract is
adjudicated at the anchor (Campaign v3, finding 1) — fixture runs
are not bit-comparable across runs (unseeded random start vector).
The cut-point CI E values are mid-Davidson snapshots (a_t1's is
5.6e-4 above the floor — cut within call 1, S=4 subspace only); the
usable digits are the first-iteration E+rPT2 (all within ~1e-3 of
the prof record) and the census.

(The plan's arm (b) "OMP_NUM_THREADS=1 with nthreads_davidson at the
target" is not launchable as written: `nthreads_davidson` is computed
from `nproc` (the OMP default, `diagonalization_hs2_dressed.irp.f:219,263`)
and there is no keyword to override it — it is not a
`QP_NTHREADS_PT2`-style knob.  Serial arm (a_t1) is the closest
realisable analogue; the bit reference for any thread sweep is the
anchor, not a fixture run — finding 1 of Campaign v3.)

**Verdict → item 0.4:** `VECLIB_MAXIMUM_THREADS` 1→8 is CPU-neutral
(0.1 %) and ~7 % in wall-thread time over adjacent windows (Campaign
v3, finding 2) — not a lever for this kernel; the OMP team size and
the BLAS pool are the threading levers. The `c_funcs.c` pin is
pt2-only and `#ifdef __linux__` (never compiles here), so no
platform-independent `QP3_PIN_BLAS`-style gate exists on this OS;
the operational signature of the unpinned-pool spin waste was
measured at fixture scale (under the load spike, 4 threads cost more
CPU-s/apply than 1 thread — finding 3). The clean pin A/B (object:
the fill dgemm) is **deferred to the Phase-1 spec**: ~1 arm, ~30
min, on a quiet machine or the dell.

## Item 0.5 — memory-model audit (report only)

The driver now prints, next to its `r1` estimate (U/W/S/u_t/v_t/s_t and
the per-thread buffers):

| term | formula | 114-MO | 264-MO |
|---|---|---|---|
| H-apply slab | 8·M³·norb_block, norb_block=8 | 0.0883 GB | 1.0967 GB |
| big_array_coulomb/exchange | 2·M³·8 B | 0.0221 GB | 0.2742 GB |
| ezfio window cache | 2^(4·shift)·8 B | 0.125 GB (shift 6) | 0.125 GB (shift 6) |

Cross-checked against observed peak RSS (max printed banner):
114-MO **3.450 GB** (b_on — the S=256 call's U/W/S workspace sits
~1 GB above the S=128 peak of 2.363 GB in b_off, which died before
the S=256 call; b_on's tail value of 2.408 GB is at the kill),
264-MO **2.976 GB** (t_on). The printed terms cover the H-apply
static part (0.235 / 1.496 GB); the remainder is pair-cache
residency (caps 2.56 / 2.05 GB, partially touched), the
Davidson workspace (scales with N_st_diag·itermax) and the
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
4. Fixture digits: t_on CI E -230.9153852202491 matches the prof
   denominator (-230.9153852201646) to Δ = 8.5e-11 (ULP level,
   same representation). 114-MO first-iteration E+rPT2: b2_off
   -231.50588378561073 ± 9.0e-4, a_v1 -231.50621446873146 ± 1.5e-3,
   a_t1 -231.50684889450741 ± 1.7e-3 vs the prof record
   -231.50609618676967 ± 5.9e-4 — all within the error bars (Δ
   2.1e-4 / 1.2e-4 / 7.5e-4); a_v8 excluded (cut-snapshot state
   ordering, S² 0.57 in slot 1). Per-iteration energies are **not**
   bit-comparable across runs — unseeded random start vector
   (Campaign v3, finding 1); the bit contract is adjudicated by the
   anchor.
5. Stats on/off bit-identical: anchor PASS (both binaries, all
   thread counts). At fixture level an on/off bit comparison is
   untestable in principle for the same reason (finding 1); the
   fixture-level evidence is that on/off runs track each other at
   the run-to-run scatter floor (call-1 iter-1 within ~3e-9;
   first-iteration E+rPT2 within ~1.5e-3).
6. Phase record appended: this file.

**Verdict: PASS — Phase 0 closed (26 November 2026).**

## Decision rules (from the Phase-0 numbers)

* fill < 50 % of apply and walk > 50 % → Phase 2 first; Phase 5 stays off.
* fill ≥ 50 % → Phase 1 first; Phase 5.1 armed.
* Filled (this campaign): 114-MO/46880-dets — walk (ab+ss) 96–97 %
  of post-first-apply → Phase 2 first, Phase 5 stays off;
  264-MO/1428-dets — fill 75–95 % → Phase 1 first, Phase 5.1 armed.
  The split is fixture-regime dependent: at 264-MO/1428-dets the
  fill is pair-count bound (3364 (k,l) pairs, c̄_a 10.4), at
  114-MO/46880-dets the walk is determinant-count bound (c̄_a 98.2).
  At production sze ≥ 1e6 the walk re-dominates at both M
  (extrapolated) → Phase 2 is the common lever; order Phase 1
  (compact + persist) first, then Phase 2.
