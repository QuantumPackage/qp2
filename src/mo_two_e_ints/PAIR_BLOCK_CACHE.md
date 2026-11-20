# Static hot-block cache for the Cholesky MO-integral block (PT2 selection)

Spec for implementation. Stage 0 is a measurement with no behaviour change and
decides whether Stage 1 is written at all.

## The cost being attacked

The PT2 selection reaches `mo_two_e_ints` through two Cholesky entry points, and
in production (`mo_two_e_ints/do_mo_cholesky = True`, `mo_num` >
2^`mo_integrals_cache_shift`) both contract the stored vectors:

* `get_mo_two_e_integrals_ij(k,l,...)` computes the whole `mo_num x mo_num`
  block `<ij|kl>` for all `i,j`: one dgemm, `2*mo_num^2*cholesky_mo_num` flops
  (`map_integrals.irp.f:666-739`).
* `get_mo_two_e_integrals(j,k,l,...)` computes **one column** of that same
  block, `<ij|kl>` for all `i`: one dgemv over the MOs that the window does not
  cover, on each side of it (`map_integrals.irp.f:370-389`, `:539-571`).

Both return `sum_P transp(P,i,k)*transp(P,j,l)` for the same `(k,l)`: column `j`
of block `(k,l)` **is** the vector the second routine returns. The two entry
points share one key, and one cached block serves both.

Measured on the cc-pVDZ probe below (114 MOs, 2694 vectors, window `(1,64)`, one
thread, 20 s of `sample` at 1 ms, 10622 samples inside the walk):

```
  DGEMV in get_mo_two_e_integrals        8209   77.3 % of the walk
  DGEMM in get_mo_two_e_integrals_ij     1390   13.1 %
  everything else                        1023    9.6 %
```

So the column stream is the wall time: five to six times the block GEMM. In the
probe a single dgemv is `cholesky_mo_num x (mo_num - cache_max)` = 2694 x 50
doubles: it moves 1.08 MB to deliver 50 doubles, 400 B - an arithmetic intensity
of 0.25 flop/byte, i.e. pure bandwidth. That is why the lever is a cache and not
a better kernel: the same numbers served from RAM are a 400 B copy, and at 300
MOs with ~1500 vectors and a level-7 window (1-128) the dgemv moves 2.1 MB to
deliver 172 doubles. The bytes grow as `nvec*(mo_num - 2^shift-1)` while the
useful payload grows with `mo_num` alone, so the ratio gets worse, not better,
with system size.

Where the two come from, which is why they are counted separately:

* `mo_integrals_cache` covers only a window of 2^shift orbitals centred on the
  Fermi level (`map_integrals.irp.f:55-62`; level 7 = 128 orbitals, 2 GB), so
  both routines take a Cholesky branch (`map_integrals.irp.f:662-695`). The
  sparse MO map is **not built** in this mode
  (`mo_bi_integrals.irp.f:65-82`), so every scalar integral outside the window
  is a ddot over the Cholesky vectors (`map_integrals.irp.f:250-292`) and every
  vector call is a dgemv (`:370-389`, `:539-571`).
* the block is called by `get_d0` alone, once per invocation: two static sites,
  one per branch of the `sp == 3` test (`cipsi/selection.irp.f:1551` for AB,
  `:1576` for AA/BB), and `get_d0` is entered only for `nt == 2` selectors
  (`selection.irp.f:1035-1036`). The column is called by `get_d1` (nt == 3),
  eight times per invocation (`:1305-1440`); `get_d2` (nt == 4) calls neither.
  The two classes therefore come from different selector subsets, and the ~6:1
  ratio is a property of the determinant set: it must be re-measured at
  production size, not extrapolated from `mo_num`.
* the block GEMM itself is `2 * mo_num^2 * cholesky_mo_num` flops. At 300 MOs
  with 1000-2000 vectors that is 0.18-0.36 GFlop: 4-9 ms on one core at
  40-60 GFlop/s. At ~45 flop/byte (4-8 MB of operands for 720 KB of output) it
  is the one genuinely compute-bound piece of the walk; the column is its
  opposite at 0.16 flop/byte. The `mo_num^2` scatter the block feeds costs
  ~0.8 ms at the same size (`tau ~ 8.5 ns` per element), so the GEMM is ~85 %
  of an `nt == 2` call.
* The key `(k,l)` is the **virtual pair of the excitation mask -> det**: `p` is
  the list of the det's orbitals that are not in the mask
  (`selection.irp.f:1017-1018`), and it is passed reversed - `(h2,h1)` for AB,
  `(p2,p1)` for AA/BB. `bitstring_to_list` yields ascending orbitals, so a
  same-spin key always has `k > l` while an AB key is `(beta,alpha)`. The
  Cholesky vectors carry no spin index, so a block depends on the two orbital
  indices alone and the same entry serves either spin case.
* Key space ~`n_a*n_b + 2*C(n_s,2)`, far larger than any affordable cache
  (720 KB per block at 300 MOs), so reuse must be measured, not assumed.

The dgemm has a further inefficiency that is *not* part of this spec: it
recomputes the entries already held in `mo_integrals_cache` (only the rows below
`mo_integrals_cache_min` are missing). Restricting the block to those rows is a
separate, non-bit-identical change (different dgemm shape).

Alternative with no code change, and why it does not work: raise
`mo_integrals_cache_shift` until the window covers `mo_num`, which removes the
GEMM branch entirely (a full window makes every path read the dense cache, and
that cache is built from the Cholesky vectors, `map_integrals.irp.f:79-98`).
Level 8 is 32 GB and covers 256 MOs; level 9 is 512 GB and covers 512. At 300
MOs neither is affordable per process, so the block path stays.

## Stage 0 - counters and the top-M coverage curve (no behaviour change)

**Implemented** in `mo_two_e_ints/pair_block_stats.f90`, a plain module (picked
up by `qp_create_ninja` as a static source, `:475-479`), plus one hook in
`map_integrals.irp.f` and two in `fci/pt2.irp.f`. A plain module rather than
IRPF90 providers: the counters are written from inside the OpenMP walk, so as
providers they would join the provider graph, be lazily built on first touch
inside a parallel region, and be broadcast in every program that links this
module.

Two hook pairs, one per entry point, both keyed on `(k,l)`:

```
  ! get_mo_two_e_integrals_ij, do_mo_cholesky branch (map_integrals.irp.f:666)
  call wall_time(t_block0)
  <the existing code, unchanged: GPU probe or dgemm/sgemm>   ! :672-739
  call wall_time(t_block1)
  call pair_block_stats_add(k,l,cholesky_mo_num,t_block1-t_block0)

  ! get_mo_two_e_integrals, after the banned_excitation early return
  stat_timed = pair_block_stats_on .and. do_mo_cholesky
  if (stat_timed) call wall_time(stat_t0)
  <the existing code, unchanged: dgemv, sgemv or window copy>
  if (stat_timed) then
    call wall_time(stat_t1)
    call pair_block_stats_add_col(k,l,stat_t1-stat_t0)
  endif
```

`nblk` and `ncol` are flat arrays indexed by `key = k + mo_num*(l-1)`, the same
key Stage 1 hashes, so no second encoding exists. Per-pair *times* are
deliberately not recorded: every block of a run is the same dgemm and every
column request the same dgemv, so the two class means (`tim_blk/tot_blk`,
`tim_col/tot_col`) are the per-pair costs, and one more `mo_num^2` array of
doubles would buy nothing. Two counters per pair, 8 bytes in all.

The block timing covers the whole `if (mo_cholesky_double)` body on purpose: a
Stage 1 hit replaces the entire branch, so the branch - GPU probe and copy-back
and the fp32 `out_array_sp -> out_array` conversion loop included - is the cost
to divide by the hit rate, not the dgemm alone. The column timing starts after
the `banned_excitation(j,l)` early return, which returns zeros without touching
an integral: those calls are not requests and must not inflate the reuse.

`wall_time` is **wall clock** (`utils/util.irp.f:323`, "the equivalent of
cpu_time, but for the wall time"), so the accumulated times are thread-seconds,
not CPU seconds; the two coincide at one thread. It is already called this way
from this directory with no `use` (`mo_bi_integrals.irp.f:62`, external
procedure). Two clock reads per millisecond of work is free.

The gate is `QP_PAIRBLOCK_STATS=1`, read **once by the driver**
(`pair_block_stats_init(mo_num, trim(ezfio_filename))` from `program pt2`, before
`call run`), not at the call site. `get_mo_two_e_integrals_ij` has nine other
call sites - `davidson/u0_hs2_u0`, `mp2`, `mu_of_r/f_psi*`, `casscf_cipsi/bielec`
- and `get_mo_two_e_integrals` has more; none of them reaches the report, so a
call-site env read would arm the counters for programs that cannot print them.
`pt2` itself is a standalone OpenMP program (no MPI in `fci/pt2.irp.f`), so a
per-process counter is the whole story. The counters are written under one
`!$OMP CRITICAL(qp_pair_block_stats)`, which is also where the periodic snapshot
runs, so a snapshot always reads a consistent set of counters and no other
thread is inside Fortran I/O.

Report: `pair_block_stats_report(walk_time, nthreads)` called from `run` after
the walk, with `walk_time` measured by `wall_time` around the
`pt2_stochastic`/`omp_selection` dispatch (`:96-102`) and
`nthreads = max(nthreads_pt2,1)` - the team size the walk used, which is the
right denominator. When the counters are off, the column hook is one logical
test per call and no clock read; the block hook still pays two `wall_time`
calls - ~40 ms over a walk's ~2x10^5 block calls, 0.0006 % - and
`pair_block_stats_add` returns at once. The map-mode run below is the
acceptance that neither perturbs the result.

It prints, and nothing more:

```
 PT2 pair-block stats: mo_num=<n> cholesky_mo_num=<m> (Stage 0 of ...)
   block  calls=<n> distinct=<n> thread-s=<t> = <s>% of the walk, <t> s/call
   column calls=<n> distinct=<n> thread-s=<t> = <s>% of the walk, <t> s/call
   union  pairs=<n> both=<n> block-only=<n> memory if all cached=<MB> MB
   already paid: <n> column calls land on a pair whose block is computed anyway; <t> s
     caching those <MB> MB costs nothing beyond the block the walk computes regardless
   column-only: <n> calls on <n> pairs; a pair needs <h> hits to pay for its block (...)
   table: <ezfio>/work/mo_pair_cache.txt
```

The "already paid" line is the cheap half of the decision: those column requests
cost nothing extra once the block the walk computes anyway is kept, and each is
served by a strided read of one column. The "column-only" line is the other
half: a pair that never sees a block request must repay its own dgemm, at
`ceil(tim_blk/tim_col)` hits.

It also writes `<ezfio_filename>/work/mo_pair_cache.txt`, **every
`QP_PAIRBLOCK_DUMP` seconds (default 60) during the walk**, not only at the end:
at 114 MOs the forced-Cholesky probe does not finish its first window of 256
positions in 8 minutes, so a run that is killed must still leave the curve
measured up to the kill. The file is `#`-commented header lines (including the
totals), then one line per pair with `nblk + ncol > 0` (`k l n_block n_col`),
then `# END`. It is rewritten in place, so a kill landing inside the write costs
the tail of the table but never the header; the reader takes what is there.
`QP_PAIRBLOCK_DUMP=0` restores a single dump at the end of the walk. It is the
**full** list, not truncated to a cap, so one measuring run serves any future
cap, and the ranking is left to the offline reader: there is no ordering
constraint between the cap and the measuring run. Text, not binary: it is read
once and must be inspectable.

Acceptance: nothing arithmetic was touched, so the printed PT2/rPT2/E+PT2 must
be bit-identical to the reference row below. `pt2` prints `no request` when the
cache window covers all MOs (level 7 with `mo_num <= 128`): the Cholesky
contraction never runs and there is nothing to cache at that size. That case was
run and the row stayed bit-identical, with `QP_PAIRBLOCK_STATS=1` set.

### Building and running Stage 0

The `.rc` is bash-only (`${BASH_SOURCE}[0]`); under zsh it falls into the wrong
`case` arm and leaves `QP_ROOT` unset. A new static source also needs the graph
regenerated, since `ninja` alone will not see it:

```
bash -c 'cd $QP_ROOT && source quantum_package.rc &&
         scripts/compilation/qp_create_ninja update &&
         ninja -f config/build.ninja $QP_ROOT/src/fci/pt2'
```

(the real targets live in `config/build.ninja`; the top-level `build.ninja` is a
stub). Then, for a serial reference run:

```
OMP_NUM_THREADS=1 QP_NTHREADS_PT2=1 QP_PAIRBLOCK_STATS=1 \
  $QP_ROOT/src/fci/pt2 <fixture>.ezfio
```

`QP_ROOT` must be exported even though the binary is invoked by absolute path:
after the walk, `nuclei.irp.f:219` opens `$QP_ROOT/data/list_element.txt`, and
with `QP_ROOT` unset the run dies *after* the whole stochastic walk and *before*
the `PT2`/`rPT2` row is printed, leaving `WALL`/`CPU`, a sample table, and
`rc=2`. Twenty-five `IRPF90_temp` copies of that routine sit under `$QP_ROOT/src`,
so grepping the tree for the error message finds a stale copy, not the one the
`fci` binary was built from. Export
`QP_ROOT` in whatever shell launches the run; `source quantum_package.rc` must be
run under `bash` for this, and sourcing it is not otherwise needed.

For a walk that cannot finish, kill it and read the snapshots:

* copy the ezfio first (`cp -Rc` on APFS, a 300 MB clone), and set
  `QP_PAIRBLOCK_DUMP=30` so the curve appears within a minute;
* check `<ezfio>/work/mo_pair_cache.txt` exists after ~90 s. If it does not,
  the counters were never armed (env var not exported, `QP_PAIRBLOCK_STATS`
  read by the driver only) or `work/` is missing - nothing else can fail
  silently, and a 4-hour cluster slot is too expensive to discover it at the
  end;
* `QP_NTHREADS_PT2` (default `nproc`) is the team size the report divides by;
  it must equal `OMP_NUM_THREADS`, so set both. More threads covers more walk
  positions per wall minute and sharpens the multiplicities, at the cost of
  memory-bandwidth contention inflating both `s/call`. Take the call counts
  and the pair-set sizes from the many-thread run, and the `t_block/t_column`
  ratio from a one-thread run;
* the multiplicities transfer to production because the `(k,l)` request
  sequence depends on the determinant set and the mask structure, not on
  `mo_num`; the pair-set *size* does not, so only the production run sizes the
  cache.

## Stage 1 - static hot-block cache (only if the Stage 0 curve pays)

New EZFIO variable, `mo_two_e_ints/EZFIO.cfg`, after `[mo_cholesky_double]`:

```
[mo_cholesky_pair_cache_mb]
type: integer
doc: Memory cap (MB) for a static cache of MO integral blocks keyed by the
     virtual orbital pair of the zero-order PT2 determinants. 0 disables.
interface: ezfio, provider, ocaml
default: 0
```

Providers in `map_integrals.irp.f` (chained `&BEGIN_PROVIDER`, same block):

```
&BEGIN_PROVIDER [ integer, n_pair_cache ]          ! 0 = disabled
&BEGIN_PROVIDER [ integer, n_pair_cache_padded ]   ! max(1,n_pair_cache)
&BEGIN_PROVIDER [ integer, pair_cache_nhash ]      ! 2*n_pair_cache_padded+1
&BEGIN_PROVIDER [ double precision, pair_cache, (mo_num, mo_num, n_pair_cache_padded) ]
&BEGIN_PROVIDER [ integer, pair_cache_key, (pair_cache_nhash) ]    ! encoded key, 0 = empty
&BEGIN_PROVIDER [ integer, pair_cache_slot, (pair_cache_nhash) ]   ! 1..n_pair_cache
```

Build block: `n_pair_cache = 0` unless `do_mo_cholesky` and
`mo_cholesky_pair_cache_mb > 0` and `work/mo_pair_cache.select` exists.

That file is **not** the Stage 0 snapshot. The snapshot is raw counts in index
order; the selection is a fitted artifact, produced offline by the 30-line
ranking in the decision rule below and written as one `k l` line per pair,
best first. Keeping the ranking out of Fortran costs nothing and buys three
things: no sort in the library, no marginal-value formula duplicated in two
languages, and a fitted set that can be inspected, diffed and regenerated
without recompiling. Stage 1 reads `n` lines, takes
`n = min(nlines, mo_cholesky_pair_cache_mb*131072/(mo_num*mo_num))` so a
mismatched cap or a stale file cannot overflow the array, and refuses
(`n_pair_cache = 0`) if any line is malformed. That bound is computed *before*
the array's dimension expression, which is the whole point of the chained
providers: the cap is enforced by arithmetic, not by catching an allocation
failure. IRPF90 allocates `pair_cache` from its declared bounds and there is no
`stat=` hook, so an over-large cap aborts like any other QP3 allocation - the
remedy is to lower it, and the startup line below says which of the two numbers
bound it.

The generator is one call, so the artifact is reproducible from the snapshot:

```
pair_curve.py SNAPSHOT --threads N [--prod-mo M] --cap-mb C \
             --emit-select $EZFIO/work/mo_pair_cache.select
```

`--emit-select` is the only addition the tool needs (about six lines): write the
ranked `(k,l)` pairs it already computes, best first, one per line, and nothing
else - no header, no comment, `mo_num` is already known to the reader. It emits
the first `slots` entries of `ranked` (`slots = cap_mb*131072/(mo_num*mo_num)`,
the same expression Stage 1 uses), so the file's *order* is the value order and
writing more pairs than the cap can hold is harmless: Stage 1's truncation then
drops exactly the least valuable ones. With `--prod-mo` it ranks on the projected
production multiplicities, with the `t_col`/`t_blk` extrapolations the decision
rule describes; without it, on the raw snapshot counts.

Filling, in two phases:

1. `!$OMP PARALLEL DO` over the `n` slots: `pair_cache(1,1,i)` computed through
   the uncached routine. Each iteration reads slot `i` of the selection file and
   writes only its own slice, so the phase is race-free by construction;
2. a **serial** loop over `i=1,n` inserting `key = k + mo_num*(l-1)` into the
   open-addressing table (linear probing, `slot = mod(key, pair_cache_nhash)+1`).
   Serial on purpose: two threads inserting into the same probe chain can both
   see the same bucket empty and one key is silently lost. The table build is `n`
   probes, microseconds against the GEMMs of phase 1, and serial order makes the
   final layout a function of the file alone - the build is deterministic and
   happens once.

The fill passes a locally declared `type(map_type) :: dummy_map` to the uncached
block routine. Naming `mo_integrals_map` there would ask IRPF90 to build the
sparse MO map, which Cholesky mode never builds and does not need; the map
argument is only read in the `do_mo_cholesky = False` branch, and the provider
already requires `do_mo_cholesky`.

Stage 0 writes the full table, so the cap and the fitted set can be chosen after
the measuring run. Nothing in Stage 1 depends on the snapshot format.

Where the fill is triggered, and where it must not be:

* the fill is paid by an explicit `PROVIDE pair_cache` in the pre-loop block of
  `fci/pt2.irp.f`, immediately after `PROVIDE psi_energy` (`:64`). That region is
  serial, before `call run`, so the build is outside the walk clock that Stage 0
  reports and inside `WALL` - which is where it belongs, since the A/B compares
  walks and the fill is a fixed startup cost;
* `pt2_omp_task` also gets `PROVIDE pair_cache` in its pre-loop block
  (`cipsi_utils/pt2_omp.irp.f:272`, next to `PROVIDE psi_det_hii ...`). The
  generated `provide_*` routines are lock-protected, so this is correct for any
  caller; in the `pt2` program it is a no-op because the pre-provide already ran.
  It matters when the cap is 0 (the default): then the body allocates one dummy
  slot, reads no file, and returns at once;
* the phase-1 `!$OMP PARALLEL DO` carries `IF(parallel_fill)` with
  `parallel_fill = .not. omp_in_parallel()`, computed before the directive. If a
  future caller with `n_pair_cache > 0` reaches the lazy path from inside a
  `!$OMP PARALLEL`, nested parallelism is off by default and the fill would
  silently run serial - correct but slow, and it would be slow exactly once, in a
  way nobody attributes to the cache. The `IF` clause makes that explicit rather
  than mysterious. It is one logical variable, not a second code path.

`pt2` is a standalone OpenMP program: the cache is process-local and pays no MPI
broadcast. The other executables that link `get_mo_two_e_integrals_ij` never
provide `pair_cache`, and with the default cap of 0 the split routines behave
exactly as today.

The provider prints one line when it is enabled, and nothing when the cap is 0 -
the same contract as Stage 0's report, and the only way to tell afterwards that
the run you are looking at actually had the cache it claims:

```
Pair cache: 1422 of 1422 blocks, 977 MiB, fill 5.4 s, 8144 pairs in the file
```

Four numbers: what was requested, what fits, what was allocated and how long the
fill took, plus the file's size so a stale selection is visible. `fill` is
measured around the two phases (both of them, so the serial insert is inside the
number rather than hidden), and it is *not* part of the Stage 0 walk time.

Out of scope by design: the iterative CIPSI driver. It does not fill, and
`n_pair_cache` stays 0 there unless someone sets the cap by hand. The block set
is fitted offline from one snapshot of one determinant list, and in a CIPSI run
that list grows every iteration, so the set would be stale after the first one;
re-fitting per iteration would put an I/O and a ranking inside the SCF-like loop
to serve a single walk whose reuse is already captured by the deterministic
hit-anyway order within that walk. Stage 1 is a measured experiment on a frozen
wavefunction, and the production-size measurement is a `pt2` run on a
production wavefunction - which is the same program, so nothing is lost by not
touching the driver.

Split both routines:

* `get_mo_two_e_integrals_ij_uncached(k,l,sze,out_array,map)` = today's body,
  unchanged, all branches (GPU, dgemm, sgemm, map).
* `get_mo_two_e_integrals_ij(k,l,sze,out_array,map)` = `PROVIDE pair_cache`,
  then: probe; hit -> `out_array(1:sze,1:sze) = pair_cache(:,:,slot)`; miss ->
  call the uncached routine. `sze == mo_num` is already asserted (`:658`).
  The `banned_excitation` test is not in this routine; it stays where it is.
* `get_mo_two_e_integrals_uncached(j,k,l,sze,out_val,map)` = today's body of the
  column routine, unchanged, with its `banned_excitation(j,l)` early return
  first (it returns zeros and needs no block).
* `get_mo_two_e_integrals(j,k,l,sze,out_val,map)` = `banned_excitation(j,l)`
  test **first**, returning zeros exactly as today, then the probe; hit ->
  `out_val(1:sze) = pair_cache(1:sze,j,slot)` - a strided read of one column,
  `mo_num` cache lines, 19 KB at 300 MOs; miss -> the uncached routine.

The banned test is duplicated on purpose. It must run before the probe, because
a banned `(j,l)` returns zeros and *not* column `j` of the block: a cache hit
that skipped it would return the true integral. Keeping it in the uncached body
too is what lets that body stay a verbatim copy of today's routine; the test is
two array loads and a branch.

Both routines use the same `slot` table, so a block computed for `get_d0` is a
hit for the `get_d1` columns that follow it in the same mask, and vice versa.

Properties:

* Bit-identical by construction in both precision modes: the block is always
  produced by the one existing routine, with the same operands and shape. No
  second implementation of the integral computation. The fill runs inside the
  serial table-build region before the walk - the same place as
  `PROVIDE mo_two_e_integrals_jj` (`fci/pt2.irp.f:59`) - and every fill dgemm is
  issued from one thread with the BLAS thread count pinned to 1 by
  `qp_blas_env_init`, exactly like a walk dgemm; a hit therefore returns the bits
  a miss would have returned. The `!$OMP PARALLEL DO` of phase 1 is harmless
  under the launch contract (`OMP_NUM_THREADS=1`, walk team via
  `QP_NTHREADS_PT2`), where it degenerates to one thread anyway.

  "A hit returns what a miss would have returned" compares a *block* gemm with a
  *block* gemm, so it holds by construction. What it does not settle is whether
  the walk's **column** calls, which today stream a dgemv, see the same bits when
  they read the block instead. That is the one non-trivial claim in Stage 1, it
  is measurable on today's binary, and it should be measured before any code is
  written: `mo_integrals_cache_shift = 7` turns the existing window into exactly
  such a cache (every column call served from a gemm-filled block), so a
  completed walk at shift 7 against a completed walk at shift 6 is the A/B. See
  Validation for the first half of that measurement.
* The array is read-only after the build: no locks, no eviction.
* On a block hit the cost is one `8*mo_num^2` byte copy: 720 KB at 300 MOs,
  ~50 us, about 1-2 % of the 3.6 ms dgemm it replaces (`t_blk` measured 0.93 ms
  at 114 MOs, scaling as `mo_num^2*nvec`). On a column hit it is a strided read
  of one column - `mo_num` cache lines, 912 B but 114 lines here, 19 KB at 300
  MOs - so 3-10 us: a hit saves 85-95 % of the 36 us call it replaces, not all
  of it. The Stage 0 saving table is an upper bound by that margin, and that is
  the one number Stage 1 should re-measure rather than assume.
* A miss is exactly today's cost plus a hash probe. A hit is a hash probe, one
  branch and a copy: no locks, no atomics, no per-thread state.

Not in scope: LRU/eviction, a traced key sequence, a runtime-adaptive hot set.
If the static set turns out to cover too little, the fork is not more caching
machinery - it is the batched device GEMM discussed under the decision rule.

## Stage 0 result (benzene cc-pVDZ probe, 20 November 2026)

Four runs of the forced-Cholesky fixture (`benzene_chol.ezfio`: 114 MOs, rank
2694, window 1-64), same binary, `QP_PAIRBLOCK_STATS=1`,
`QP_PAIRBLOCK_DUMP=30`, each stopped at a fixed time budget rather than at the
end of the walk (the walk does finish: ~1.5-2 h at shift 6 and ~28 min at shift
7, one thread - see Validation). The 8-thread runs maximise the number of walk
positions covered; the 1-thread runs give uncontended per-call costs.

| wall | thr | block calls | block share | column calls | column share | pairs | both |
|---|---|---|---|---|---|---|---|
| 60 s  | 8 | 11 222  | 11.1 % | 1 635 761  | 81.7 % | 7 979 | 3 414 |
| 570 s | 8 | 79 149  | 10.2 % | 13 667 186 | 84.1 % | 8 144 | 6 132 |
| 480 s | 1 | 46 679  |  9.0 % | 10 973 097 | 83.2 % | 8 144 | 5 470 |
| 870 s | 1 | 90 299  | 10.5 % | 17 944 922 | 78.5 % | 8 144 | 6 412 |

The last row ran alongside two other single-thread probes for its final minutes,
which is why its column share is lower and "everything else" higher; its counts
are the deepest of the four and the per-call costs match the 480 s row.

Per-call cost, one thread: `t_blk` 0.93-1.55 ms (70 MFlop, 45-75 GFlop/s),
`t_col` 36 us (1.08 MB, 30 GB/s), ratio **21-26**. With eight threads the block
inflates 4x and the column 8x (bandwidth), so the ratio stays near 21: the
shares are robust, the absolute times are not.

What the numbers say:

* the `sample` profile is confirmed by an independent instrument: 9-11 % block,
  82-84 % column, i.e. 92 % of the walk in the two contractions;
* the pair set **saturates**: 7979 pairs after 60 s, 8144 after 480 s, and it
  does not grow afterwards while the call counts keep climbing. 8144 is 63 % of
  the 12 996 ordered pairs of 114 orbitals. Coverage therefore *converges* with
  walk time instead of decaying - the opposite of the staleness worry, at least
  within one wavefunction;
* 60-80 % of column calls land on a pair whose block the walk computes anyway
  (69 % at 60 s, 80 % at 570 s, 83 % at 870 s), for 3.4-6.4 k blocks;
* `block-only = 0` in all four snapshots: the block key set is a strict subset
  of the column key set. No cached block serves a pair that only the block
  routine wanted, and no column pair is unreachable by the block routine;
* marginal-value curve, corrected accounting (the block term included), at 1
  thread after 870 s, as fractions of the **integral** time: `M = 1000` (104 MB)
  removes 35.5 %, `M = 2000` (208 MB) 53.9 %, the whole set (847 MB) 88.1 %. The
  deeper eight-thread run gives the converged version: 41 %, 60 %, 93 %. The
  walk-time fractions are these times `T_int/T_walk`, about 0.80 (derived
  below), i.e. 28 %, 43 % and 70 %. The full-set figure is the consistency
  check: it must equal column share + (1 - pairs/blocks) x block share = 78.5 +
  (1 - 6412/90299) x 10.5 = 88.3 %, which it does.

The one input that curve had to *assume* - that a cache hit is free - is now
measured on the unmodified binary, because `mo_integrals_cache_shift = 7` turns
the existing dense window into exactly a Stage 1 cache on this fixture (all 114
orbitals inside the window, so every column request is a read out of a
gemm-filled array). One thread, fp64, 420 s into a walk:

```
 shift 6   column  17 944 922 req / 683.10 s = 38.07 us   block 90 299 req / 91.19 s = 1.010 ms
 shift 7   column  40 174 979 req /  16.04 s =  0.40 us   block      0 req
```

A hit is **0.40 us, 1 % of the streamed column**, and it is a 912-byte memcpy
from a prefilled `mo_num x mo_num` block - the same operation a sparse cache
would do, since column `j` of a stored block is contiguous. The shift-7 run also
pays the whole pre-fill up front - 114 x 114 = 12 996 blocks, the same
114 x 114 x 2694 product per block as the walk's own block calls, so ~1.010 ms
each, ~13 s, once. So the curve's 88 % is not an upper bound any more: an
all-resident cache really does remove 88 % of the integral time, block term
included. It removes nothing else, which is why the walk-time fraction is lower -
and the completed-walk proxy measures that remainder directly, at 1530 s spent
on ~93 % of the calls with only 44 s of it in the column path. Same fixture,
same binary, same call sequence: the residual is not-integral work, and it is
what `T_walk - T_int` means below.

The arithmetic that ties the two measurements together, one thread, shift 6:

```
T_walk = C_total x t_col + B_total x t_blk + T_other
       = 1.3x10^8 x 38.07 us + 6.6x10^5 x 1.010 ms + T_other
       = 4950 s              + 660 s               + T_other
T_other, from the shift-7 proxy (1530 s, 93 % of the calls, 44 s column path) = 1510 s
T_walk ~ 7100 s (2.0 h),   T_int/T_walk ~ 0.79
```

`C_total` is the measured 1.3x10^8 (the shift-7 run issued 1.278x10^8 for 93 % of
the positions), `B_total` follows from the head's one-block-per-198-columns
ratio, and `T_other` is the independent number. The head's own rate gives a
lower `T_walk`, ~6300 s, and the gap is explained by `T_other` being mostly per
walk *position*: a head rate weighs the call-rich positions more than the whole
walk does. Arm A settles it and reports its own split, turning `T_int/T_walk`
into a measurement. The curve conversion above uses 0.8 and is insensitive to
the choice: 35.5 % of the integral time is 28-32 % of the walk either way.

What does not transfer, and must be measured at production size:

* memory per block scales as `mo_num^2`: 104 KB here, 720 KB at 300 MOs, so the
  same `M` costs 7x more;
* the pair set is a property of the determinant set plus `mo_num`, not of the
  window, so 8144 says nothing about benzene/cc-pVTZ;
* the window covers 64 of 114 orbitals here and 128 of 300 at production. That
  does not change *which* calls take the dgemv branch - every column call needs
  the MOs outside the window, whichever side they are on - it changes how many
  MOs are outside: 50 here, 172 at 300 MOs with `shift = 7`. So `t_col` grows as
  `nvec*(mo_num - 2^shift)` while `t_blk` grows as `nvec*mo_num^2`: the ratio
  `t_blk/t_col` rises roughly with `mo_num`, 25 here and ~50 at 300 MOs. The
  column share is `r/(r + t_blk/t_col)` in the call ratio `r`, so at the head's
  `r = 201` (measured: 17.5x10^6 column calls and 87 044 block calls) the 89 %
  becomes ~80 % at production, before any word about the determinant set. The
  production run must therefore report both the share and the call ratio; the
  share alone would hide which of the two moved.

The decision this supports: at a 1 GB cap (1422 blocks at 300 MOs) Stage 1
removes 35-55 % of the integral time - 28-43 % of the walk - *if* the production
pair set is as concentrated as this one; the production curve is the only thing
that settles it. Nothing here justifies a larger cap than ~2 GB per process, and
nothing justifies eviction.

One model did not survive contact with the data, and the design does not need it.
The obvious byte model - `t_col` linear in the number of MOs outside the window,
`t_col = 0.35 us + 0.755 us x n_out` from the two long runs (0 and 50 columns
streamed) - predicts 62 us at shift 5 (82 outside) and 74 us at shift 4 (98), and
the measured means are 82 and 66: the wrong way round.
The means are not comparable. The walk visits generators in `pt2_J` order, i.e.
most expensive first, so a fixed *wall time* samples a different *prefix* at each
shift, and those two cases only walked 150 s (the probe script subtracts its
300 s snapshot wait from the budget, whatever the snapshot actually took). Only
the two long runs - 870 s and 1530 s - average over comparable stretches.
So `t_col` at production size stays an extrapolation (~73 us at 300 MOs,
`shift = 7`), not a measurement. It does not gate the design: the cap ranks
pairs by `n_col*t_col + (n_block-1)*t_blk` with `t_col >> t_hit`, and both
per-call costs scale the same way with size, so the *order* is insensitive; only
the headline percentage moves, and the end-to-end A/B bounds that directly.

For the record, the 8-thread and 1-thread means differ by 4-8x while the shares
agree to 2 points: any future Stage 0 run can be done at full thread count for
coverage, and only the `t_blk/t_col` ratio needs a serial run.

## Stage 0 at production size: where, what it costs, what it decides

Two of the measured quantities do not transfer. The pair-set size is a property
of the determinant set and `mo_num`, and the per-call ratio is a property of
`cholesky_mo_num`. What does transfer, and what made the 114-MO probe worth
running, is the *shape* of the request sequence: the `(k,l)` multiplicities
depend on the mask structure, not on the basis.

The fixtures are cluster-only. `cholesky_mo_transp` at 510 MOs and ~11 500
vectors is 23.9 GB (12 GB in the fp32 twin), and 4.2 GB at 264 MOs; the local
machine has 8 GiB of RAM. So the run happens where the fixture lives.

Naive scaling says nothing useful about the split: per call the block grows as
`nvec*mo_num^2` and the column as `nvec*mo_num`, but the two classes come from
different selector subsets (`nt == 2` vs `nt == 3`), and their counts move
independently with the determinant set. At 114 MOs the ratio of the growth
factors already points downward for the block; whether it does is the first
thing the production run prints.

Both routes start with a ten-minute inventory on the login node (no job):

```
ls -l ~/qp3/src/fci/pt2 ~/qp2/src/fci/pt2          # is there a binary already?
for f in mo_num cholesky_mo_num do_mo_cholesky do_pt2 ; do
  echo $f $(cat <fixture>/mo_two_e_ints/$f) ; done  # which fixture is usable
which perf ; cat /proc/sys/kernel/perf_event_paranoid
du -sh <fixture> ; df -h $WORK
```

Route A, the proxy: if a `pt2` binary already exists, profile it on the 264-MO
fixture with `perf record -g` and read the self-time share of `dgemv` and
`dgemm`. No build, under an hour. It answers the one question that can kill the
cache thesis - does the column stream still dominate at 264 MOs - and gives no
pair sets.

Route B, the definitive one: build QP3 with the in-tree ARM config and run the
counters. `config/gfortran_armpl.cfg` already targets ARM + armpl
(`-larmpl_lp64_mp`). Three traps: its `-Ofast` enables `-ffast-math`, so such a
build is not bit-comparable with the laptop row (use `-O2` when that matters,
and it does for any `pt2_data_I` comparison); HOME is a 10 GB quota, so the
tree goes in `$WORK`; and compute nodes are offline, so the clone happens on
the login node. Carry the branch as `git archive pt2-pair-block` rather than
pushing it.

Then one slot:

```
OMP_NUM_THREADS=1 QP_NTHREADS_PT2=40 QP_MAXMEM=400 \
QP_PAIRBLOCK_STATS=1 QP_PAIRBLOCK_DUMP=60 $QP_ROOT/src/fci/pt2 <fixture>
```

`QP_NTHREADS_PT2` is the walk team and is free (`src/fci/BENCHMARK.md`); the OMP
default must stay 1. Kill the job once a snapshot covers a few million column
calls - the counters are rates - and read it with `pair_curve.py` from the
measurement scratch (carry the file along; it is not in the tree).

The run decides, in this order:

1. the column share at 264 and 510 MOs. Below ~40 % the cache thesis is a
   benzene artefact and Stage 1 is dropped, together with the bandwidth half of
   the device case;
2. the union pair-set size and its marginal-value curve, which fix `M` and the
   cap;
3. `t_blk/t_col` at production `nvec`, the break-even of the column-only pairs.

One caveat about the times, from the probe: at eight threads both classes
inflate almost linearly in the team size (36 us becomes 280 us), which is the
signature of a saturated aggregate bandwidth - 8 x 1.08 MB / 0.28 ms = 31 GB/s,
the whole machine. On a production node the column stream is therefore not just
a per-thread latency but a node-wide byte budget, and a hit that reads 19 KB
instead of 2 MB removes bytes rather than hiding latency. That is the strongest
version of the cache argument, and it is also why the host cache beats sending
the same stream to a device.

## The no-memory alternative: group the column calls by `k`

The cache is not the only way to cut the column stream, and at production size
it may not be the best one. `get_d1` issues its eight `get_mo_two_e_integrals`
calls with **one and the same column index**, `j = hfix` (`selection.irp.f:1305,
1306, 1349, 1350, 1402, 1403, 1439, 1440`), and with `(k,l)` drawn from ordered
pairs of the three orbitals `{p1, p2, pfix}`:

```
  (p1,p2) (p2,p1)      (pfix,p1) (pfix,p2)      (p1,pfix) (p2,pfix)
```

Six distinct keys, all with the same `j`, over **three** distinct `k`-slabs,
eight times through them. (The `sp /= 3` branch issues them from the loop
`do i=1,3` at `:1399-1403`; the `sp == 3` branch from the six explicit sites.)
The dgemv operand `cholesky_mo_transp(1:nvec, max+1:mo_num, k)` depends on `k`
only, so calls sharing `k` stream the same contiguous slab - 1.08 MB at the
cc-pVDZ probe, ~2 MB at 300 MOs with ~1500 vectors - to produce 50 and 172
doubles respectively.

Grouped by `k`, the same six columns cost three gemm passes with two columns
each (`out(i,s) = sum_P transp(P,i,k)*transp(P,j,l_s)`), a 2x traffic cut on the
77 % stream, with no memory, no fit, no staleness, and no cap to choose. It
applies to every column call, so it multiplies with Stage 1 instead of competing
with it.

Two costs, both real:

* the arithmetic of each output element is unchanged (a dot product over `P`,
  the same cancellation), but BLAS does not promise that a gemm and a gemv
  return bit-identical results - the `K` blocking differs. That caveat is now
  measured rather than assumed, and it came out favourably: on Accelerate at 114
  MOs, a walk served entirely from gemm-filled blocks reproduces the E0 of the
  window-gemm-plus-streamed-dgemv walk to all 16 digits (see Validation). Same
  caveat, but no longer a reason to prefer one fork on faith.
* eight call sites become three, and the results are needed in
  `hij_cache1(1:sze,1:2)` with the two columns of a `k`-group not adjacent in
  today's code, so the plumbing of `hij_cache1` has to be re-indexed. That is a
  refactor of a 200-line unrolled kernel, not a wrapper.

A zero-risk partial version exists and is worth measuring before either fork:
reorder the eight call sites so that calls sharing `k` are adjacent, and let the
slab stay in L2/L3. That needs no arithmetic change at all. It is only
plausible because the slab is 1-2 MB per thread; at the production thread count
with 8 MB of slabs in flight it is exactly the kind of assumption that stops
holding on the machine you did not test on, so it must be measured, not assumed.

## Decision rule

Stage 0 reports, per class (`block`, `column`), the number of requests, the
number of distinct pairs, the thread-seconds and hence the share of the walk;
plus the size of the shared key set. The decision is then bookkeeping:

* A pair with `n_block >= 1`: the fill computes the block once, and from then on
  every block call is a copy, so caching saves `(n_block-1)*t_blk + n_col*t_col`
  for one slot. Never negative, and not small: the probe's block class averages
  **12.9 calls per pair** (79 149 calls on 6 132 pairs), so e.g. 11 ms of GEMM
  per pair before the first column is counted. Today's walk pays the block
  `n_block` times; the cache pays it once.
* A pair with `n_block = 0` (column-only): caching it means computing a dgemm
  the walk never asked for, which pays back after `ceil(t_blk/(t_col-t_hit))`
  hits - 27 on the probe, where `t_hit = 0.40 us` is now measured rather than
  assumed (Validation), so the earlier `t_blk/t_col` form overstates the
  threshold by 1 %. The 2674 column-only pairs of that run average 1647 calls
  each, so the threshold is not close: in this regime column-only pairs are as
  profitable as the others.
* The marginal-value curve is therefore: rank the pairs by
  `n_col*t_col + (n_block-1)*t_blk` for `n_block >= 1` and `n_col*t_col - t_blk`
  for `n_block = 0`, walk the ranking, keep the prefix that fits
  `cap/(8*mo_num^2)` blocks, and the saving is the sum over that prefix. The
  report prints the ingredients; `pair_curve.py` (in the measurement scratch,
  not in the tree) is the 60-line version of that ranking. No eviction, no
  adaptivity, no runtime decision.
* `block-only = 0` on the probe means the two key sets nest; expect that to
  hold, but the rule does not depend on it.
* If the union of the production key set fits the cap outright, there is no
  threshold to choose: cache everything and the design shrinks further.
* If both classes are small on the production run (a few % of the walk), stop
  regardless of the curve: the target is elsewhere (per-mask `mat = 0`, the
  `banned(:,:,1..2)` copies, the `get_double_excitation` re-derivation inside
  `diag_H_mat_elem_fock`).

If the curve instead shows reuse ~1 in both classes (every mask a fresh pair),
the contractions are genuinely irreducible *and* the case for a device becomes
the strongest one available - but not for the reason this spec started with. The
quantities are:

* the column stream is bandwidth-bound **on the host**: 2 MB read per call at
  300 MOs, ~30 GB/s on one core, ~15 GB/s per core when eight share the socket;
* `cholesky_mo_transp` is `cholesky_mo_num x mo_num x mo_num` doubles
  (`cholesky.irp.f:47`): **1.08 GB** at 300 MOs and 1500 vectors, 0.54 GB in its
  existing fp32 twin (`:225`). It fits in one A100-80 GB, and in the 12 GB of
  the laptop's RTX 3500 Ada;
* so a device does not accelerate the column stream, it **removes** it: the
  vectors are uploaded once per process (0.5 ms at 2 TB/s) and the 2 MB read
  becomes an on-chip read. What is left per call is the kernel launch, a few us
  against the ~70 us the host spends;
* the block GEMM is the opposite case and the one a device genuinely accelerates:
  0.18-0.36 GFlop, compute-bound, 200x the fp64 throughput of one host core.

Two consequences for the fork, both worth writing down before opening it. A host
cache and a resident device copy attack the same bandwidth problem, so they are
alternatives rather than complements; only the block GEMM is unique to the
device. And Vulkan - the portability route - has no fp64: Apple's GPUs have none
at all, and most consumer devices expose `shaderFloat64` as false. The Vulkan
path is therefore fp32, which is exactly what `mo_cholesky_double = False`
already selects (`EZFIO.cfg:38-42`), and why an fp32 host build is needed as the
oracle before any device number is believed.

Caveat for the report: the hot set is fitted to one wavefunction. In a CIPSI
loop it must be refreshed by re-running Stage 0, and the curve's decay with `M`
is the only evidence about how fast such a set goes stale.

## Validation (and why benzene is not a Cholesky test bench)

The cc-pVDZ benzene fixture carries 2985 AO Cholesky vectors (rank 2694 in the
MO basis) for 114 MOs, and it is *map* mode (`do_mo_cholesky = False`,
`mo_integrals_cache_shift = 7`, window 1-128, so it covers all 114 MOs). Forcing
the Cholesky path on it is not the hopeless case an early estimate suggested -
"~80x map mode, the first 256 walk positions do not complete in 8 minutes" came
from a cold window that was still paying the provider chain and the cache fill,
and it survives in no table here. Measured, one thread, fp64, on this fixture: **20 600 column calls/s at shift 6**
(870 s of walk) and 83 500/s at shift 7. The shift-7 rate *decays* along the walk
(93 900/s at 600 s, 87 700 at 1110 s) because `pt2_J` orders the generators most
expensive first, so the only trustworthy completion number is the direct one:
killed at 1530 s, that walk had issued 1.278x10^8 calls and was ~7 % short of the
end - `C_total > 1.3x10^8`, ~1700 s. At shift 6 the same count at the
head-of-walk rate is ~2.0 h; budget 4 h, expect 1.5-2.2 h. Both fit on an idle
laptop overnight, so the fixture validates mechanics, coverage *and* end-to-end
bit-identity without production size; the coverage curve transfers to production
because the `(k,l)` request sequence does not depend on the window. Only the
per-block memory and the pair set do.

The same 1530 s spent 44.1 s in the column path - **2.9 %** - so with free
integrals the walk is bounded below by its non-integral work, and the Stage 0
shares (89 % at shift 6) are shares of the expensive *head*, not of the walk.
The end-to-end saving is therefore capped in the 75-78 % range, not 88 %: the
88 % is a fraction of the integral time, and ~22 % of the walk is not integral
work at all. That gap is exactly what the completed A/B measures, and why it is
worth two overnight runs rather than a model.

The acceptance tests that follow from all this:

* Stage 0 mechanics: any run with requests > 0 (fixture above, killed after a
  few minutes), then inspect `work/mo_pair_cache.txt` and the summary block.
  `no request` at level 7 is expected on 114 MOs and is not a bug.
* Stage 1 bit-identity, mechanical: in the hit path, **on the first hit only**,
  recompute the block through `*_uncached` into a scratch buffer and compare bit
  patterns (`transfer()` to `integer(8)`); on mismatch `qp_bug(irp_here,...)`.
  One extra dgemm per run, and it executes at production scale where it matters.
* Stage 1 end-to-end: run with the cap at 0, then with the cap set, on the
  same binary and fixture; compare the printed row (`PT2`, `rPT2`, `E+PT2`,
  `E+rPT2` with their errors) to all digits. The `PT2 Energy denominator` line
  must agree too, but it cannot move: it is `<Psi|H|Psi>` computed from
  unchanged providers before any fill, so a difference there means the cache
  corrupted shared memory, not that the cache changed an integral - which is
  still worth failing on.
  That is the only artifact that exists today: `pt2_data_I` is a local
  allocatable inside `pt2_stochastic` (`pt2_omp.irp.f:246`, allocated `:334`)
  and is deallocated at `:482` without ever being written out. A
  per-determinant dump - ~10 lines before that `deallocate`, behind an
  environment variable - would turn the comparison from an aggregate into a
  per-generator one; it is worth writing only if the row ever disagrees, since
  a 1-ulp change anywhere in the walk moves one of the 17 printed digits. The
  walk is schedule-independent by contract (`pt2_omp.irp.f` header), so the
  A/B can be done at the production thread count.

The walk cannot be shortened by moving `pt2_relative_error`: `do_exit` needs
`c >= 10` combs *and* convergence (`pt2_omp.irp.f:725`, `:800`), and the
certified run ended at the last generator with `conv = 2.1x10^-3` against its
2x10^-3 target - so a looser target saves at most 2 of its 12 combs. The walk
length is fixed by the comb construction (`pt2_R`), not by the target, which is
also why the same `C_total` can be assumed in both window settings.

### The one non-trivial claim, half-answered on the existing binary

Whether a column call served from a gemm-filled block returns the bits of the
streamed dgemv it replaces is the only Stage 1 claim that is not true by
construction. The probe fixture answers half of it without any new code: at
`mo_integrals_cache_shift = 7` the window covers all 114 MOs, so every column
request reads a block the dgemm fill produced - a Stage 1 hit for the whole walk.

```
 shift 6, window 1-64    E0 = -231.2007184830427    (window gemm + streamed dgemv)
 shift 7, window 1-114   E0 = -231.2007184830427    (every value from a gemm block)
```

What the printed line is matters, because it is the only integral-derived number
both runs produce before their walks. `PT2 Energy denominator` prints
`pt2_E0_denominator + nuclear_repulsion` (`cipsi/energy.irp.f:9-30`), and with
`h0_type = EN` that array is `psi_energy` - which the run *recomputes* as
`<Psi|H|Psi>` from the stored coefficients (`davidson/u0_hs2_u0.irp.f:3`), it is
not read from EZFIO. The proof is the representation gap: the same frozen
wavefunction gives -231.2015993197768 in map mode and -231.2007184830427 here.
So the agreement above is a genuine statement about integrals:

* at shift 6 the determinant-pair contraction reads pairs on both sides of the
  window - inside, the window gemm's values; outside, the streamed dgemv of
  `get_mo_two_e_integrals_i1j1`'s partial-window branch
  (`map_integrals.irp.f:788-812`);
* at shift 7 the branch is not taken and *every* pair is read out of the window
  array that `dgemm('T','N',114,114,2694)` filled, at `ldc = 16384`
  (`map_integrals.irp.f:76-92`).

Identical to all 16 printed digits, same fixture, same flags, fp64, one thread.
The terms of that contraction are O(10^2) Ha with alternating signs and the
result is -231.2, so the sum is heavily cancelling and a *systematic* or
*widespread* 1-ulp difference in the streamed columns would reach the last digit.
A single 1-ulp change in one integral would not - an aggregate is a weak
per-value test, which is why the 30-line program below stays on the list. What
the aggregate does rule out is the failure mode that actually threatens the
design: a wrong index or a misaligned column, i.e. an O(1) error on some pairs.

### The completed A/B, and the two numbers it produces

A completed walk looked out of reach when this section was first written - ~26 h,
from a window-based extrapolation that the measured per-call rates replaced. One
thread now needs ~2 h at shift 6 and ~28 min at shift 7: two overnight runs on an
idle laptop rather than a research project.

The walk issues the same call sequence in both window settings - the integrals
are bit-identical, so the determinant trajectory is too - and `C_total > 1.3x10^8`
column calls.

* **bit-identity.** The printed `PT2`/`rPT2` row must agree to all digits. That
  is much sharper than the E0 line above: PT2 is a nonlinear function of the
  ~10^8 matrix elements the walk evaluates - squares and `E0 - H_jj` denominators
  - and every one of them is an integral or a diagonal element built from
  integrals, so any per-value disagreement anywhere is amplified into the
  printed digits rather than averaged away.
* **the cache's worth.** `WALL` arm A against `WALL` arm B. The floor is
  measured, not modelled: arm B runs 1530 s for ~93 % of the calls with 44.1 s
  in the column path, so a full walk with free integrals costs ~1600 s, and no
  integral-side change can go below it. The ceiling saving is therefore
  `1 - 1600/T_A`, where `T_A` is arm A's time. Two extrapolations disagree about
  `T_A` and arm A settles it:
  - head-of-walk rate: 20 900 calls/s over 840 s at one thread, i.e. ~6300 s for
    `C_total = 1.3x10^8`;
  - integral + floor: `C_total*t_col + B_total*t_blk` = 4950 + 660 s, plus the
    1510 s of non-integral work the shift-7 arm measured, i.e. ~7100 s.

  They differ because the non-integral cost is mostly per walk *position*
  (46 824 of them) while `C_total` counts calls: the head, where the expensive
  generators issue many calls each, amortises that cost better than the tail
  does, so a rate measured in the head overstates the whole walk. Either way the
  ceiling is 75-78 % of arm A, and the "88 % of the two contractions" figure is
  unreachable *whatever* the cache does - it counted only the contractions.

What the A/B cannot give is a per-*value* statement, because it aggregates. For
that, a 30-line Fortran program against the real slab shape (`nvec` 2694, 114
columns, stride `nvec`): one dgemv for column `j`, one dgemm for the block,
compare `transfer()`ed bit patterns. Seconds to run, no QP3 dependency, on the
BLAS that will run the walk. It is the only test that compares the two numbers
one at a time.

Two shape details keep that A/B honest, and one of them is a design rule:

* at 114 MOs, shift 7 makes the window fill a `dgemm('T','N', mo_num, mo_num,
  cholesky_mo_num)` - the same M, N and K as the block routine's own call
  (`map_integrals.irp.f:712`), so arm B is the Stage 1 fill *shape*, not an
  analogy. The leading dimensions differ (`nwin**2 = 12996` for the window,
  `sze = 114` for `out_array`), and a BLAS is free to pick a different kernel
  for a different `ldc`: nothing forbids a 1-ulp difference, so the A/B reports
  what this BLAS does and cannot guarantee it on another;
* therefore the Stage 1 fill must not duplicate a dgemm - it calls the walk's
  own `get_mo_two_e_integrals_ij`, so the cached block *is* the array the walk's
  block call returns. Block identity is then by construction on any BLAS, and
  the only claim left is column-of-block versus streamed dgemv, which is exactly
  what the A/B and the 30-line program test.

One design consequence worth stating, because it explains why the key is the
pair and not the column: a cache keyed by `(k,l,j)` would be bit-identical *by
construction* (it stores the dgemv result it computed), needs no offline fit and
no fill GEMM, and is `mo_num` times smaller per entry. It is rejected because
the reuse is not there: 10.97x10^6 column calls at one thread land on 8144
distinct pairs, i.e. ~1350 calls per pair, while a given `(k,l,j)` triple is
requested by the masks that share `j = hfix` - the reuse is in the pair, which is
exactly what `n_block` and `both` measure. Stage 0 does not count the triple set,
so this is an inference from `both = 80 %`, not a measurement; adding a
`(k,l,j)` counter is one hash and would close it.

Do not compare across `do_mo_cholesky` modes. Map mode on this same fixture
gives E0 = -231.2015993197768 against -231.2007184830427 for either Cholesky
number, a 3.8x10^-6 relative difference with bit-identical orbitals, AO basis
and Cholesky vectors (md5-checked). It is a difference between integral
*representations* - one contraction plus a 1e-15 threshold and a symmetry key,
versus dgemm inside the window and dgemv outside - not an error of either, and
8.8e-4 absolute is comparable to the PT2 statistical error of the certified run
(+/- 6.4e-4). Cause not tracked down; the threshold cannot explain it. The
consequence is narrow and already respected below: bit-identity and the fp32
oracle are both *within* one representation, so the acceptance for Stage 1 is
cap-0 against cap-M on Cholesky mode, never Cholesky against the map-mode
reference row.

## Reference row (benzene, 1 thread)

First completed serial run of the fixture (`mo_num` 114, `n_det` 46880,
`do_mo_cholesky` False so the full level-7 window applies and no block GEMM
runs, `mo_cholesky_double` True, `h0_type` EN, `OMP_NUM_THREADS=1`,
`QP_NTHREADS_PT2=1`), with the `is_cfg` hoist applied:

```
 N_det = 46880   N_cfg = 9150
 PT2    = -0.32980666295167710   +/- 6.4074183210651380E-004
 rPT2   = -0.30532326708442498   +/- 5.9317597705745831E-004
 E+PT2  = -231.53140598272844
 E+rPT2 = -231.50692258686118
 WALL 1453.5 s, CPU 1067.0 s
```

Stage 1 is accepted when a run with the cap set reproduces `PT2`, `rPT2`,
`E+PT2`, `E+rPT2` printed to full precision against a run with cap 0 on this
same fixture. Two caveats: the fixture has `h0_type = EN`, so it never enters
the CFG branches - it validates mechanics and arithmetic, not the CFG path; and
a single run does not prove the walk deterministic, so the cap-0 row must be
regenerated by the same binary in the same session as the cap-on row rather
than trusted across sessions.
