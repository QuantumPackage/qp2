# Static hot-block cache for the Cholesky MO-integral block (PT2 selection)

Spec for implementation. Stage 0 is a measurement with no behaviour change and
decides whether Stage 1 is written at all. It said yes: Stage 1 is written,
built, and gated (*Stage 1 gate at 264 MOs*).

## Where this stands (21 November 2026)

* **Stage 0 is done and committed** (`pair_block_stats.f90`, opt-in via
  `QP_PAIRBLOCK_STATS=1`), and it has been run to the *end* of a walk, not just
  into its head: 114 MOs, one thread, fp64, `shift = 6` -> 7673.7 s of walk,
  1 447 862 block calls (1273.7 s, 16.6 %), 139 365 752 column calls (4777.2 s,
  62.3 %), 8144 pairs of which 8100 are in both classes, `block-only = 0`.
  **78.9 % of the walk is these two contractions**, and 139 349 290 of the
  column calls land on a pair the walk fills a block for regardless.
* **The ceiling is measured twice, and it is 80 %.** The same binary with
  `mo_integrals_cache_shift = 7` - every column served from a gemm-filled block,
  i.e. Stage 1 with an infinite cap - completes in 1524 s against 7739 s:
  **5.08x, 80.3 % saved**, with 1438 s of non-integral work that no integral-side
  change touches. See *Validation*.
* **At this fixture size that one EZFIO line is the whole win**, and Stage 1
  lands on it rather than above it: `pair_cache(i,j,slot)` has `i` as the fastest
  index, so a column hit is a *contiguous* read of `mo_num` doubles - the same
  copy, at the same width, that the shift-7 path does, and the measured cost of
  that path whole-call is 0.293 us. Predicted: 1438 s floor + 139.4M hits at
  0.29-0.40 us + ~7 s of fill = 1485-1500 s, which is the shift-7 walk's own
  number. Stage 1 earns its code only
  where the window cannot cover `mo_num`, which is production: at 300 MOs,
  `shift = 7` is 550 GB (see *Stage 0 at production size*).
* **The 264-MO probe changed the design, and simplified it.** At 264 MOs the
  walk asks for 235 888 blocks covering 1833 distinct pairs - **129 recomputations
  per pair**, 99.2 % of the block time - and 4.68M columns over 1844 pairs. The
  whole pair set is settled after 135 s of a 4226 s walk and it is 975 MiB, so
  the cache fills itself on first miss, with no selection file and no pre-walk
  pass, and the cap is the only knob (*Stage 1*, revised below). One caveat came
  out of the same run: the walk settles its *membership* early but says nothing
  about *ordering*, so an offline ranking is still the right answer if a
  production cap cannot hold the whole set.
* **The loop-swap alternative is answered, and it loses.** A block's key is a
  property of the interesting determinant, not of any loop that iterates it, so
  there is no invariant block to hoist; inverting onto the pair moves the cost
  from an 803 MiB table of blocks to a resident set of partial `mat` arrays of the
  same order, and buys only the 47 % block half. The same read did turn up a
  loose end worth ten lines: `get_d0` mallocs a full 545 KiB block on each of its
  235 888 calls (*The one hoist that does pay*) - bounded at 0.2 % of the walk, so
  hygiene rather than a rung on the way to the cache's 8.6x.
* **Replacement, not freeze, is specified for the capped regime** (*When the cap
  binds*), and it is the last thing to implement: it is a thrash bound, not the
  value ranking, and it costs a seqlock on the hit path.
* **Stage 1 is implemented, built and gated at 264 MOs** (`map_integrals.irp.f`,
  cap deduced from `qp_max_mem`, `QP_PAIR_CACHE_MB` overriding it): ~9x CPU and
  8.6-9.9x wall on the 264-MO fixture,
  `rPT2` one ulp from the cap-0 row against a 1e-9 gate, both hit classes
  asserted in-process, the freeze branch exercised by a 512 MiB cap, and the
  deduced cap reproducing the best capped arm **bit for bit** - all nine printed
  numbers, `diff`-clean. What the laptop still cannot supply is the production
  pair set - the pairs touched by a
  production determinant list at production `mo_num` - and how much of the block
  and column mass fits a 1-2 GB cap. That needs a Stage 0 run on a production
  wavefunction, which needs a QP3 build and a fixture on a cluster - a separate
  effort, deliberately set aside. What the 264-MO run does supply is that the
  *shape* of the answer does not depend on the window: the block path is 129x
  redundant there, against 179x at 114 MOs, and it is the larger half of the
  integral time.
* **The cap is deduced, not configured** (21 November 2026). `pair_cache_mb` is
  half of `qp_max_mem`, QP's own memory target, which is the share
  `cholesky.irp.f` already gives its temporary AO buffer; the arrays are not
  zeroed by the allocate, so the cap bounds the fills and does not reserve the
  pages. The deduction reads the budget and never the memory free at that
  instant: below the demand the admitted set depends on the slot count, and a
  slot count that moved with the machine's load would move the PT2 row in its
  last ulp between two runs of the same input. `QP_PAIR_CACHE_MB` remains the
  per-run override, and 0 disables the cache, which is what keeps the three gate
  arms reproducible.
* The device fork that this spec opened was costed on an A100 and **closed by
  measurement**; the numbers are in `GPU_ROUTE.md` in the measurement scratch,
  and the short version is under *Decision rule*.

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
   block  calls=<n> distinct=<n> thread-s=<t> = <s>% of the walk, <us> us/call
   column calls=<n> distinct=<n> thread-s=<t> = <s>% of the walk, <us> us/call
   union  pairs=<n> both=<n> block-only=<n> memory if all cached=<MB> MB
   already paid: <n> column calls land on a pair whose block is computed anyway; <t> s
     caching those <MB> MB costs nothing beyond the block the walk computes regardless
   column-only: <n> calls on <n> pairs; a pair needs <h> hits to pay for its block (...)
   table: <ezfio>/work/mo_pair_cache.txt
```

The "already paid" line is the cheap half of the decision: those column requests
cost nothing extra once the block the walk computes anyway is kept, and each is
served by one contiguous read of a column. The "column-only" line is the other
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
  memory-bandwidth contention inflating both per-call means. Take the call counts
  and the pair-set sizes from the many-thread run, and the `t_block/t_column`
  ratio from a one-thread run;
* the multiplicities transfer to production because the `(k,l)` request
  sequence depends on the determinant set and the mask structure, not on
  `mo_num`; the pair-set *size* does not, so only the production run sizes the
  cache.

## Stage 1 - a lazily filled block cache (only if the Stage 0 curve pays)

The cap is deduced, and `QP_PAIR_CACHE_MB` only overrides it. By default
`pair_cache_mb` is half of `qp_max_mem` - QP's own target, physical RAM minus
1 GB, or `QP_MAXMEM` when set - which is the same half `cholesky.irp.f`
already gives its temporary AO buffer (the two do not overlap in time). The
deduction reads the budget and never the memory free at that instant: once the
pool is smaller than the demand, the admitted set depends on the slot count, so
a slot count that moved with the machine's load would move the PT2 row in its
last ulp from run to run. `QP_MAXMEM` is how a shared machine is told to take
less; `QP_PAIR_CACHE_MB` overrides the cap per run, and 0 disables the cache.

`qp_max_mem` is `MemTotal - 1 GB` on Linux, so on a cgroup-limited node it is
the node's RAM and not the job's allocation, exactly as it already is for the
AO Cholesky buffer and every Davidson sized from it. A job under a memory limit
sets `QP_MAXMEM` to that limit, the cache follows, and nothing new is
introduced. What the cache adds under that misreading is virtual address space
and not resident memory: the pool is not zeroed and its pages are touched only
by the blocks it actually fills. If the commit fails anyway, the `allocated()`
guard turns it into `qp_bug` with a line naming the cache, which is the one
failure mode worth paying for here.

It was an EZFIO keyword (`mo_cholesky_pair_cache_mb`) in the first version of
this section, and the generated interface is what killed it: the provider that
`scripts/ezfio_interface/ei_handler.py` writes for a key of interface `provider`
does not fall back to `default:` when the key is absent. It prints
`mo_two_e_ints/<key> not found in EZFIO file` and `stop 1` - see the
`mo_cholesky_double` provider in the generated
`src/mo_two_e_ints/ezfio_interface.irp.f`, which is the template. A new key
would therefore refuse every EZFIO directory written before it, including the
3.38 GB TZ fixture, until it is rewritten - and a cache that sizes itself needs
no key at all. The override remains what the A/B wants: the two arms differ in
the cap and in nothing else, in the same directory, with the same binary.
`QP_PAIRBLOCK_STATS` (Stage 0) and `QP_NTHREADS_PT2` are the same pattern. The
startup line is the record of what a given run used; nothing else is needed.

Providers in `map_integrals.irp.f`, in **four** blocks - the cap, then the
scalars, then each array, so that no dimension is read before the block that
writes it has run:

```
BEGIN_PROVIDER [ integer, pair_cache_mb ]
BEGIN_PROVIDER [ integer, n_pair_cache ]                            ! 0 = disabled
&BEGIN_PROVIDER [ integer, n_pair_cache_padded ]                    ! max(1,n_pair_cache)
&BEGIN_PROVIDER [ integer, pair_cache_nhash ]                       ! 2*n_pair_cache_padded+1
&BEGIN_PROVIDER [ integer, pair_cache_used ]                        ! 0..n_pair_cache
&BEGIN_PROVIDER [ integer, pair_cache_check, (2) ]                  ! 1 = block check done, 2 = column
BEGIN_PROVIDER [ double precision, pair_cache, (mo_num, mo_num, n_pair_cache_padded) ]
BEGIN_PROVIDER [ integer, pair_cache_key, (pair_cache_nhash) ]      ! 0 empty, -1 filling, >0 ready
&BEGIN_PROVIDER [ integer, pair_cache_slot, (pair_cache_nhash) ]    ! 1..n_pair_cache
```

The split is not cosmetic. IRPF90 emits the allocations of a provider block
*before* that block's body runs, so a dimension written in the same block as the
array it sizes is read while it is still 0; the scalars must live in a block the
array block depends on. And only the first line of a block carries no `&`: a
leading `&` is a continuation marker, so a block whose first line starts with one
has no opener at all - `irpf90` then fails with `Beginning of block not matched`
on the `END_PROVIDER` that closes it.

Build block: `n_pair_cache = 0` unless `do_mo_cholesky` **and**
`mo_cholesky_double` and `pair_cache_mb > 0`, else
`n_pair_cache = min(pair_cache_mb*131072/(mo_num*mo_num), mo_num*mo_num)` slots -
the cap in doubles divided by one block, clamped to the number of keys
`k + mo_num*(l-1)` that can exist. Past that no slot can be written, and the
only thing an excess cap buys is the zeroed `2*n+1`-integer bucket table, which
*is* real memory at startup: without the clamp a deduced cap gives a
30-orbital system gigabytes of it. The slot bound is computed *before* the
array's dimension expression, which is the whole point of the separate provider
blocks: the cap is enforced by arithmetic, and the `allocated()` test in the
`pair_cache` body turns IRPF90's silent `stat=` failure into a `qp_bug` instead
of a segfault three minutes into the walk. IRPF90 emits a plain
`allocate(...,stat=irp_err)` and on failure only prints `Allocation failed:
pair_cache` and carries on with an unallocated array
(`external/irpf90/src/variable.py:581-597`). So the cap must be chosen to fit,
the `allocated()` test is the only thing between a failed allocation and a
segfault, and the startup line below says what the cap bought. The
dimension is `max(1, n_pair_cache)` because IRPF90 guards every allocation with
`dim > 0`.

`mo_cholesky_double` belongs in that condition because the fp32 mode was never
measured: there the block path sums each element with `sgemm` and the column path
with `sgemv`, so a column hit would return the `sgemm` bits where today's walk
returns the `sgemv` ones. The gemm-against-gemv difference was measured in double
only (2105/2112 elements, ~1e-11 relative), and nothing says the fp32 sums differ
by the same amount - `sgemm` kernels are frequently blocked differently from their
dp counterparts. Putting an unmeasured rounding difference into a mode that is
dead by default is not a trade worth making, so the cache is not armed there. A
run that set `QP_PAIR_CACHE_MB` anyway says why; the deduced cap stays silent,
since the integral mode is a choice and not a misconfiguration:

```
Pair cache: cap 2048 MiB (QP_PAIR_CACHE_MB)
Pair cache: cap 2048 MiB ignored: needs do_mo_cholesky and mo_cholesky_double
Pair cache: QP_PAIR_CACHE_MB="2G" is not an integer; cache disabled
```

The third line exists because `read` of a malformed value returns a nonzero
`iostat`, and a silent fall back to 0 would turn a typo into the 8.6x and leave
nothing in the log to find it by.

Allocation is not residency: the emitted `allocate` does not zero the array, so
its pages stay untouched until a fill writes them, and RSS follows the blocks
actually cached rather than the cap. Do not lean on that - with over-commit an
over-large cap succeeds silently and surfaces later as swap.

### Why the cache fills lazily, with no selection file and no pre-walk pass

The first version of this section fitted a hot set offline (`pair_curve.py
--emit-select`), read it back at startup and filled it with a `!$OMP PARALLEL DO`
before the walk. The 264-MO run retired that machinery - not because it was
wrong about *what* to fill, but because the walk already knows that; what it
cannot say is *which to keep* if the cap binds (last bullet).

| probe | MOs | rank | block calls / distinct pairs | column calls / distinct pairs | integrals of walk |
|---|---|---|---|---|---|
| cc-pVDZ, 46 880 det, complete, 1 thread | 114 | 2694 | 1 447 862 / 8100 = **179x** | 139 365 752 / 8144 | 80 % |
| cc-pVTZ, 1428 det, complete, 4 threads | 264 | 6064 | 235 888 / 1833 = **129x** | 4 684 818 / 1844 = **2541x** | 84 % |

* **the block path recomputes the same GEMM 129 to 179 times per walk.** A block
  is a `dgemm` of `mo_num x mo_num x rank`, identical on every call for a given
  `(k, l)`. At 264 MOs it measures 33.9 ms of thread time and the walk spends
  8008 thread-s issuing it 235 888 times to produce 1833 distinct matrices:
  **99.2 % of the block time is recomputation** - 62 thread-s of fill against
  8008. The redundancy, not the window, is what Stage 1 removes; at 114 MOs the
  same accounting gives 1377 of 1385 block thread-s wasted.
* **the split is 57:43 and it does not drift.** 47.4 % of the four-thread budget
  went to blocks, 36.2 % to columns, 16.4 % to everything else (4226 s x 4 =
  16 904 thread-s). The same walk cut at 60 % of its
  length reports the same 57:43 split of the integral time. The two timers are disjoint - `get_mo_two_e_integrals_ij`'s wraps only
  a `dgemm` in the Cholesky branch, and the column calls it appears to contain
  belong to the map-mode `else` branches, unreachable when `do_mo_cholesky`
  (`map_integrals.irp.f:687-770`) - so the percentages add up rather than
  overlap. They are wall-clock windows, though, so on a loaded box they also
  contain the time a thread spent off-CPU: they rank the work, they do not
  partition its CPU time (measured in *Stage 1 gate at 264 MOs*, below).
* **the whole pair set is settled in the first 3 % of the walk, but arrival order
  is not a ranking.** At 264 MOs the snapshot holds 1656 of the final 1844 pairs
  after 45 s, 1804 after 90 s and all 1844 after 135 s, against a 4226 s walk:
  the mask structure touches nearly the entire set before the walk has revealed
  anything about multiplicity. Membership is therefore free, *ordering* is not.
  Measured on the final table, the top 155 pairs carry 50 % of the value, the top
  530 80 % and the top 1073 95 % - but a cap that keeps the first 530 to arrive
  keeps the same *number* of pairs with no reason for them to be the valuable
  ones. This is the one place the offline ranking still earns its keep: if the
  cap holds the whole set (975 MiB here, 803 MiB at 114 MOs) the question
  disappears, which is the case this spec is written for. If a production cap
  binds, rank offline - no online admission filter can recover the ordering, be-
  cause the pairs are all seen within the first 1 % of the walk, long before
  their multiplicities diverge.
* **a miss is rare: 0.8 % of block calls, 0.03 % of column calls.** One lazy fill
  per pair, paid by the first request that needs it, costs exactly the work that
  request was going to do anyway. Nothing to schedule, nothing to warm up,
  nothing to go stale when the determinant list changes. 11 pairs (1381 of the
  4.68M column calls) never see a block call at all, which is why they are never
  cached.
* **the cache is small**: 1833 x 0.53 MiB = 975 MiB at 264 MOs, 8100 x 0.10 MiB =
  803 MiB at 114. Both fit the 8 GiB laptop with the fixture in RAM (the TZ
  fixture peaks at 4.2 GB, so 5.2 GB with the cache).
* **what that buys at 264 MOs**: the walk's 16 904 thread-s split into 8008 block
  + 6119 column + 2777 other. Stage 1 keeps 62 s of fill (priced at the contended
  33.9 ms, i.e. conservative - filled one or two at a time it measures 9.75 ms),
  ~12 s of block copies, ~16 s of column reads and *all* of the other, so the
  walk becomes ~2870 thread-s = **5.9x** (wall ~715 s), against 4.8x measured at
  114 MOs. While the cap holds the set, the win grows with `mo_num` because the
  block path scales as `mo_num^2 * rank` and the column path as `rank * mo_num`;
  at 114 MOs the block path is 16 % of the walk, at 264 it is 47 %. The gate
  measured **8x**, not 5.9x: the 16.4 % residual was a thread-s figure taken on a
  contended machine, while the walk's actual bookkeeping is ~1300 s of CPU, and
  one number closes both the cap-0 and the capped arm (*Stage 1 gate at 264
  MOs*, below). The procedure that conserves the residual is the one validated
  at 114 MOs, where the residual *is* real bookkeeping and the projection landed
  at 1747 s against 1748 s measured.
* **the column path is at the bandwidth wall**: 1.386 ms of thread time to
  stream 9.7 MB per call is 7.0 GB/s per thread and 28 GB/s for four - the
  socket's DRAM bandwidth. Threads cannot reduce it (the one-thread calibration
  `t = 0.276 us + 18.3 us/MB` predicts 177 us for the same 9.7 MB), and a
  narrower window makes every miss read a bigger slab. So the column side is not
  short of arithmetic or of threads; it is short of a copy of the data, which is
  what a block cache is.

So Stage 1 has one knob, the cap, and no artifact - and no reason to exclude the
CIPSI driver, where a disabled cache (`QP_PAIR_CACHE_MB=0`) leaves today's
behaviour untouched.

Filling, on miss, with single-flight:

The table is open-addressing over `pair_cache_nhash = 2*n_pair_cache+1` buckets,
linear probing, `slot = mod(key, pair_cache_nhash)+1` with
`key = k + mo_num*(l-1)`. `pair_cache_key` carries three states: `0` empty, `-1`
being filled by some thread, `> 0` ready (an empty bucket ends a probe). A block
request is:

1. **lock-free probe** - plain reads along the chain. `> 0`: `!$OMP
   FLUSH(pair_cache_key)`, copy the block out, return. `-1`: call the uncached
   routine and return - **never wait**. This is the path 99.2 % of the calls
   take at 264 MOs, where 235 888 block calls cover 1 833 pairs;
2. on a miss, `!$OMP CRITICAL(pair_cache_probe)` and re-probe, because the state
   may have changed in between. Ready: leave, copy out, return. `-1`: leave,
   uncached, return. **Slot pool exhausted** (`pair_cache_used >= n_pair_cache`):
   leave, uncached, return - this is the freeze. Otherwise take the next slot
   (`pair_cache_used = pair_cache_used + 1`), write `pair_cache_slot(bucket) =
   slot` and `pair_cache_key(bucket) = -1`, leave;
3. compute the block through the uncached routine into the caller's `out_array` -
   outside the critical, so a 34 ms GEMM does not serialize the other walkers;
4. copy `out_array` into the claimed slot, `!$OMP FLUSH`, then write `key`. Only
   the claimer writes that slot, so no critical and no atomics are needed on the
   publish; the flush is what orders the data before the key.

So the walk pays one critical section per *pair*, 1844 of them at 264 MOs, and
none per hit: a hit is a probe, one flush and a copy. A thread that loses the
race for the same block computes its own copy, exactly as today, and does not
store it - the cost of the race is bounded by today's cost, and the bits are the
same either way, because both threads run the same routine.

**Only the block routine fills.** The column routine probes and, on a miss,
takes today's streamed path; it never writes to the cache. A column call has no
block to store, and the pairs that only ever see column calls carry 0.01 % of
the column calls at 114 MOs (16 462 of 139.4M, 44 pairs) and 0.03 % at 264
(1 381 of 4.68M, 11 pairs), so the case does not deserve a second fill path.

One hazard remains, and it is the only one: a reader could see `key > 0` and
stale data. The flush pair closes it - the writer's flush after the copy and
before the key, the reader's flush after the key and before the copy - and a
reader that sees `-1` or `0` takes the uncached path, which is always correct.
The reader's flush is a `dmb ish` on arm64, tens of cycles against a 0.3-3 us
hit.

When the slot pool is empty the table freezes: later pairs are computed uncached
and the run stays correct, only slower. That is the default, and the only policy
the four steps above implement. It is the right one whenever the cap holds the
whole pair set, which is true of both completed benches; where the cap binds it
is wrong, because arrival order is not a ranking, and *When the cap binds* below
gives the replacement that fixes it together with its price.

The freeze test is the slot pool, not the bucket table, and the difference
matters. With `nhash = 2n+1` buckets and at most `n` published keys, at least
`n+1` buckets stay empty for the life of the process, so a probe always
terminates at one: "no empty bucket on the path" is unreachable and would be dead
code. What does happen is an empty bucket on the path while all `n` slots are
already claimed, and that is the case the check must catch, or the first claim
past the pool writes `pair_cache(:,:,n_pair_cache+1)` out of bounds.
`pair_cache_used` is a provider so that the arming `PROVIDE` covers it and it is
zeroed with the table.

The fill passes a locally declared `type(map_type) :: dummy_map` to the uncached
block routine. Naming `mo_integrals_map` there would ask IRPF90 to build the
sparse MO map, which Cholesky mode never builds and does not need; the map
argument is only read in the `do_mo_cholesky = False` branch, and the provider
already requires `do_mo_cholesky`.

Stage 0 writes the full table, so the cap can be chosen after the measuring run.
Nothing in Stage 1 depends on the snapshot format.

Where the cache is armed, and by whom:

* the block and column routines run inside the walk's `!$OMP PARALLEL`, so the
  provider is armed before it, twice for good measure: `PROVIDE pair_cache` in
  `pt2`'s pre-loop block (`fci/pt2.irp.f`, next to
  `PROVIDE mo_two_e_integrals_jj`, `:59`) and in `pt2_omp_task`'s
  (`cipsi_utils/pt2_omp.irp.f:272`, next to `PROVIDE banned_excitation`). With
  the cache disabled (`QP_PAIR_CACHE_MB=0`) the body allocates one slot, prints
  nothing, and returns; with the deduced cap it prints its line and allocates
  the pool;
* that is about *when*, not *whether*. IRPF90 already inserts a guarded
  `provide_pair_cache` at the top of every routine that references any of the
  entities - the generated `get_mo_two_e_integrals_ij` opens with
  `if (.not.irp_is_built(512)) call provide_mo_num`, and the same for each
  dependency - so a caller cannot read an uninitialised `n_pair_cache` whatever
  path reached it. This matters because these routines have callers outside PT2
  (`mp2`, `davidson`, `mu_of_r`, `casscf_cipsi`, `gamma_of_r`, `selection`).
  The generated guard is a nest lock from `irp_locks`, initialised in every
  `main`, so four threads reaching it at once collapse to one build; the two
  `PROVIDE` statements above only make sure that build happens *before* the
  parallel region, in one thread, instead of inside it.

`pt2` is a standalone OpenMP program: the cache is process-local, shared by the
walk threads, and pays no MPI broadcast.

The provider prints one line when it is enabled, and nothing when the cap is 0 -
the same contract as Stage 0's report, and the only way to tell afterwards that
the run you are looking at actually had the cache it claims:

```
Pair cache: cap 2048 MiB -> 3851 blocks of 557568 B, 2047 MiB allocated
```

What the cap bought, in blocks, with the block size in bytes so the arithmetic
can be checked against `mo_num` (`8*mo_num^2`); nothing is printed when the cap
is 0. Hits and misses are Stage 0's own counters, which already separate the two
classes per `(k, l)`, so a Stage 1 run needs no new instrumentation - but their
meaning changes, and it has to be read correctly:

* the Stage 0 timers sit in the `*_uncached` bodies, so with the cache armed
  `block_thread_s` is *fill* time and `block_requests` is the *miss* count: 1833
  distinct pairs plus a handful of races, against 235 888 block calls in the
  cap-0 walk of the same fixture;
* the hit rate is therefore `1 - block_requests(cap)/block_requests(0)`, both
  from completed runs of the same walk, and the same for the column class;
* a hit is not timed at all. That is deliberate: what Stage 1 is judged on is the
  wall time of the completed walk, and two `clock_gettime` calls on a 0.3 us hit
  path would measure the instrumentation rather than the cache.

Out of scope by design: the iterative CIPSI driver needs no change. With the cap
at 0 it is untouched; with a cap it fills lazily, per process, and a growing
determinant list is harmless - a block is a function of `(k, l)`, not of the
determinants, so blocks cached in one iteration stay valid in the next, and the
cap simply stops admitting new ones. The production-size measurement remains a
`pt2` run on a production wavefunction, which is the same program.

Split both routines:

* `get_mo_two_e_integrals_ij_uncached(k,l,sze,out_array,map)` = today's body,
  unchanged, all branches (GPU, dgemm, sgemm, map).
* `get_mo_two_e_integrals_ij(k,l,sze,out_array,map)` = `PROVIDE pair_cache`,
  then `if (n_pair_cache == 0 .or. sze /= mo_num)` -> uncached; else probe; hit ->
  `out_array(1:sze,1:sze) = pair_cache(1:sze,1:sze,slot)`; miss ->
  the single-flight fill of *Filling* above, which ends by returning the block it
  just stored.
  The `banned_excitation` test is not in this routine; it stays where it is.
* `get_mo_two_e_integrals_uncached(j,k,l,sze,out_val,map)` = today's body of the
  column routine, unchanged, with its `banned_excitation(j,l)` early return
  first (it returns zeros and needs no block).
* `get_mo_two_e_integrals(j,k,l,sze,out_val,map)` = `PROVIDE pair_cache`, the
  same `n_pair_cache == 0 .or. sze /= mo_num` short circuit, then the
  `banned_excitation(j,l)` test **before** the probe, returning zeros exactly as
  today; then probe; hit -> `out_val(1:sze) = pair_cache(1:sze,j,slot)` - one
  contiguous read of `mo_num` doubles, 2.1 KB at 264 MOs, because `i` is the
  fastest index of `pair_cache` exactly as it is of `out_array`; miss -> the
  uncached routine.

The probe is one shared helper, `pair_cache_probe(k,l,slot)`, returning the slot
or 0. Both wrappers call it, so the open-addressing walk exists once.

Every wrapper arms the cache with its own `PROVIDE pair_cache` rather than
trusting a caller. IRPF90 does not initialise a provider variable, so reading
`n_pair_cache` without building it reads whatever was in memory - and these are
not PT2-private routines: `get_mo_two_e_integrals_ij` is called from `mp2`,
`davidson`, `mu_of_r` and `casscf_cipsi`, and the column form from `mp2`,
`gamma_of_r` and `selection`. Providing the same provider from four threads is
safe: the generated guard is a nest lock from `irp_locks`, whose 740 entries are
initialised by `irp_init_locks_*()` at line 10 of every program's `main`, so the
second and later `PROVIDE`s are a no-op after the first.

The banned test is duplicated on purpose. It must run before the probe, because
a banned `(j,l)` returns zeros and *not* column `j` of the block: a cache hit
that skipped it would return the true integral. Keeping it in the uncached body
too is what lets that body stay a verbatim copy of today's routine; the test is
two array loads and a branch.

Both routines use the same `slot` table, so a block computed for `get_d0` is a
hit for the `get_d1` columns that follow it in the same mask, and vice versa.

Properties:

* The block class is bit-identical by construction in both precision modes: the
  block is always produced by the one existing routine, with the same operands
  and shape. No second implementation of the integral computation. A fill dgemm
  is issued from a walk thread, inside the parallel region, with the BLAS thread
  count pinned to 1 by `qp_blas_env_init`, exactly like the walk dgemm it
  replaces - and that equivalence is measured, not assumed: arm A of the closed
  A/B issued 1.45M block dgemms from four threads, arm B read the same matrices
  from an array the fill had built, and the two printed rows are identical to
  all 17 digits. So the bits do not depend on which thread runs the dgemm or on
  when it runs, and a block hit returns exactly what a block miss returns.

* The **column** class is not bit-identical, and that was measured before the
  fill was written. The walk gets `<ij|kl>` through two entry points that do not
  agree in the last bits: the block routine computes the whole matrix with one
  `dgemm` (`map_integrals.irp.f:951`), the column routine computes its column
  with a `dgemv` (`:640`). A 30-line program on this machine's BLAS (Accelerate,
  what `-llapack -lblas` resolves to in this tree), at the fixture's real shape
  `chol = 2694`, `n = 264`:

  ```
  gemm vs gemv:             2105 / 2112 elements differ in bits, max rel 1.8e-11
  gemv full vs gemv window:    0 / 512  elements differ in bits
  ```

  A column hit therefore returns the gemm bits of a value that is computed by a
  gemv today. The difference is ~1e-11 relative and not 1 ulp, because the sums
  have 2694 terms and cancel. (`~/OpenCode/tmp/prof/gemm_vs_gemv.f90`; the
  fixture's real Cholesky vectors may cancel harder still, which is why the
  check below is norm-based and not element-wise.)

  The second line is the control, and it corrects the record: the shift-6 vs
  shift-7 A/B compared one dgemv against three shorter dgemvs, which is
  bit-invariant, so that A/B was a *null* test for this question and its
  byte-identical row says nothing about gemm against gemv. An earlier version of
  this section proposed it for exactly that purpose; it cannot answer it,
  because at shift 7 the column routine still calls dgemv.

  Three consequences, all accepted deliberately:

  - **The instrument is restated.** Block hits are checked bit-exactly against
    `get_mo_two_e_integrals_ij_uncached` (dgemm against dgemm, zero tolerance -
    this is what catches a transposed, shifted or mis-strided block). Column
    hits are checked against `get_mo_two_e_integrals_uncached` on the 2-norm of
    the column, `||d|| / ||y|| <= 1e-9`, and the measured ratio is printed once
    per process. A layout error gives O(1) - nine orders of margin; a rounding
    difference gives ~1e-11 and passes.
  - **The printed `pt2` row is not bit-identical between cap 0 and cap M.** It
    agrees to ~1e-10 relative, and the deviations are reported with the row.
    That is what the arithmetic allows on 17 printed digits; the Cholesky
    threshold of these runs is 1e-6 to 1e-8, so no physical quantity lives in
    the difference.
  - **Determinism.** In the uncapped regime the result still does not depend on
    the thread count: every pair is eventually filled, and a pair's bits do not
    depend on who filled it or when. In the *capped* regime, where the freeze
    cuts different pairs in different runs, which pairs are hits varies with
    arrival order and the row varies at 1e-11. That caveat belongs to the capped
    regime, it does not exist when the cap holds the whole pair set - which the
    deduced cap does for every bench whose demand is under half the budget -
    and the replacement policy below only ships if it fixes admission
    deterministically.
* By default the array is append-only: a slot is written once, by the thread that
  claimed it, and never moved or freed. No eviction, no resize. The bounded mode
  below is the one exception and changes exactly this sentence.
* On a block hit the cost is one `8*mo_num^2` byte copy: 0.53 MiB at 264 MOs,
  ~50 us, **0.15 %** of the 33.9 ms dgemm it replaces. On a column hit it is one
  contiguous read of `mo_num` doubles - 912 B at 114 MOs, 2.1 KB at 264 - for the
  same reason the block is contiguous: `pair_cache(i,j,slot)` has `i` as the
  fastest index, so column `j` is a run of `mo_num` doubles, not a stride. The
  cost is therefore the *measured* 0.293 us of the shift-7 copy path, not a
  bound: that path does exactly this copy (912-byte contiguous read plus the
  probe) 139 365 752 times in 40.9 s. At 264 MOs the read is 2.3x wider and a hit
  is 0.3-0.4 us against the 1.31 ms dgemv it replaces, **>99.9 %** of the call;
  at 114 MOs, 0.29 us against 38 us, 99 %.
* A miss is exactly today's cost plus a hash probe and two critical sections. A
  hit is a probe, one `FLUSH`, one branch and a copy: no locks, no atomics, no
  per-thread state.

Not in scope: a traced key sequence, a runtime-adaptive policy. Replacement is in
scope only as the bounded mode below, with its price stated there. If the cap
turns out to be too small for the production pair set, the first answer is that
mode; the second is the batched device GEMM discussed under the decision rule.

### Implementation order, and the gate of each step

Three commits, each with a gate that is cheap enough to run before the next one
starts. None of it needs a cluster, and every gate is a completed TZ walk at
`shift = 6`, 4 threads, on `tz_basis_test.ezfio` (`tz_fixture_notes.md`), because
the aggregate row is the only artifact that can see an O(1) index error and the
per-value assertion below is the only one that can see a 1-ulp column error.

1. **The `get_d0` workspace hoist** - ten lines, the table in *The one hoist that
   does pay*. It comes first and alone so that the cache diff is pure. Gate:
   `ninja` clean, then the cap-0 row is the archived one to all digits
   (`PT2 = -0.86579718276647843E+00`, `rPT2 = -0.67308215287443873E+00`,
   `E+PT2 = -231.78118240293102`, `E+rPT2 = -231.58846737303898`,
   `E = -230.91538522016455`, `PT norm = 0.53508620642076787`,
   `overlap = 0.28631724830176858`) and `block_thread_s/block_requests` is
   reported, not gated.
2. **The cap, the providers and the split, cache still at 0.** No `EZFIO.cfg`
   edit and no interface regeneration: the cap is `QP_PAIR_CACHE_MB`, read by
   `pair_cache_mb`. The four provider blocks of `map_integrals.irp.f`; the
   `*_uncached` bodies copied verbatim; each wrapper starts with
   `if (n_pair_cache == 0 .or. sze /= mo_num)` -> uncached, so a cap-0 run pays
   one integer test and one comparison per call and nothing else, and the
   `PROVIDE pair_cache` in `fci/pt2.irp.f` and `cipsi_utils/pt2_omp.irp.f`
   allocates the single padding slot before the walk starts. Gate: `ninja` clean,
   the wrappers read against the original bodies line by line, and the cap-0 row
   unchanged - which the step-3 gate re-checks on the final binary, so this one
   can be a diff review plus a short run to first output.
3. **The fill**: probe, claim, compute, publish, the startup line, and the two
   one-shot assertions of the acceptance list. Gate: cap 2048 MiB, 4 threads,
   every block hit bit-exact and every column hit within 1e-9 on the column
   norm (both asserted in-process, the measured column ratio printed), the eight
   row numbers within 1e-9 relative of step 1's, and the snapshot's
   `block_requests` against the 1833 distinct pairs giving the hit rate. Bit
   equality of the row is *not* the gate and cannot be - see the column bullet
   in *Properties*. The wall time is the headline (target ~715 s against the
   4279 s cap-0 walk) and is reported, not gated: a contended laptop has already
   produced a 4 % swing between two runs of an identical arm.

**Outcome.** Steps 1-3 are implemented in that order in the branch
`pt2-pair-block`, and the gate above ran on the final binary: step 1's cap-0 row
is the archived one (*Reference row, 264 MOs*), step 3's 2048 MiB arm gives zero
block misses and one column call with both assertions silent, and the 512 MiB arm
was added because a cap above the union never exercises the freeze branch or the
miss counters. The wall target of ~715 s was met with room: 428.4 s on an idle
box, 8.63x against the cap-0 arm of the same day. *Stage 1 gate at 264 MOs* has
the record.

The assertions are once per *process* per class, not per thread, and the hot path
pays one plain integer read for them. The wrapper tests the flag lock-free -
after the first hit it is 1 and the branch is never taken - and only a thread
that sees 0 enters the verifier, which re-tests and sets the flag under the same
`CRITICAL(pair_cache_probe)` that already exists, keeping the claim atomic so a
walk that starts with four threads firing at once does not run four identical
34 ms dgemms and print the same line four times. The comparison itself runs
outside the critical into a local `allocatable` scratch (`transfer()` to
`integer(8)`), then `qp_bug` on any mismatch; the flag is set before the
comparison, so a failure is reported by the thread that claimed it and the
others stay out of the way. The block hit is checked against
`get_mo_two_e_integrals_ij_uncached` (one 34 ms dgemm) and the column hit against
`get_mo_two_e_integrals_uncached` (one 1.3 ms dgemv). They are always compiled in
and never gated on an environment variable: a branch on a zero flag per hit is
free, and the check is the per-value guarantee at production shape.

### When the cap binds: replace instead of freeze

Freeze is wrong at production size, and the numbers above say why: the cap holds
1031 blocks at 510 MOs against a pair set already 1833 at 264, and the pairs the
walk goes on to ask for again are not the ones that arrived first - the two are
uncorrelated in the measurement. Freezing hands the first 1031 arrivals a
permanent monopoly on the table.

In-place replacement is not free here, and this table makes it cost two things:

* **probe integrity.** Linear probing looks a key up along `mod(key,nhash)+1`,
  `+2`, ... and stops at `key == 0`, so every key must sit on its own probe path.
  A replacement may only claim a bucket that lies on the *new* key's path -
  possible only when that path holds no empty bucket, which is exactly the case
  that reaches replacement - so the new key is still found. The evicted key then
  vanishes from its own path, which reads as a miss by construction, because its
  probe walks past the foreign key to the `0` terminator. Replacing at a hand
  pointer instead - classic clock - would strand the new key off its own path and
  make it unfindable.
* **the reader may be copying the block being overwritten.** The `-1` state is
  written before the data is touched, so a reader that has not yet passed the
  probe takes the uncached path; a reader already copying cannot. It must
  therefore validate: re-read the key after the copy and discard the copy if it
  changed. That is sound except for ABA - the evicted key re-admitted into the
  same bucket during the ~50 us copy - which needs a sequence number:
  `pair_cache_stamp(nhash)`, odd while a slot is being rewritten, read before and
  after the copy. Because a block is bit-identical every time it is computed, the
  key alone would do if keys were never re-admitted, and the stamp is what buys
  the right to re-admit.

Cost: one `integer` array of `nhash` (31 KB at 264 MOs with the 2 GB cap), two
loads and two `FLUSH`es on the hit path, and a rare lost hit, which recomputes at
today's cost.
What it buys is a *thrash bound*: two colliding pairs take turns instead of one
starving for the rest of the walk. It does **not** buy the value ranking - with
1031 slots for 1833+ pairs the concentrated pairs still have to be chosen, not
discovered - so if the production dump comes back with the 264-MO shape (top 530
pairs = 80 % of the value) then the ranking is what pays, and this mode is not a
substitute for it.

Implement it last, behind the same keyword, and only after the cap-0 arm and the
append-only arm have both been timed on the same binary. Until then the case for
replacement is an argument, not a number.

## Stage 0 result (benzene cc-pVDZ probe, 20 November 2026)

Four runs of the forced-Cholesky fixture (`benzene_chol.ezfio`: 114 MOs, rank
2694, window 1-64), same binary, `QP_PAIRBLOCK_STATS=1`,
`QP_PAIRBLOCK_DUMP=30`, each stopped at a fixed time budget rather than at the
end of the walk (the walk does finish: 2 h 09 at shift 6 and 26 min at shift 7,
one thread, measured - see Validation). The 8-thread runs maximise the number of walk
positions covered; the 1-thread runs give uncontended per-call costs.

| wall | thr | block calls | block share | column calls | column share | pairs | both |
|---|---|---|---|---|---|---|---|
| 60 s  | 8 | 11 222  | 11.1 % | 1 635 761  | 81.7 % | 7 979 | 3 414 |
| 570 s | 8 | 79 149  | 10.2 % | 13 667 186 | 84.1 % | 8 144 | 6 132 |
| 480 s | 1 | 46 679  |  9.0 % | 10 973 097 | 83.2 % | 8 144 | 5 470 |
| 870 s | 1 | 90 299  | 10.5 % | 17 944 922 | 78.5 % | 8 144 | 6 412 |
| **7674 s** | **1** | **1 447 862** | **16.6 %** | **139 365 752** | **62.3 %** | **8 144** | **8 100** |

The last three rows ran alongside other single-thread probes for their final
minutes, which is why their column share is lower and "everything else" higher;
their counts are the deeper ones and the per-call costs match the 480 s row.

The bold row is the **completed walk** (`full_s6_T.log`, flags `QP_PAIRBLOCK_STATS=1`,
no dump), the only unbiased row in the table: the others are 60-870 s of a 7674 s
walk, so they sample the head, where `pt2_J` puts the generators that issue the
most calls. Over the whole walk the two contractions are **78.9 %** (16.6 + 62.3),
not the 92-94 % the head shows, and the per-call means of that row (0.8796 ms
block, 34.28 us column) are what the ceiling arithmetic below uses.

Per-call cost, one thread: `t_blk` 0.93-1.55 ms (70 MFlop, 45-75 GFlop/s),
`t_col` 36 us (1.08 MB, 30 GB/s), ratio **21-26**. With eight threads the block
inflates 4x and the column 8x (bandwidth), so the ratio stays near 21: the
shares are robust, the absolute times are not.

What the numbers say:

* the `sample` profile is confirmed by an independent instrument: 9-11 % block,
  82-84 % column in the head, i.e. 92 % of the *head* in the two contractions -
  and 78.9 % of the whole walk, which is the number that bounds any change;
* the pair set **saturates**: 7979 pairs after 60 s, 8144 after 480 s, and it
  does not grow afterwards while the call counts keep climbing. 8144 is 63 % of
  the 12 996 ordered pairs of 114 orbitals. Coverage therefore *converges* with
  walk time instead of decaying - the opposite of the staleness worry, at least
  within one wavefunction;
* 60-80 % of column calls land on a pair whose block the walk computes anyway
  (69 % at 60 s, 80 % at 570 s, 83 % at 870 s), for 3.4-6.4 k blocks. The
  completed walk closes this: 139 349 290 of its 139 365 752 column calls land on
  a pair the walk fills a block for anyway - **all but 2462 calls, on 44 of the
  8144 pairs**. There is nothing left for a smarter admission policy to find at
  this size: the set is not merely saturating, it is closed;
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
and the completed walk measures that remainder directly: at shift 7 the whole
walk takes 1478.9 s with 40.9 s of it in the column path and no block call at
all. Same fixture, same binary, same call sequence (verified pair by pair): the
residual is not-integral work, and it is what `T_walk - T_int` means below.

The arithmetic that ties the two measurements together, one thread, shift 6:

```
T_walk = C_total x t_col + B_total x t_blk + T_other
       = 139 365 752 x 34.28 us + 1 447 862 x 0.8796 ms + T_other
       = 4777 s               + 1274 s                + T_other
T_other, measured on the completed shift-7 walk       = 1438 s   (1478.9 - 40.9)
T_walk, measured                                      = 7673.7 s (residual 1623 s)
T_int/T_walk = 6051/7674 ~ 0.79
```

The residual is 185 s larger in arm A than in arm B. The streamed column evicts
the walk's own working set - 1.08 MB through the caches per call - and that costs
13 % on the work no counter attributes. So the integral share is 79 % of the
walk and the ceiling for any integral-side change is 80 %, above the 75-78 % the
two extrapolations bracketed.

`C_total` and `B_total` are counted rather than inferred now (139 365 752 and
1 447 862), and `t_col`/`t_blk` are the class means the same run reports. The
head's own rate gave a lower `T_walk`, ~6300 s, and the gap is `T_other` being
mostly per walk *position*: a head rate weighs the call-rich positions more than
the whole walk does. The curve conversion above uses 0.8 and is insensitive to
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
that settles it. Nothing here justifies a larger cap than ~2 GB per process.
Eviction never pays while the set fits; its only justification is the capped
regime, which is exactly the regime this paragraph cannot see.

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
default must stay 1. `pt2_J` puts the call-rich generators first, so the
**pair set and its multiplicities saturate in the first 3 % of the walk** and a
two-minute snapshot is enough for coverage. The **split is the opposite** and the
264-MO run measures how badly it drifts:

| fraction of walk | 1 % | 2 % | 4 % | 10 % | 31 % | 36 % | 47 % | 63 % | 100 % |
|---|---|---|---|---|---|---|---|---|---|
| block share of integral time | 18.5 % | 14.5 % | 24.1 % | 31.7 % | 26.2 % | 57.8 % | 57.2 % | 57.0 % | 56.7 % |
| integrals / walk | 23.5 % | 23.6 % | 41.2 % | 68.8 % | 48.6 % | 87.5 % | 86.4 % | 85.2 % | 83.6 % |

It converges at ~36 % of the walk and then holds to within one point for the
remaining 64 %, but the head says 18 % block when the answer is 57 %. So: kill
the job early for the pair set, **do not kill it early for the split** - and never
for a per-call cost, which is contention-dependent by a factor of 3.5 (the block
dgemm measures 9.75 ms/call at 1 % of the walk when one thread at a time is
issuing it, and 34 ms/call once four are). Walk-time prediction is what the
completed A/B is for.

The contention is not noise, it is the regime: a block call streams
`2*mo_num*cholesky_mo_num*8` bytes (25.6 MB at 264 MOs) to do `2*mo_num^2*rank`
flops, so its arithmetic intensity is `mo_num/8` flop/byte - 14 at 114 MOs,
33 at 264, 64 at 510 - against a single-core machine balance near 40. The block
path is bandwidth-bound at 114, at balance at 264, and compute-bound beyond;
that, not the call counts, is what makes its share rise with `mo_num`.

### Both benches are in the uncapped regime, and production is not

The one conclusion the two completed walks cannot support by extrapolation is the
headline. A cap buys `cap*131072/(mo_num*mo_num)` slots, so the slot count falls
as `mo_num^-2` while the pair set grows with the determinant list:

| `mo_num` | bytes per block | slots in a 2048 MiB cap | pair set measured |
|---|---|---|---|
| 114 | 104 KB | 20 600 | 8 100 (2.5x the set, holds) |
| 264 | 545 KB | 3 851 | 1 833 (2.1x the set, holds) |
| 510 | 2.0 MB | 1 031 | unmeasured - and already 1 833 at 264 with only 1 428 determinants |

At the 510-MO fixture the cap therefore binds with certainty, and by the 264-MO
concentration curve the slice that fits carries well under its share of the
value. So **4.8x at 114 and 8.6x at 264 are the ceiling of the uncapped regime,
not rungs on the way to a bigger production number**: the win grows with `mo_num`
only while the whole set fits, and beyond that every extra MO costs a slot and
buys nothing. What the production run must print is the pair set at 510 MOs and
its curve - the two numbers that decide whether the capped regime is worth 2x or
worth nothing. That is why the cap, and not the cache, is the switch: at 264 MOs
the deduced cap holds the whole set and the cache is a free 8.6x, at 510 it is
an open question, and `QP_PAIR_CACHE_MB=0` turns it off while the production run
answers it.

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

## The other no-memory alternative: swap the loops so the block stays hot

The cache keeps *data*. The obvious alternative keeps nothing: turn the walk
inside out so a block is materialised once and every consumer of it runs while it
is still in L1-L3. That is the right shape of idea, and it does not work here, for
a reason that is four lines of `get_d0` (`selection.irp.f:1551`):

```
if(sp == 3) then ! AB
  h1 = p(1,1)
  h2 = p(1,2)
  call get_mo_two_e_integrals_ij(h2,h1,mo_num,hij_cache,mo_integrals_map)
```

`p` is refilled for every `i` of `splash_pq`'s loop over `interesting`
(`:1004-1013`), from `mobMask = not(mask) .and. det(i)`. **The block key is a
property of the interesting determinant, not of the loop that iterates it.** So
there is no loop whose body is "the block is invariant": one call of `splash_pq` -
one (generator, hole pair) unit - issues a fetch (or a set of column reads) per
interesting determinant that reaches `get_d0`, and the 235 888 fetches over 1833
pairs are spread over the units whose interesting lists contain the pair. There is
nothing to hoist.

To make a block hot you must invert onto the pair and run the units that fetch it
back to back - 129 fetches per pair on average. That collides with `mat`, the
unit's carrier: allocated once per generator (`:338`), zeroed at the top of
`splash_pq` (`:973`), accumulated over the **whole** interesting list, and consumed
only at the end by `fill_buffer_double` (`:540-541`). A unit cannot be re-entered
piecewise - the `mo_num^2` sweep needs the complete `mat` - so splitting one
across pair iterations means keeping its partial `mat` alive, and each open unit
also holds `banned`/`bannedOrb` (2 x `mo_num^2` logicals, 136 KB at 264 MOs,
rebuilt per unit from `banned_excitation`, `:506-507`) plus its `interesting`
list.

The two orders cannot both be optimal, which is the whole argument:

* iterate **units**, and the block is refetched per unit - that is today's
  235 888 calls for 1833 blocks;
* iterate **pairs**, and you pay 0.56 MB of partial `mat` (2.1 MB at 510 MOs) per
  open unit. A unit can only close after its *last* pair, so the open set is not
  the consumers of the current block: it is every unit whose first pair precedes
  this point and whose last follows it. At a minimum it holds the units that fetch
  a single pair - 129 fetches per pair on average, so 72 MB at 264 MOs, and more
  for the hot pairs - and in practice, with the pair graph a single connected
  component (58 MOs, 55 % dense at 264 MOs), most of the walk is open at once. So
  the resident set becomes accumulators of the same order as the table it
  replaces, and it grows with hole pairs per generator while the cache does not.

So the menu is: hold the per-consumer state (partial accumulators) or the shared
object (blocks). The cache is the cheaper encoding of the same idea, and it
leaves the control flow, the dedup arrays and the accumulation order untouched,
which is what makes it bit-identical instead of a rewrite of the selection. Its
ceiling is higher too: the loop swap buys the block half (47 %, 1.9x) and does
nothing for the 36 % column half, which the same table serves at no extra cost.

### The one hoist that does pay: the block allocation in `get_d0`

The same read turned up a loose end that has nothing to do with caching. `get_d0`
allocates and frees a full block on every call - `:1545`
`allocate(hij_cache(mo_num,mo_num))`, `:1600` `deallocate` - so at 264 MOs the
walk runs 235 888 malloc/free pairs of 545 KiB, **131 GB of allocator traffic**,
for a buffer that is the same size every time.

The chain is three deep with one call site each, so the buffer has an obvious
owner:

```
select_singles_and_doubles  (:160)   already owns banned/bannedOrb/mat, :337-338
  -> splash_pq              (:953)   called once, :540
       -> get_d0            (:1514)  called once, :1036
```

Ten lines, all in `cipsi/selection.irp.f`:

| line | now | after |
|---|---|---|
| `:197` | `double precision, allocatable :: mat(:,:,:)` | add `double precision, allocatable :: hij_cache(:,:)` |
| `:338` | `allocate(mat(N_states, mo_num, mo_num))` | add `allocate(hij_cache(mo_num, mo_num))` |
| `:540` | `call splash_pq(..., banned, mat, interesting)` | `..., banned, mat, hij_cache, interesting)` |
| `:556` | `deallocate(banned, bannedOrb, mat)` | `deallocate(banned, bannedOrb, mat, hij_cache)` |
| `:953` | `subroutine splash_pq(..., banned, mat, interesting)` | `hij_cache` after `mat` |
| `:966` | `double precision, intent(inout) :: mat(N_states, mo_num, mo_num)` | add `double precision, intent(inout) :: hij_cache(mo_num, mo_num)` |
| `:1036` | `call get_d0(det(1,1,i), ..., banned, mat, mask, ...)` | `..., banned, mat, hij_cache, mask, ...` |
| `:1514` | `subroutine get_d0(..., banned, mat, mask, ...)` | `hij_cache` after `mat` |
| `:1541` | `double precision, allocatable :: hij_cache(:,:)` | `double precision, intent(inout) :: hij_cache(mo_num, mo_num)` |
| `:1545`, `:1600` | `allocate(hij_cache(mo_num,mo_num))`, `deallocate(hij_cache)` | both deleted |

The result is bit-identical, not merely equivalent: both branches of `get_d0`
call `get_mo_two_e_integrals_ij`, whose `out_array` is `intent(out)` and written
whole on both of its branches, and the dgemm runs with `beta = 0`, so the buffer
is overwritten before every read. The allocation sits inside the parallel-do over
tasks (one task, one call of `select_singles_and_doubles`), so the buffer is
private by construction and needs no clause; and there is no early `return`
between `:338` and `:556`, so there is no leak.

Per task, not per thread. A per-thread buffer means threading it through
`omp_selection_task` and `select_connected` - both callers,
`selection_omp.irp.f:325` and `pt2_omp.irp.f:889` - to turn 1 428 allocations
into 40, and up to 80 MB of resident scratch at 510 MOs. Not worth it. Not a
module variable either: a provider-shaped global sized by `mo_num` is the
lazy-provider trap the header of this subtree warns about.

**What it is worth, and where to look.** Less than the call count suggests, and
the faults are inside the timed region, not outside it. The buffer *is* the
dgemm's C matrix, so the first touch happens between `t_block0` and `t_block1`
and lands in **`block_thread_s`**; only the allocator call itself lands in
`other`. The bound is 545 KiB = 34 pages at 16 KiB (arm64) to 136 at 4 KiB, at
~1 us a fault: 34-136 us against the measured 33.95 ms per call, **0.1-0.4 % of
the block path and at most 0.2 % of the walk**. On the cluster at 80 threads the
`munmap` TLB shootdowns across the team are the part that can grow, which is the
reason to do it there. On the laptop: hygiene, not a speedup - it makes a cache
miss exactly "probe + dgemm + copy" with no allocator in the loop, which is what
makes hit-versus-miss a clean comparison. Expect the cap-0 A/B at 264 MOs to land
within noise, and that is the expected outcome: the timings elsewhere in this
file do not need re-measuring for it.

**Two things beside it that must not be touched.** `get_d1` looks the same and is
not: `:1280-1281` allocates 4.2 KB and 0.5 KB, both under any mmap threshold so
both free-list hits, and 4.68 M x ~100 ns is 0.5 s out of 6 119 thread-s.
Cosmetic. The `lbanned = bannedOrb` copy beside it is semantics, not scratch, and
stays. And the fp32 branch is the one per-call allocation that would actually
hurt: `get_mo_two_e_integrals_ij` allocates `real :: out_array_sp(sze,sze)` at
`:722` and copies it back at `:750`, inside the block timer, 272 KiB at 264 MOs -
dead code under the default `mo_cholesky_double = True`. If fp32 is ever used,
hoist that, not `hij_cache`.

## Decision rule

Stage 0 reports, per class (`block`, `column`), the number of requests, the
number of distinct pairs, the thread-seconds and hence the share of the walk;
plus the size of the shared key set. The decision is then bookkeeping:

* A pair with `n_block >= 1`: the fill computes the block once, and from then on
  every block call is a copy, so caching saves `(n_block-1)*t_blk + n_col*t_col`
  for one slot. Never negative, and not small: the completed walks average
  **179 block calls per pair at 114 MOs and 129 at 264**, and at 264 a block is
  33.9 ms against a 1.3 ms column, so one slot pays for itself 26 columns in.
  Today's walk pays the block `n_block` times; the cache pays it once.
* A pair with `n_block = 0` (column-only): caching it means computing a dgemm
  the walk never asked for, which pays back after `ceil(t_blk/(t_col-t_hit))`
  hits - 27 on the probe, where `t_hit = 0.40 us` is now measured rather than
  assumed (Validation), so the earlier `t_blk/t_col` form overstates the
  threshold by 1 %. Where those pairs exist they are usually cheap to justify:
  in the 870 s snapshot 2674 of them average 1647 calls each, far past the
  threshold. The completed walks leave **44** column-only pairs carrying 2462
  calls at 114 MOs and **11** carrying 1381 at 264 (0.03 % of the columns): over
  a whole walk this class is noise, and the ranking below needs it only for
  honesty, not for value.
* The marginal-value curve is therefore: rank the pairs by
  `n_col*t_col + (n_block-1)*t_blk` for `n_block >= 1` and `n_col*t_col - t_blk`
  for `n_block = 0`, walk the ranking, keep the prefix that fits
  `cap/(8*mo_num^2)` blocks, and the saving is the sum over that prefix. That is
  an **offline projection** - what a cap of C would have bought on this
  snapshot - and at 264 MOs the curve is steep: 155 pairs (8 %) carry 50 % of
  the value, 530 carry 80 %, 1073 carry 95 %. Stage 1 fills in arrival order
  instead, which is free and exactly right when the cap holds the whole set (975
  MiB at 264 MOs, 803 MiB at 114) and is *not* a substitute for the ranking when
  it does not: the walk touches 90 % of its pairs in the first 45 s of 4226,
  before any multiplicity exists, so arrival order carries no value information
  at all. Rank offline if the production cap binds. The report prints the
  ingredients; `pair_curve.py` (in the measurement scratch, not in the tree) is
  the 60-line version of that ranking.
* `block-only = 0` on the probe means the two key sets nest; expect that to
  hold, but the rule does not depend on it.
* If the union of the production key set fits the cap outright, there is no
  threshold to choose: cache everything and the design shrinks further.
* If both classes are small on the production run (a few % of the walk), stop
  regardless of the curve: the target is elsewhere (per-mask `mat = 0`, the
  `banned(:,:,1..2)` copies, the `get_double_excitation` re-derivation inside
  `diag_H_mat_elem_fock`).

If the curve instead shows reuse ~1 in both classes (every mask a fresh pair),
the contractions are genuinely irreducible and the cache has nothing to hoist.
That was the branch on which a device was worth costing - and the device was
costed and closed by measurement before this spec left draft: a per-call round
trip on an A100 is 26.8 us, one serialized server, against 68.5 us on one host
core but 40 host threads answering 15x more columns per second, and the only
shape the device wins is the whole `(k,l)` block in one gemm - which is the
cache's fill. The probe, its numbers and the rejected variants are in
`GPU_ROUTE.md` in the measurement scratch (`~/OpenCode/tmp/prof/`, not in the
tree). What survives of the argument, for the record:

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

Two consequences were written down before opening the fork, and the measurement
confirmed both. A host cache and a resident device copy attack the same bandwidth
problem, so they are alternatives rather than complements; only the block GEMM is
unique to the device. And Vulkan - the portability route - has no fp64: Apple's
GPUs have none at all, and most consumer devices expose `shaderFloat64` as false.
The Vulkan path is therefore fp32, which is exactly what
`mo_cholesky_double = False` already selects (`EZFIO.cfg:38-42`), and why an fp32
host build would have been needed as the oracle. Neither is built: with the
device closed by the per-call round trip and the cache taking the host walk to
~80 % saved, the fp32 build lost its only client.

Caveat for the report: the reuse structure is a property of one determinant list.
A CIPSI loop changes that list every iteration, and the cache - now filled in
arrival order rather than from a fitted set - simply keeps whatever it admitted
and stops admitting when the cap is reached. The cap is therefore the thing to
choose from a production snapshot, not a set.

## Validation (and why benzene is not a Cholesky test bench)

The cc-pVDZ benzene fixture carries 2985 AO Cholesky vectors (rank 2694 in the
MO basis) for 114 MOs, and it is *map* mode (`do_mo_cholesky = False`,
`mo_integrals_cache_shift = 7`, window 1-128, so it covers all 114 MOs). Forcing
the Cholesky path on it is not the hopeless case an early estimate suggested -
"~80x map mode, the first 256 walk positions do not complete in 8 minutes" came
from a cold window that was still paying the provider chain and the cache fill,
and it survives in no table here. Measured, one thread, fp64, on this fixture: **18 162 column calls/s at shift 6**
(7673.7 s of walk for 139 365 752 calls) and 94 234/s at shift 7 (1478.9 s for
the same calls). The shift-7 rate *decays* along the walk (93 900/s at 600 s,
87 700 at 1110 s) because `pt2_J` orders the generators most expensive first, so
a rate measured in the head overstates the whole walk; the completion numbers
are the direct ones, 2 h 09 at shift 6 and 26 min at shift 7. Both fit on an idle
laptop overnight, so the fixture validates mechanics, coverage *and* end-to-end
identity without production size; the coverage curve transfers to production
because the `(k,l)` request sequence does not depend on the window. Only the
per-block memory and the pair set do.

The shift-7 walk spent 40.9 s of its 1478.9 s in the column path - **2.8 %** -
and 1438.1 s in work no integral-side change touches. So the end-to-end saving
is capped at **80 %**, not 88 %: the 88 % is a fraction of the integral time, and
21 % of the walk is not integral work at all. That gap is exactly what the
completed A/B measures, and why it was worth two overnight runs rather than a
model.

The acceptance tests that follow from all this:

* Stage 0 mechanics: any run with requests > 0 (fixture above, killed after a
  few minutes), then inspect `work/mo_pair_cache.txt` and the summary block.
  `no request` at level 7 is expected on 114 MOs and is not a bug.
* Stage 1 arithmetic, mechanical, **on the first hit of each class only**:
  recompute through the matching `*_uncached` routine into a scratch buffer and
  compare. One extra call per class per process, and it executes at production
  scale where it matters.
  - *block hit*: bit patterns (`transfer()` to `integer(8)`) against
    `get_mo_two_e_integrals_ij_uncached` - dgemm against dgemm, so zero
    tolerance, and on mismatch `qp_bug(irp_here,...)`. This is the check that
    sees a transposed, shifted or mis-strided block.
  - *column hit*: the column norm, `||cached - dgemv|| / ||dgemv|| <= 1e-9`,
    `qp_bug` above it, and the ratio printed either way. Tolerance and not bit
    equality because a column hit returns the gemm roundings of a value the
    walk computes today with a dgemv, and the two differ at ~1e-11 relative on
    this BLAS (the measurement is in *Properties*). A layout error gives O(1)
    on the same norm, so the check keeps its power with nine orders of margin.

  This replaces the standalone 30-line program the first version of this
  section owed, and it replaces the claim that the 114-MO shift-6/7 A/B had
  already made the column comparison: it had not - both its arms call dgemv.
* Stage 1 end-to-end: run with the cap at 0, then with the cap set, on the
  same binary and fixture; compare the printed row (`PT2`, `rPT2`, `E+PT2`,
  `E+rPT2` with their errors) within 1e-9 relative, printing every deviation.
  This is an aggregate statement whose job is to catch an O(1) error, not a
  rounding one: the column class differs at ~1e-11 by construction, so the row
  cannot be required to match digit for digit. The error columns printed next to
  each number are the reference for how far apart two legitimate runs can be.
  The `PT2 Energy denominator` line
  must agree to the last digit too, but it cannot move at all: it is `<Psi|H|Psi>`
  computed from
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

Identical to all 17 printed digits, same fixture, same flags, fp64, one thread.
The terms of that contraction are O(10^2) Ha with alternating signs and the
result is -231.2, so the sum is heavily cancelling and a *systematic* or
*widespread* 1-ulp difference in the streamed columns would reach the last digit.
A single 1-ulp change in one integral would not - an aggregate is a weak
per-value test, and the per-value test is the one-shot assertion in the Stage 1
hit path (see *Bit-identity* below, and the acceptance list). What the aggregate
does rule out is the failure mode that actually threatens the design: a wrong
index or a misaligned column, i.e. an O(1) error on some pairs.

### The completed A/B: 5.08x, and a row identical to the last digit

Both arms ran to the last sample, one thread, fp64, same binary, same 46 880
determinants, `QP_PAIRBLOCK_STATS=1`. Two campaigns: v1 timed the arms on an
idle machine, v2 was launched only for the printed row - and v2's arm B shared
the machine with a 264-MO probe for ~13 min, so v1 remains the timing and v2 is
the row. v1's arm A also died in the report (see below), which is why v1 has no
row at all.

| campaign | arm | `shift` | window | walk | `PT2` wall | block | column | residual |
|---|---|---|---|---|---|---|---|---|
| v1 | A | 6 | MOs 1-64 of 114 | 7673.7 s | 7739.3 s | 1273.7 s | 4777.2 s | 1622.8 s |
| v1 | B | 7 | all | 1478.9 s | 1524.4 s | 0 | 40.9 s | 1438.1 s |
| v2 | A | 6 | MOs 1-64 of 114 | 8405.4 s | 8461.3 s | 1385.2 s | 5333.3 s | 1682.8 s |
| v2 | B | 7 | all | 1673.6 s | 1748.0 s | 0 | 52.3 s | 1621.3 s |

* **the cache's worth: 5.08x** on the `PT2` wall, 80.3 % saved, 5.19x on the
  walk. v2 measures 4.84x on a contended arm B, which bounds the same quantity
  from below; the idle-run 5.08x stands. That is the ceiling for every
  integral-side change at this fixture, and the v2 arm-B breakdown says why:
  of its 1748 s, 52.3 s are integrals and **1621 s are the walk itself** -
  mask construction, the `out_val` protocol, the comb bookkeeping. Stage 1
  should land at arm B or a little above it, and nothing that speeds up
  integrals further can go below that floor.
* **the books close to 4 %.** Price arm A's counters with the model Stage 1
  implements - 8100 distinct blocks at 0.96 ms (7.8 s of fill), 139.4M hits at
  the 0.37 us arm B measures (52 s), and arm A's own non-integral residual
  (1687 s) - and the walk comes to 1747 s against arm B's measured 1673.6 s. The
  73 s gap is arm A's per-call overhead on 1.45M block calls, which arm B never
  pays and which arm A books outside the block timer. Three measured terms, no
  fitted parameter, so the infinite-cap cache is fully described - and the
  264-MO accounting below is the same with every term re-measured at 264: it
  projected 5.9x and the gate measured 8.6x, the bookkeeping term being the one
  that moved.
* **the row is identical on all eight printed numbers**, not just the four
  energies (`ab_compare.py full_s6_T.log full_s7_T.log`, exit 0):

  ```
   PT2      -0.32987731016744182   +/-  6.4086257836338200E-004
   rPT2     -0.30537770372701645   +/-  5.9326645559783155E-004
   E+PT2    -231.53059579321010
   E+rPT2   -231.50609618676967
  ```

  Both arms issued exactly 139 365 752 column calls. Arm B served every one of
  them from a block the fill `dgemm('T','N',114,114,2694)` produced; arm A
  streamed them all through dgemv. What that excludes is the failure mode Stage
  1 can actually commit - a misindexed or misaligned column, or a *systematic*
  kernel difference between "read the block" and "stream the column" - because
  a difference of that kind is O(1) relative on many pairs and the sum being
  printed is a heavily cancelling -231.2 built from O(10^2) terms. What it does
  not exclude is a 1-ulp difference on one matrix element, which no aggregate
  can see. That case is covered by the one-shot assertion in the Stage 1 hit
  path (acceptance tests above): the first column hit recomputes the column
  through the streamed `dgemv` and compares bits - the per-value check, at the
  real shape and on the real data.
* **the per-pair counts are identical too**, and that is a separate statement.
  Both arms wrote `<ezfio>/work/mo_pair_cache.txt`, and on the `(k, l, n_col)`
  projection the two tables are byte-identical: 8144 pairs, 139 365 752 column
  calls, zero differing lines. Every pair was requested the same number of
  times. A swapped index changes those counts and the walk with them; a 1-ulp
  difference in one matrix element changes neither. So the two tests are
  complementary rather than redundant.
* **v1's missing row was a Stage 0 bug, not arithmetic.** Arm A died with
  `Fortran runtime error: Expected REAL for item 8 in formatted transfer, got
  CHARACTER` - the `column-only` line of the report, one `A` descriptor short
  (fixed). The walk had already finished and its counters were intact; the abort
  was in the report, before `print_summary`. `ab_compare.py` reporting the row
  as missing is the correct reading of that log, and v2 is the re-run.
* **a snapshot cosmetic worth knowing**: v2 arm B's header reports
  `cholesky_mo_num: 0`, because `nvec_stats` is only set by a block call and arm
  B issues none (`block_requests: 0`). The field is the first `nvec` seen, not
  the run's `cholesky_mo_num`; an all-window walk leaves it at 0.

What the A/B cannot give is a per-*value* statement, because it aggregates. The
per-value statement is now made by the code rather than by a side program: on the
first column hit in a process, recompute that column through the streamed `dgemv`
into a scratch vector and compare `transfer()`ed bit patterns. That is the
30-line Fortran program the spec used to owe - one dgemv, one dgemm, compare bits
- except it runs on the real `nvec` and the real stride, on the BLAS that will
run the walk, on data the walk actually asked for, and it cannot drift from the
production code path because it *is* the production code path. Cost: one extra
dgemv and one extra dgemm per process, 1.3 ms and 34 ms at 264 MOs.

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
  what the A/B tests in aggregate and what the first-column-hit assertion tests
  per value.

One design consequence worth stating, because it explains why the key is the
pair and not the column: a cache keyed by `(k,l,j)` would return the dgemv
result it computed, so column hits would be **bit-identical by construction** and
would even survive the capped regime's arrival-order dependence. It is not the
design, for two reasons. It cannot help the block class at all - that is 47.4 %
of the walk's thread-seconds and needs the whole matrix - so a block table would
exist anyway and this would be a *second* table, at 2.1 KB an entry against
2 MB for the block, which is the complexity this stage does not need. And the
reuse is in the pair, not the column: over the completed walk 139 365 752 column
calls land on 8144 distinct pairs, **17 113 calls per pair**, while a given
`(k,l,j)` triple is requested only by the masks that share `j = hfix`. Stage 0
does not count the triple set, so the per-triple reuse is an inference from
`both = 99.98 %` (all but 2462 of the completed walk's column calls), not a
measurement; adding a `(k,l,j)` counter is one hash and would close it. That
counter is the first thing to run if the 1e-11 column difference ever becomes
unacceptable, because it is the only route to bit-identical column hits without
paying a second gemm per pair.

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
`E+PT2`, `E+rPT2` printed to full precision against a run with cap 0 **on the
same fixture and in the same precision mode**, to within 1e-9 relative on every
number. Not to the last digit, and it cannot be: a column hit returns the dgemm
roundings of a value today's walk computes with a dgemv, and the two differ at
~1e-11 relative (measured in *Properties*; the block class is bit-exact and
asserted in-process, the column class is the one that moves). The gate is 1e-9
relative, which on a 0.87 Ha PT2 energy is 9e-10 Ha: three orders below the
1e-6 Cholesky threshold of these runs, two orders above the measured
deviation, and nine orders below the O(1) signature of a transposed or shifted
block. `ab_tol.py` prints every deviation rather than a verdict alone. Never
compare
across
`do_mo_cholesky` modes (see above: map mode differs in the 6th digit of E0, at
either precision). Which cap-0 row depends on the mode being tested: this
map-mode row is the wrong one for Stage 1, whose whole point is the Cholesky
branch - the Cholesky cap-0 row is `full_s6_T.log` (same 46 880 determinants,
`shift = 6`, one thread). Two caveats on this row: the fixture has
`h0_type = EN`, so it never enters the CFG branches - it validates mechanics and
arithmetic, not the CFG path; and one run does not prove the walk deterministic,
which is why the acceptance must regenerate its cap-0 row with the same binary
rather than trust a row from an older session. The walk itself is now known to
be deterministic given the integrals: two completed walks that differ only in
the window produced byte-identical per-pair call counts on all 8144 pairs (the
A/B above), which is a much stricter statement than repeating a printed row.

## Reference row, 264 MOs (cc-pVTZ, 4 threads)

The 264-MO fixture is the second Cholesky bench and the first one at a `mo_num`
where the block path outweighs the column path. Same contract: `pt2` on
`tz_basis_test.ezfio` (`mo_num` 264, `ao_num` 300, rank 6064, `n_det` 1428,
`do_mo_cholesky` True, `mo_cholesky_double` True, `mo_integrals_cache_shift` 6
so the window is (1,64), `h0_type` EN, `OMP_NUM_THREADS=4`,
`QP_NTHREADS_PT2=4`, `io_mo_cholesky = Read` so no AO Cholesky vector is in
memory):

```
 N_det = 1428    N_cfg = 352    < S^2 > = 9.08e-11
 E       = -230.91538522016455      (also the "PT2 Energy denominator" line)
 Variance=  3.4123389912803850      +/- 8.1684089150127015E-004
 PT norm =  0.53508620642076787     +/- 4.8551242274436647E-005
 PT2     = -0.86579718276647843     +/- 1.6078615206390628E-004
 rPT2    = -0.67308215287443873     +/- 1.2499727596452631E-004
 E+PT2   = -231.78118240293102      +/- 1.6078615206390628E-004
 E+rPT2  = -231.58846737303898      +/- 1.2499727596452631E-004
 WALL 4279.1 s, CPU 10871.3 s, peak RSS 4.2 GB, walk clock 4226.0 s
```

This is the cap-0 row for a Stage 1 A/B at 264 MOs. Read it with the tree's own
caveat: TZ total energies are good to ~4e-4 Ha (the basis-normalization residual
of that tree, see the fixture notes), which is the scale of the PT2 error bar, so
the row is a *mechanical* reference - it is the printed digits that Stage 1 must
reproduce, not a number to compare with any other calculation.

### The warm build, and what the counters can see

`pt2` arms the cache at `fci/pt2.irp.f:61` (`PROVIDE pair_cache`), calls
`PROVIDE psi_energy` at `:65`, and arms the counters at `:66`
(`pair_block_stats_init`). `psi_energy` (`davidson/u0_hs2_u0.irp.f:3`)
builds `h_matrix_all_dets` through `u0_hs2_u0`, whose `!$OMP DO` at `:962` calls
`get_mo_two_e_integrals_ij` for every `(k,l)` whose two orbitals are neither
always occupied nor always virtual - a superset of the pairs the walk asks for,
each computed once, in parallel, 27 s into the run.

Those dgemms are real work and they are in neither the 235888 nor the 8007.76
thread-s: they precede the counters. A run therefore performs more block dgemms
than Stage 0 counted. The recomputation fraction is unchanged, because a warm
block is computed once and the walk's is computed `235888/1833 = 129` times.

The warm build is also why the capped arm's stats block looks degenerate. It
asks for the whole union (1028 MB of blocks against the gate's 2048 MiB cap), so
by the time the walk starts every block and every column it will ever ask for is
already in the table. The walk then misses zero times in 235888 block requests
and once in 4.68M column requests - the once being the column assertion's own
dgemv, the only call in the run that reaches an `*_uncached` body after
`pair_block_stats_init`. The block assertion fires earlier, during the warm build
itself, and its reference dgemm is uncounted like the rest. The report reads

```
 PT2 pair-block stats: mo_num=264 cholesky_mo_num=0
   block  calls=0 distinct=0 ...      column calls=1 distinct=1 ...
```

and that is a correct measurement of a pure-hit walk, not a broken counter: the
counters count *misses*, so a run with none reports none. It also settles the
fill question on this fixture. The warm build computes those 1844 blocks whether
or not the cache exists (arm A pays them too), so the cache's marginal cost here
is the copy, not the fill; at production the fill replaces the first of
`235888/1833 = 129` recomputations. What the cached arm spends that arm A does
not is the block copy - 235888 x 557568 B = 131 GB of memcpy, plus 4.68M column
reads of 2112 B - which the mid-walk sample puts at ~3 % of the capped walk's
1369 s of CPU: the hit path is not what a capped run spends its time on.

The cap must sit *below* the union for the counters to say anything, because the
rate they expose is a miss rate. At 512 MiB this fixture gets `n_pair_cache`
962 slots for 1833 blocks: the pool freezes during the warm build, the walk
misses and recomputes, and `block calls` / `block_thread_s` become the miss count
and the fill time. That arm, and the two it is compared against, are the record
below.

### Stage 1 gate at 264 MOs

Four arms, two binaries - the three of 21 November 11:55 and the deduced cap on
the build that carries it - one fixture, `shift = 6`, four threads,
`QP_PAIRBLOCK_STATS=1`, `tz_basis_test.ezfio`:

```
 cap        slots   block req. (miss)   column req. (miss)      WALL        CPU    vs cap 0
 0             -        235 888             4 684 818        3698.2 s   11 927.3 s     1
 512         962         10 078             1 808 581        1174.7 s    3 330.8 s   3.15x
 2048       3851              0                     1         428.4 s    1 368.9 s   8.63x
 auto(3584) 6740              0                     1         372.6 s    1 300.0 s   9.93x
```

With the cache off every request is a call; with it on, a request that misses is
the only one the counters see. The row, against the 1e-9 gate (`ab_tol.py`
prints every deviation; none is suppressed):

```
                          cap 0                   cap M     rel. dev
denom        -230.9153852201646      -230.9153852201646     0.00e+00
PT2        -0.86579718276647843    -0.86579718276647843     0.00e+00
rPT2       -0.67308215287443873    -0.67308215287443884     1.65e-16
E+PT2       -231.78118240293102     -231.78118240293102     0.00e+00
E+rPT2     -231.58846737303898     -231.58846737303898     0.00e+00
```

Seven of those nine numbers are bit-identical - the denominator, `PT2`, `E+PT2`,
`E+rPT2` and their four `+/-` - and both of `rPT2`'s move one ulp: the value by
1.65e-16 and its `+/-` by 2.17e-16 (...2631E-004 against ...2633E-004). That is
1.6e-16 against a 1e-9 gate, seven orders of margin, and it is the column class
moving as predicted, not a layout error. Both in-process assertions fired at the
production shape and stayed silent:

```
Pair cache: first block hit, (7,8), bit-identical to the uncached block
Pair cache: first column hit, (20,21,23), relative deviation 6.8296E-16
```

**The deduced cap is the same run, not a nearby one.** With `QP_PAIR_CACHE_MB`
unset the provider gives 3584 MiB on this 7 GB box - 6740 slots against a
1844-pair demand - so the pool never freezes, and all nine printed numbers are
bit-identical to the 2048 MiB row: a `diff` of the five row lines is empty.
That includes `rPT2`'s ulp, which *is* the column-class difference and is
therefore absent between two arms that both serve every column from a block.
The deduction was checked to follow the budget rather than return a constant:
`QP_MAXMEM=3` gives 1536 MiB and 2888 slots, still above the demand, and both
arms fill and hit with the block assertion silent. A malformed value is
announced - `QP_PAIR_CACHE_MB=2G` prints `QP_PAIR_CACHE_MB="2G" is not an
integer; cache disabled` - because a silent fall back to 0 turns a typo into the
whole 8.6x and leaves nothing in the log to find it by.

**The freeze branch executed.** It is not announced, and it does not need to be:
a pair's first miss claims a slot and publishes its key, so a pair that misses a
*second* time was refused at `pair_cache_claim`'s `:503`. The 512 MiB arm has 801
pairs with two or more block misses (max 103) and 106 with exactly one, against
962 slots, so the branch ran and the pool stopped admitting partway through the
warm build.

**A cap worth 52 % of the union buys 3.15x of the 8.63x.** 962 slots against
1833 pairs keeps 95.7 % of the block requests: 926 pairs hit every time and
absorb 225 810 of the 235 888, because the warm build admits in its `(k,l)` order
and the walk's block demand is concentrated there. The residual is column work
the pool cannot serve - a column is a slice of a block, so a pair whose block was
never admitted pays for every one of its column calls. The 907 uncached pairs
carry 1.81M column calls at 1009 us, which is 1825 s of that arm's 3331 s of CPU.
This is the measured case for both remedies already specified: admission order,
and the replacement policy of *When the cap binds*.

**Arrival-order admission did not reach the row.** The two capped arms admit
different sets (962 slots against 3851) and serve 2.87M against 4.68M columns
through different arithmetic, and all nine printed numbers - the denominator,
the four values and their four errors - are bit-identical (`ab_tol.py` reports
0.00e+00 on every line of that pairing). Cap 0 differs from both by the one ulp
above. The caveat in *Properties* is therefore unexercised on this fixture, not
disproved, and the 1e-9 gate is what carries the claim.

**The accounting, corrected.** The class timers are wall-clock windows, so on a
loaded box they contain the time a thread spent off-CPU and over-count its
running time. One number settles all three arms: the walk's own bookkeeping. The
2048 MiB arm removes both integral classes (zero block misses, one column call)
and its walk still costs 1369 s of CPU, of which the hit path is ~3 % - a
mid-walk `sample(1)` puts 109 + 107 of 7770 samples in the two wrappers and the
rest in `fill_buffer_double` and `splash_pq` - so the copies are ~50 s and the
bookkeeping ~1300 s. With that, and one free parameter, arm A closes: 11 % of its
11 866 thread-s of timers spent preempted leaves 10 600 s of BLAS, which with
1300 s of bookkeeping and ~30 s of startup is 11 927 s of CPU, 0.3 % from
measured. The shares of *Stage 0* rank the work; they do not partition CPU time.
And what the cache left behind is now the walk: at 264 MOs it is
bookkeeping-bound, so further speedup on this fixture comes from
`fill_buffer_double`, not from more caching.

**Not gated: wall time.** `WALL` is the clock, and this box was under load 37 for
the cap-0 arm and nearly idle for the capped ones, so every ratio below is an
upper bound: the reference is the inflated arm. The scatter within the capped
side is large enough to matter - three runs of the zero-miss walk took 428.4,
440.1 and 372.6 s of wall and 1368.9, 1437.3 and 1300.0 s of CPU, for work that
computes no integral element differently (18 % on wall, 10 % on CPU). CPU time is
the steadier number and the one to quote: **9.18x** against the cap-0 arm with
the deduced cap, 8.71x with the 2048 MiB cap, so ~9x with a few per cent of
slack. The wall column's 9.93x is the same cache measured on a quieter box.
