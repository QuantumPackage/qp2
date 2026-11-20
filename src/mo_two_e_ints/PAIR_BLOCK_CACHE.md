# Static hot-block cache for the Cholesky MO-integral block (PT2 selection)

Spec for implementation. Stage 0 is a measurement with no behaviour change and
decides whether Stage 1 is written at all.

## The cost being attacked

In production (`mo_two_e_ints/do_mo_cholesky = True`, `mo_num` > 2^`mo_integrals_cache_shift`):

* `mo_integrals_cache` covers only a window of 2^shift orbitals centred on the
  Fermi level (`map_integrals.irp.f:55-62`; level 7 = 128 orbitals, 2 GB), so
  `get_mo_two_e_integrals_ij` takes the GEMM branch (`map_integrals.irp.f:662-695`).
  The sparse MO map is **not built** in this mode (`mo_bi_integrals.irp.f:65-82`),
  so every scalar integral outside the window is a ddot over the Cholesky
  vectors (`map_integrals.irp.f:250-292`) and every vector call is a dgemv
  (`:370-389`, `:539-571`).
* The GEMM branch computes a full `mo_num x mo_num` block:
  `2 * mo_num^2 * cholesky_mo_num` flops. At 300 MOs with the production count
  of 1000-2000 vectors that is 0.18-0.36 GFlop: 4-9 ms on one core at
  40-60 GFlop/s, and compute-bound (4-8 MB of operands, ~45 flop/byte). The
  `mo_num^2` scatter the block feeds costs ~0.8 ms at the same size
  (`tau ~ 8.5 ns` per element), so the GEMM is ~85 % of an `nt == 2` call.
  Stage 0 measures that share rather than trusting this estimate.
* Only `get_d0` calls it, once per invocation: two static sites, one per branch
  of the `sp == 3` test (`cipsi/selection.irp.f:1551` for AB, `:1576` for
  AA/BB). `get_d0` is entered only for `nt == 2` selectors
  (`selection.irp.f:1035-1036`). By contrast `get_d1` (nt == 3) issues eight
  `get_mo_two_e_integrals` dgemvs (`:1305-1440`) and `get_d2` (nt == 4) calls
  neither routine. So the block GEMM is the exclusive property of the
  **nt == 2 subset**: `pair_block_total` *is* the nt == 2 selector count, and
  no second counter is needed to know that share.
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

File: `map_integrals.irp.f`, inside `get_mo_two_e_integrals_ij`, in the
`do_mo_cholesky` GEMM branch only (`:664-695`), i.e. after the
`sze /= mo_num` check.

State, as IRPF90 providers in the same file (zero-initialised in their block):

```
&BEGIN_PROVIDER [ logical, pair_block_stats_on ]      ! env QP_PAIRBLOCK_STATS
&BEGIN_PROVIDER [ integer, pair_block_count, (mo_num*mo_num) ]
&BEGIN_PROVIDER [ integer(8), pair_block_total ]
&BEGIN_PROVIDER [ double precision, pair_block_time ]  ! CPU seconds in the block
```

`pair_block_count` is a flat array indexed by `key = k + mo_num*(l-1)`, the same
key Stage 1 hashes, so no second encoding exists. Inside the branch:

```
    call wall_time(t0)
    !$OMP ATOMIC UPDATE
    pair_block_count(k + mo_num*(l-1)) = pair_block_count(k + mo_num*(l-1)) + 1
    !$OMP ATOMIC UPDATE
    pair_block_total = pair_block_total + 1
    <the existing code, unchanged: GPU probe or dgemm/sgemm>
    call wall_time(t1)
    !$OMP ATOMIC UPDATE
    pair_block_time = pair_block_time + (t1-t0)
```

`wall_time` is already called this way from this directory
(`mo_bi_integrals.irp.f:62`, no `use` needed). The timing covers the whole
`if (mo_cholesky_double)` body so the GPU and CPU paths are comparable, and
`pair_block_time` is **CPU** time summed over threads: it yields the achieved
GFlop/s directly, and dividing it by `walk_time * nthreads` gives the share of
the walk. Two clock reads per millisecond of work is free.

`pair_block_stats_on` is set once from `get_environment_variable`; when false the
branch skips the atomics (`call wall_time` included). The gate is there so a
debugging run is untouched - the branch only executes in the production regime,
where a call costs milliseconds.

Report: `print_pair_block_stats(walk_time)` in `map_integrals.irp.f`, called from
`fci/pt2.irp.f` guarded by `pair_block_stats_on`, with `walk_time` measured by
`wall_time` around the `pt2_stochastic` call (`:97-100`).

It prints, and nothing more:

```
PT2 pair-block stats: calls = <total>, distinct = <#count>0>, reuse = <R>
  block GEMM: <pair_block_time> CPU-s = <share>% of walk CPU-s (<nthreads> threads)
              mo_num = <n>, cholesky_mo_num = <m>, <F> GFlop, <Q> GFlop/s
  top-M      coverage      memory
    100        ..%          <100*8*n^2 MB>
    ...
```

where coverage(M) = sum of the M largest counts / total. Sort a scratch copy of
`pair_block_count` with `isort_c(A, iorder, isize)` from
`utils/qsort_module.f90`: it sorts `A` **in place** and returns the permutation
in `iorder`, which must be initialised to `1..isize` first (`src/utils/qsort.c`).
Walk `m = isize` down to `isize-M+1` for the largest counts and recover the pair
as `key = iorder(m)-1; k = mod(key,mo_num)+1; l = key/mo_num+1`.

It also writes `<ezfio_filename>/work/mo_pair_cache.txt`: one line per pair with
count > 0, `k l count`, ordered by decreasing count, truncated to the largest M
that fits `mo_cholesky_pair_cache_mb` (see Stage 1). Text, not binary: it is
read once and must be inspectable. Only the master rank reaches the report, so
no MPI collective and no file race.

Acceptance: `pt2_data_I*` and the printed energies are bit-identical with and
without `QP_PAIRBLOCK_STATS=1` on the same fixture. A report of
`calls = 0` is not a bug: it means the cache window already covers all MOs
(level 7 with `mo_num <= 128`), so the GEMM branch never runs and there is
nothing to cache at that size.

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
`mo_cholesky_pair_cache_mb > 0` and `work/mo_pair_cache.txt` exists. Read the
first `n = min(lines, mo_cholesky_pair_cache_mb*131072/(mo_num*mo_num))` lines
(1 MB = 131072 doubles), fill `pair_cache(1,1,i)` with `!$OMP PARALLEL DO` over
i, insert `key = k + mo_num*(l-1)` into the open-addressing table (linear
probing, `slot = mod(key, pair_cache_nhash)+1`). The build is deterministic and
happens once.

Note the ordering constraint: Stage 0 truncates the file to the cap that was
set *during that run*, so the cap must already hold the intended value when the
measuring run is launched - otherwise the file is empty (cap 0) or shorter than
wanted. `n_pair_cache` is then re-checked against the cap at build time, so a
stale or oversized file cannot overflow the array.

The cache should be PROVIDEd in the pre-loop block of `pt2_omp_task`
(`cipsi_utils/pt2_omp.irp.f:272`, next to `PROVIDE psi_det_hii ...`). The
generated `provide_*` routines are lock-protected (`call irp_lock(...)`), so a
lazy build would still be correct, but it would stall every worker at the first
mask and hide the build cost: pre-provide it so it is paid once, before the
walk, and shows up in the log. IRPF90 broadcasts providers across MPI ranks, so
an MPI-launched run pays one broadcast of the cache; the PT2 walk itself is
OpenMP-parallel in one rank.

Split the routine:

* `get_mo_two_e_integrals_ij_uncached(k,l,sze,out_array,map)` = today's body,
  unchanged, all branches (GPU, dgemm, sgemm, map).
* `get_mo_two_e_integrals_ij(k,l,sze,out_array,map)` = `PROVIDE pair_cache`,
  then: probe; hit -> `out_array(1:sze,1:sze) = pair_cache(:,:,slot)`; miss ->
  call the uncached routine. `sze == mo_num` is already asserted (`:658`).

Properties:

* Bit-identical by construction in both precision modes: the block is always
  produced by the one existing routine, with the same operands and shape. No
  second implementation of the integral computation.
* The array is read-only after the build: no locks, no eviction.
* On a hit the cost is one `8*mo_num^2` byte copy: 720 KB at 300 MOs, ~50 us,
  about 1 % of the 4-9 ms dgemm it replaces. Every hit pays, and it never
  matters.
* A miss is exactly today's cost plus a hash probe.

Not in scope: LRU/eviction, a traced key sequence, a runtime-adaptive hot set.
If the static set turns out to cover too little, the fork is not more caching
machinery - it is the batched device GEMM discussed under the decision rule.

## Decision rule

Stage 0 reports three numbers for the production wavefunction: `calls`,
`distinct` (reuse `R = calls/distinct`) and the CPU-seconds in the block, hence
the block's share `s` of the walk. Everything else is arithmetic:

* cached block cost = `distinct * t_gemm + calls * t_copy`, uncached =
  `calls * t_gemm`. Since `t_copy ~ 1 %` of `t_gemm`, the block term shrinks by
  ~`R`, and any reuse above 1.05 pays. **Whether to cache is not the question.**
* the question is how much of `calls` fits the memory cap, i.e. the coverage
  curve at a given `M`. A cache holding the top-M pairs leaves the block term at
  `1 - C(M)` of its uncached value, so the walk saves `s * C(M)`.
* if `s` itself is small on the production run (say < 5 %), stop regardless of
  the curve: the target is elsewhere (per-mask `mat = 0`, the
  `banned(:,:,1..2)` copies, the `get_double_excitation` re-derivation inside
  `diag_H_mat_elem_fock`).
* otherwise choose the smallest `M` whose coverage buys more than its memory
  costs, and report `M`, its GB, and the projected saving.

If the curve instead shows `R ~ 1` (every mask a fresh pair), the block is
genuinely irreducible *and* the case for batching it on a device becomes the
strongest one available: `calls` independent GEMMs of 0.18-0.36 GFlop each, one
720 KB download apiece. That is the GPU fork, and it should be opened with `s`
in hand, not before.

Caveat for the report: the hot set is fitted to one wavefunction. In a CIPSI
loop it must be refreshed by re-running Stage 0, and the curve's decay with `M`
is the only evidence about how fast such a set goes stale.

## Validation (and why benzene is not a Cholesky test bench)

The benzene fixture carries 5918 Cholesky vectors: one block dgemm is ~4 ms and
one scalar ddot ~6 us. Forcing a partial window on it
(`mo_integrals_cache_shift = 5`, window 32) turns the ~10^7 block calls and
~10^9 scalar calls of the walk into tens of CPU-hours, so it can only validate
mechanics if the run is stopped early and only the counters and the file are
read.

* Stage 0 mechanics: any short run with `calls > 0` (fixture above, stopped
  after a few minutes), then inspect `work/mo_pair_cache.txt` and the coverage
  table. `calls = 0` at level 7 is expected on 114 MOs.
* Stage 1 bit-identity, mechanical: in the hit path, **on the first hit only**,
  recompute the block through `*_uncached` into a scratch buffer and compare bit
  patterns (`transfer()` to `integer(8)`); on mismatch `qp_bug(irp_here,...)`.
  One extra dgemm per run, and it executes at production scale where it matters.
* Stage 1 end-to-end: production run with the cap at 0, then with the cap set;
  compare the printed PT2 energies and `pt2_data_I*` byte for byte. The walk is
  schedule-independent by contract (`pt2_omp.irp.f` header), so this can be done
  at the production thread count.

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
