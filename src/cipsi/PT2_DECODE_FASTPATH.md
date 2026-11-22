# Decode fast path for `diag_H_mat_elem_fock` (PT2 selection)

Implemented and gated 2026-11-21/22. Companion to
`src/mo_two_e_ints/PAIR_BLOCK_CACHE.md`.

## What it does

`fill_buffer_double` already holds the two movers `(h1,p1)` and `(h2,p2)` of the
candidate it is about to score, but handed the *determinants* to
`diag_H_mat_elem_fock`, which recovered those movers with
`get_excitation_degree` + `get_double_excitation` + `decode_exc`. On the 264-MO
fixture those calls were 32 % of the PT2 walk.

`diag_H_mat_elem_fock_exc(h1,p1,h2,p2,s1,s2,fock_diag_tmp)` computes the same
Fock diagonal from the movers. The three arithmetic branches are **copied
verbatim** from `diag_H_mat_elem_fock` (alpha/alpha, beta/beta, alpha/beta),
including the `jj_anti` / `jj` distinction. `diag_H_mat_elem_fock` is not
touched and remains the fallback.

## Contract

The caller takes the fast path only when

```fortran
ok_h .and. ok .and. p1 /= h1 .and. p1 /= h2 .and. p2 /= h1 .and. p2 /= h2
```

that is, when the hole application and the particle application both succeeded
and the resulting determinant is a genuine double excitation. Everything else
falls back to the frozen routine.

Why each term is needed:

* `p1 == h1` (or `p2 == h2`) is **reachable by design**: `bannedOrb(h1,s1)` and
  `bannedOrb(h2,s2)` are reset at `selection.irp.f:529/532`, so a hole can be
  reused as a particle. Such a candidate is a degree-1 or degree-0 excitation
  and its `Hii` must come from the frozen function's degree branch; computing
  the degree-2 formula for it moves the weight by ~1e-2, not 1e-11.
* `p1 == h2` / `p2 == h1`: the pair cancels and the candidate is degree 0,
  same story.
* For `sp == 3` the caller fixes `s1 = 1, s2 = 2`; for `sp /= 3`, `s1 == s2`.
  The same-spin branch is symmetric under 1 <-> 2, so the pair order cannot
  change the value. The routine `STOP`s on any other order rather than
  guessing.
* The `$IS_DOUBLE` guard keeps the block out of the `single` instantiation of
  the template (`fill_buffer_single` keeps calling the frozen routine);
  verified by zero `bl` to the new symbol in `_fill_buffer_single_`.

## Traps

* **`fock_diag_tmp` is `(2, mo_num+1)`.** `E0` lives in `(1, mo_num+1)`, not
  `(1, mo_num)`. Other dummies in this file declare `(mo_num)`; do not
  "harmonize" either side — a `(mo_num)` dummy reads a different `E0` (~1e-2).
* **Do not merge the twin back into the frozen function.** The frozen routine
  is the bit-exactness proof: its instruction stream must stay byte-identical,
  so the fast path cannot become an `if` inside it.
* **Do not reorder the arithmetic.** Under `-Ofast` a sum's evaluation order
  follows the operands' SSA ranks; the branches are verbatim for that reason,
  and a re-association shows up as a 1e-11-class row move.
* The last branch is a `STOP`, not a silent fall-through: an unexpected spin
  order must not produce a number.

## Gate (264 MO, `tz_basis_test.ezfio`, 4 threads, deduced pair-cache cap)

| quantity | before (`f393926a`) | after | delta (abs) |
|---|---:|---:|---:|
| denom | -230.9153852201646 | -230.9153852201646 | 0 |
| PT2 | -0.86579718276647843 | -0.86579718275149187 | 1.499e-11 |
| PT2 err | 1.6078615206390628e-04 | 1.6078615642958863e-04 | 4.366e-12 |
| rPT2 | -0.67308215287443884 | -0.67308215286269746 | 1.174e-11 |
| rPT2 err | 1.2499727596452633e-04 | 1.2499727935844855e-04 | 3.394e-12 |
| E+PT2 | -231.78118240293102 | -231.78118240291604 | 1.498e-11 |
| E+PT2 err | 1.6078615206390628e-04 | 1.6078615642958863e-04 | 4.366e-12 |
| E+rPT2 | -231.58846737303898 | -231.58846737302724 | 1.174e-11 |
| E+rPT2 err | 1.2499727596452633e-04 | 1.2499727935844855e-04 | 3.394e-12 |

`n_det = 1428` (494 alpha / 494 beta) in both arms, and each arm is
bit-identical to itself on rerun. The criterion is 1e-9 **absolute** Hartree on
all nine quantities: they span 0.86 to 231 Ha and the err bars are 1.6e-4, so a
relative line would demand 1e-13 of an err bar while accepting 1e-11 of the
energy. The move is the reassociation class expected from any edit to this
frame, below the 6.4e-11 measured for adding one argument to it.

## Speed (paired A B B A, 2026-11-22)

| run | arm | wall s | CPU s | util |
|---|---|---:|---:|---:|
| A1 | before | 237.2 | 862.5 | 3.636 |
| B1 | after | 227.8 | 754.6 | 3.313 |
| B2 | after | 287.4 | 911.8 | 3.173 |
| A2 | before | 336.5 | 1123.8 | 3.339 |

CPU ratio 1.19x on the means, adjacent pairs 1.143 and 1.232, so the claim is
the band **1.14-1.23**; wall 1.11x is derived. Utilization is 3.49 vs 3.24
(7.5 % apart). Within-arm drift was +30 % (A) and +21 % (B), larger than the
effect: the A B B A symmetry cancels its linear part, but the band, not the
point, is what the four runs support.

Structural results, noise-free:

* `_diag_h_mat_elem_fock_` inclusive 70,723 ms -> 34 ms (45.69 % -> 0.03 % of
  the process profile); `get_double_excitation` and `get_excitation_degree`
  collapse to rounding.
* the new symbol holds 7,612 ms (5.92 %).
* the caller's own self time rises 38.0 -> 60.4 s in the same 60 s windows —
  the per-candidate predicate, not spills (frame 0x4a0 -> 0x4b0, 16 bytes).
  That rise is what keeps the leg at 1.19x instead of the 1.66x naive bound.

## Re-verifying

1. **Static.** Per-symbol normalized disassembly against the parent binary. Hard
   items: `_diag_h_mat_elem_fock_` byte-identical, caller FP operation sequence
   unchanged with FPR names stripped, zero `bl` to the new symbol from the
   single instantiation. Whole-function re-registration of the callers is
   expected, not a defect.
2. **Row.** Nine-quantity table against the parent binary, 1e-9 absolute, count
   unchanged, rerun bit-identical.
3. **Speed.** A B B A with the binary sha256-verified before each run; CPU time
   primary, utilization as the contention control, no wall-clock point claim.
