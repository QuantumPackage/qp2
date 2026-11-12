# Golden corpus for the QP2 tool-layer port (P0, T0.1–T0.3)

Reference EZFIO fixtures, `qp_edit -c` transcripts, and the md5 corpus
against which the Rust port is gated. The normative byte-level contract
is [`../../CONTRACTS.md`](../../CONTRACTS.md); the normative check-mode
behaviour is [`contracts/check-mode.md`](contracts/check-mode.md).

## Layout

```
inputs/       5 geometry files (xyz / Gaussian Z-matrix)
fixtures/     6 EZFIO directories produced by the reference OCaml build
              (~/qp2 @ be4b71c1, pre-fix). Volatile keys stripped,
              .gz recompressed mtime=0 level 9. hcn has converged SCF MOs.
transcripts/  per fixture: stdout, stderr, exit of `qp_edit -c` on the
              fixture (reference build = EVIDENCE), mutated.txt (the
              mutation log), added/ (decompressed payload of every file
              the pre-fix build created; the .gz suffix is stripped).
              pre-fix-quantisation/probe.log: the coordinate-quantisation
              drift measured on the pre-fix build (contract C5).
md5/          basis-md5.jsonl: (basis, element, pseudo, n_ao, ao_md5,
              ao_extra_md5, basis_file_sha256) for 27 bases; last element
              of each basis file.
contracts/    hand-written, NEVER regenerated: check-mode.md (C1–C7) and
              check-mutated.json (per fixture, per contract post-fix /
              pre-fix: the exact set of files qp_edit -c must create).
manifest.json generator provenance: fixtures, scenarios, normalisation.
```

Fixture set: `h2` (cc-pvdz), `h3` (cc-pvdz, multiplicity 2), `c2h2`
(cc-pvdz_ecp_bfd + `-p bfd`), `b2` (6-31g, Z-matrix), `h2g` (cc-pvdz,
ghost X at `-g 5.0`), `hcn` (aug-cc-pvdz, `qp_run scf`).

## Commands

All scratch goes to `$HOME/OpenCode/tmp` (never `/tmp`).

```
H=rust/tools/golden-harness/golden_harness.py
C=rust/tests/golden

# Regenerate the corpus with a tool set (~20 s, ~30 MB scratch):
python3 $H generate --qp-root ~/qp2 --corpus $C

# Replay the manifest with two tool sets and diff the snapshots:
python3 $H run --qp-root ~/qp2 --corpus $C --out $HOME/OpenCode/tmp/golden-a
python3 $H run --qp-root <rust-build> --corpus $C --out $HOME/OpenCode/tmp/golden-b
python3 $H diff $HOME/OpenCode/tmp/golden-a $HOME/OpenCode/tmp/golden-b

# Gate a tool set against the corpus:
python3 $H verify --qp-root ~/qp2 --corpus $C                     # post-fix target
python3 $H verify --qp-root ~/qp2 --corpus $C --contract pre-fix  # reference build
python3 $H verify --qp-root ~/qp2 --corpus $C --only check        # scenarios subset
```

`--qp-root` may be any build tree (OCaml: binaries in `$QP_ROOT/ocaml`;
Rust: `$QP_ROOT/ocaml` or `$QP_ROOT/bin`); `--qp-bin` overrides the
search. `verify create` compares the recreated tree byte for byte
(decompressed); `verify check` compares the transcript and enforces the
contract's created-file list plus `added/` payloads. A snapshot stores
stdout, stderr, exit and the decompressed tree per
`(fixture, scenario)`; `diff` exits 1 on any difference.

## Normalisation policy (what is dropped, and why)

- `ezfio/creation` (date), `ezfio/user` (login), `ezfio/library` (absolute
  build path): written by every tool run, no contract. `.version` is kept.
- gzip containers: decompressed payload is the oracle (plan D6); the
  corpus and snapshots store it with mtime=0 level 9, snapshots as text.
- Machine lines from tool output (see `NOISE` in the harness): `Date:`,
  `Wall time:`, Git stamps, memory budget, progress.
- `work/` scratch directory excluded from snapshots.
- `OMP_NUM_THREADS=1` for every run (threaded BLAS reductions are not
  bit-reproducible).

## md5 collisions (expected, by construction)

The last element of `6-31g_star` and `6-31g_star_star` (KRYPTON, 70 AOs)
is byte-identical in the shipped basis data, as are the RADON tails of
`def2-svp` and `def2-tzvpp` (118 AOs). Two rows per pair share
`(ao_md5, ao_extra_md5)`: `dc7fb1…` / `6b108b…`. These are data
collisions, not a defect.

## provenance of `transcripts/<f>/added/` vs `contracts/`

`transcripts/` (including `added/`) is *evidence* recorded from the
pre-fix reference build `be4b71c1`; `contracts/` is the *normative*
target transcribed from commit `32c9899d` (post-fix). For `hcn` the
difference is 9 vs 7 created files (the two `*_qp_edit` scratch files).
`generate` must never touch `contracts/`.
