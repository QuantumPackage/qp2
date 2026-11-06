# `qp_edit -c` — normative contract

Status: **normative for the Rust port**. Transcribed from the source of
commit `32c9899d` ("qp_edit: do not rewrite the EZFIO in check mode, and
do not quantize coordinates", also `15e725aa` on `origin/qp-edit-fix`),
not inferred from the reference binary. The reference build at
`be4b71c1` predates that commit; its observed behaviour is recorded
separately in `transcripts/` and is *evidence*, not the target.

## What `qp_edit -c <ezfio>` does after the fix

1. open the EZFIO, clear `qpstop`/`qpkill` under `work/`;
2. repair the AO / AO-extra basis ordering if it is not already ordered
   (`Warning: Basis set is not properly ordered. Redordering.` on stderr,
   then `Input.Ao_basis.write`, and `Input.Mo_basis.write` when MOs exist);
3. render every section to RST with the section getter (`get`), which
   materialises missing defaults along the way;
4. read the RST back and parse every section (`set ~apply:false`): a
   malformed section is still detected, but nothing is written;
5. exit 0.

## Rules

**C1 — exit status and streams.** On a valid directory the exit status is
0, stdout is empty, stderr is empty.

**C2 — no rewrite.** No file that existed before the run is modified or
deleted: every such path has the same decompressed payload afterwards.
This includes the `ao_md5`, `ao_extra_md5` and `mo_basis/ao_md5` stamps,
which the pre-fix build rewrote on every read with identical content and
which D4 of the plan forbids the Rust port from writing at all. The
exceptions are the in-place repairs the reference readers themselves
perform: step 2 (an unordered basis is reordered and written back, a
repair, not a round trip), and a stale-value repair inside a getter —
the known one is `read_state_average_weight`, which rewrites
`state_average_weight.gz` when its length no longer equals
`determinants/n_states`. Every corpus fixture is at its fixed point, so
no repair fires in the corpus and C2 holds as stated; the repair is
pinned by the Rust check-mode suite
(`a_stale_declared_shape_is_read_and_repaired`).

**C3 — missing defaults are materialised.** Keys the getter needs but
that are absent from the directory are written, because
`bin/qp_reset:114` relies on `qp_edit --check` leaving a usable
directory after it deletes the determinants. The set is fixture
dependent: it is exactly the keys of the `Determinants_by_hand` section
(`read ~full:false`) that the fixture lacks. The other five fixtures
have no MOs, so that `read` returns `None` after `read_n_det` and
`read_n_det_qp_edit`; those two keys are already present (every fixture
was produced by `qp_create_ezfio`, which ends with its own `qp_edit -c`),
so they get nothing. For `fixtures/hcn`, which has converged SCF MOs, it
is seven files:

| path | content |
| --- | --- |
| `determinants/bit_kind` | `8` (64-bit integers) |
| `determinants/expected_s2` | `0.0` (`2S(S+1)` of the closed shell) |
| `determinants/mo_label` | `Canonical` (copied from `mo_basis/mo_label`) |
| `determinants/n_int` | `1` (from `mo_basis/mo_num` = 14) |
| `determinants/psi_coef.gz` | rank 2, `(1,1)`, one coefficient `1.0` |
| `determinants/psi_det.gz` | rank 3, `(1,2,1)`, the HF determinant (`127`, `127`) |
| `determinants/state_average_weight.gz` | rank 1, `(1)`, `1.0` |

`transcripts/<fixture>/added/` holds the reference payload of every file
the pre-fix build created, decompressed and with the `.gz` suffix
stripped; `verify` compares these bytes when it knows them.

The value written for each of these keys is not the `EZFIO.cfg`
`default:` line the editor's config copies: the generated readers take
it from `$QP_ROOT/data/ezfio_defaults/` at run time
(`Qpackage.get_ezfio_default`), and the editor ports that scan, with the
config value as a fallback when the data tree is not installed. The
corpus is at its fixed point, so the back-fill is a **no-op** on all six
fixtures and on `hcn` beyond the determinants record above. It is
observable on a fresh `qp_create_ezfio` tree: before this fix the
reference created 8 files the port did not (the `qmcchem` section, which
the port did not visit, and the three `ao_extra_basis` logicals plus its
MD5, which the handler never reached), while the other materialised
keys already agreed because their config defaults equal the installed
ones. With the fix the two trees agree file-for-file (275 vs 275) and
payload-for-payload (273 common files identical, 2 machine stamps
skipped).

**C4 — the `*_qp_edit` scratch pair is not written.**
`determinants/psi_det_qp_edit.gz` and `determinants/psi_coef_qp_edit.gz`
belong to the section writer (`write ~force:false`), which `-c` skips.
The pre-fix build created them (identical to `psi_det`/`psi_coef`); the
post-fix build must not. This is the 9 → 7 difference between the two
lists in `contracts/check-mutated.json`.

A pair that already exists is a separate case and must **survive**: no
deletion, no rewrite. The pre-fix driver called `set` on every section even
in check mode (the post-fix commit adds `apply:false`), so its determinants
writer re-ran `write ~force:false`; whenever `n_det <= n_det_qp_edit` that
re-creates the pair from the parsed `psi_det`/`psi_coef` — a stale pair is
refreshed and even a consistent one is rewritten. The post-fix build applies
nothing, so the pair keeps its bytes. The pair is read (never written) when
`n_det_qp_edit != n_det`: that is the truncated-wave-function state the
editor leaves behind, and the getter must still render its section from it.
The residual is that a pair which disagrees with the working arrays stays
stale after `-c`; it is accepted in `DEVIATIONS.md` B9, and the Rust suite
pins a consistent pair, a stale pair and the read-only branch.

**C5 — no coordinate quantisation.** Coordinates are read and written as
bohr, `Point3d.to_string` prints `%.17g` and Angstrom → bohr is a
division by `Units.bohr_to_angstrom`, not a multiplication by its
reciprocal. A perturbation of `+5e-9` bohr on one coordinate survives
`qp_edit -c` exactly. The pre-fix build printed Angstrom at `%16.8f`
(a 1.9e-8 bohr grid) and parsed that back; the measured drift is in
`transcripts/pre-fix-quantisation/probe.log`. Nuclear charges and ghost
atoms follow: `Charge.to_string` is `%.17g` (no `+%f`, no `to_int`).

**C6 — idempotence.** A second `qp_edit -c` on the result of the first
creates nothing new and changes nothing.

**C7 — malformed input is still reported.** Parsing is not skipped by
`apply=false`. There are two measured failure classes: (a) a sexplib
conversion failure (a value the fancy type rejects) prints
`Error:  (...)` on stderr (two-space prefix, value snippet; the
type-path wording is approximated) and the exit status stays **0** —
there is no `Info:` line; (b) an escaped `Invalid_argument` (e.g. an
out-of-range index) prints an `Info:` line, the retry raises again, and
the driver banners and exits **1**. Not exercised by this corpus (all
fixtures are valid); the Rust suite pins it (C7-A/C7-B, prefix and exit
status only). The original wording here claimed `Info: Read error in
<section>` + exit 0 for any bad section: wrong (see CONTRACTS §6).

## How the corpus gates this

* `transcripts/<fixture>/{stdout,stderr,exit}` — the reference build's
  streams and status. C1 and C7 are checked here (C7 vacuously).
* `contracts/check-mutated.json` — the C3/C4 file list, per fixture and
  per contract (`post-fix` = target, `pre-fix` = reference build).
* `transcripts/<fixture>/added/` — the C3 payloads.
* `verify --contract post-fix` (default) and `verify --contract pre-fix`.

C2, C5 and C6 are not observable in these transcripts because the
pre-fix build's rewrite is content-identical and the fixtures are at
their fixed point; the same holds for C4's existing-pair clauses. They
are covered by `transcripts/pre-fix-quantisation/` (C5), and by the Rust
port's own unit tests (a perturbed coordinate, a second run, a consistent
and a stale pair, the read-only branch).
