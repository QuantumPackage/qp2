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
4. read the RST back and parse every section (`set`): a malformed
   section is still detected. Since w12 the nuclei section is the one
   check mode **applies** (its parsed master-format record is written
   back, C2); every other section parses with `apply:false` (nothing
   written);
5. exit 0.

## Rules

**C1 — exit status and streams.** On a valid directory the exit status is
0, stdout is empty, stderr is empty.

**C2 — no rewrite.** No file that existed before the run is modified or
deleted: every such path has the same decompressed payload afterwards.
This includes the `ao_md5`, `ao_extra_md5` and `mo_basis/ao_md5` stamps,
which the pre-fix build rewrote on every read with identical content and
which D4 of the plan forbids the Rust port from writing at all. The
exceptions are (i) the in-place repairs the reference readers themselves
perform: step 2 (an unordered basis is reordered and written back, a
repair, not a round trip), and a stale-value repair inside a getter —
the known one is `read_state_average_weight`, which rewrites
`state_average_weight.gz` when its length no longer equals
`determinants/n_states`; and (ii) since w12 (user decision,
reproducibility) the **nuclei section**, which check mode now writes:
the parsed master-format record is applied and the four `nuclei` files
are re-emitted. The file-list gate is unchanged (no path added or
removed). The nuclei write is content-identical on the four on-grid
fixtures (`h2`, `h3`, `h2g`, `hcn`); on `c2h2` and `b2` it is the
one-time normalisation pinned below. Every corpus fixture is at its
fixed point for the reader repairs, so no repair fires in the corpus
and C2 holds as stated; the repair is pinned by the Rust check-mode
suite (`a_stale_declared_shape_is_read_and_repaired`).

**The c2h2/b2 `nucl_coord` normalisation (pinned; w12 review
2026-11-07).** The committed fixtures were produced by the reference
pre-`15e725aa` build, whose angstrom→bohr conversion multiplied by the
reciprocal (`fl(a * (1/K))`); the port divides by `K` (`fl(a / K)`) —
the `15e725aa` fix, in the port since w9. For 1.6644 Å (c2h2) and
3.0 Å (b2) the two pipelines differ by 1 ulp at the 16th significant
digit, exposed by the 16-sig-digit `%24.15e` storage text (c2h2
`±3.145260163204818e+00` → `±…817e+00`; b2
`5.669178376360522e+00` → `…521e+00`); every other corpus coordinate
renders identically in both pipelines. Since w12 the nuclei write-back
runs the divide pipeline, so the first `-c` on a committed `c2h2`/`b2`
tree normalises that one payload: display-invisible (the 8-decimal Å
strings are unchanged), one-time, idempotent (a second `-c` changes
nothing). This is a documented consequence of the w9 pipeline, first
observable in w12's nuclei write-back — pinned, not papered over:
`contracts/check-mutated.json` lists the two files under
`post-fix-changed`, the expected decompressed payloads are in
`transcripts/<name>/changed/nuclei/nucl_coord`, and `verify` gates
them in both scenarios (`check`, and `create`, whose trailing `-c`
leaves the same bytes behind) as a **one-way allowance**: a toolset
that applies the normalisation must produce exactly those bytes, and
a toolset that does not apply it (the reference OCaml build in CI,
whose multiply pipeline is the corpus's own fixed point) touches
nothing and passes; the negative side — that the port does apply it —
is enforced by the port's own test suite (`check_mode.rs`
`c2_no_pre_existing_file_is_rewritten`, which also asserts the change
is real, so a stale pin or a silent loss of the write-back fails the
cargo gate). The pre-fix contract (reference evidence) is unchanged:
the reference's multiply pipeline is the corpus's fixed point, and its
`-c` touches nothing. The fixtures are deliberately **not**
regenerated: the CI golden step runs the reference build against this
corpus, and regenerating at the port's divide values would fail the
reference's own create self-verify.

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
lists in `contracts/check-mutated.json`. Since w12 the nuclei section is
the **only** section check mode applies; the determinants section is
still skipped, so the scratch pair is still not written.

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

**C5 — coordinate handling in the nuclei display.** The nuclei block
prints the **master form** (restored by the w12 user decision,
2026-11-06, which retracts the `15e725aa` no-quantisation print for the
nuclei print): `%16.8f %16.8f %16.8f` Angstrom (a 1.9e-8 bohr grid) with
the unsigned integer charge (`Charge::to_int`, `%3d`; no `+%f`).
`Point3d::to_string` / `Charge::to_string` (`%.17g`) remain the
non-display form (e.g. `Atom`). Angstrom → bohr is a **division** by
`Units.bohr_to_angstrom`, not a multiplication by its reciprocal (the
`15e725aa` fix, kept). Because the display quantises, a `+5e-9` bohr
perturbation of a grid coordinate is **quantised away** by `qp_edit -c`
— the master display rounds it back to the grid value and the
check-mode write-back (C2) restores it; an off-grid coordinate is
changed to its grid value. The divide form is what the w12 write-back
stores; on the committed `c2h2`/`b2` fixtures that is the one-time
1-ulp normalisation pinned in C2. The pre-fix build's drift is in
`transcripts/pre-fix-quantisation/probe.log`, kept as evidence of the
pre-`15e725aa` reference behaviour.

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
  per contract (`post-fix` = target, `pre-fix` = reference build), plus
  `post-fix-changed`: the existing files the post-fix build may touch
  (the C2 pinned c2h2/b2 normalisation; nothing else).
* `transcripts/<fixture>/added/` — the C3 payloads.
* `transcripts/<fixture>/changed/` — the C2 pinned normalisation
  payloads.
* `verify --contract post-fix` (default) and `verify --contract pre-fix`;
  `verify --only create` gates the port's create output against the
  fixtures with the same pinned exceptions.

C2, C5 and C6 are not observable in the pre-fix transcripts because
the pre-fix build's rewrite is content-identical and the fixtures are
at their fixed point; the same holds for C4's existing-pair clauses.
They are covered by `transcripts/pre-fix-quantisation/` (C5), and by
the Rust port's own unit tests (a perturbed coordinate, a second run,
a consistent and a stale pair, the read-only branch). The C2
normalisation is observable in both places, split by side: `verify`
(one-way allowance) fails any unexpected touch and any pinned touch
whose payload is not the pinned payload, in both scenarios, for any
toolset; the Rust suite's `c2` pins the positive side (the
normalisation actually happens, to exactly the pinned bytes, and a pin
that matches the fixture already fails).
