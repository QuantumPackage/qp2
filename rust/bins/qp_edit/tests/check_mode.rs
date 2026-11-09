//! Contract tests C1–C7 for `qp_edit -c` (`tests/golden/contracts/check-mode.md`).
//!
//! Every test copies a golden fixture into a temporary directory and runs the
//! real binary with `QP_EDIT_CONFIG` pointing at the generated artifact; the
//! checked-in fixtures are never opened for writing.  The gate that runs the
//! same binary over the same corpus is
//! `tools/golden-harness/golden_harness.py verify --only check --contract
//! post-fix`.
//!
//! The tests are named after the contract rule they pin, so a failure names
//! the rule that regressed.

mod common;

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Read;
use std::path::Path;
use std::process::Output;

use common::{assert_silent_success, command, fixture, qp_edit};
use flate2::read::GzDecoder;
use qp_ezfio::{Array, EzFio};
use qp_input::Engine;

/// The six fixtures of the check-mode corpus, in the order of
/// `tests/golden/contracts/check-mutated.json`.
const CORPUS: [&str; 6] = ["h2", "h3", "c2h2", "b2", "h2g", "hcn"];

/// The seven files `-c` must create for `hcn` (contract C3), sorted as a
/// `BTreeMap` walks them.  `determinants/psi_{det,coef}_qp_edit.gz` are
/// deliberately absent (contract C4).
const HCN_CREATED: [&str; 7] = [
    "determinants/bit_kind",
    "determinants/expected_s2",
    "determinants/mo_label",
    "determinants/n_int",
    "determinants/psi_coef.gz",
    "determinants/psi_det.gz",
    "determinants/state_average_weight.gz",
];

/// The four `nuclei` files the check-mode nuclei write-back (w12) re-emits on
/// every `-c`: the section handler's `write` always rewrites all four, even
/// when their content is unchanged.  The three `.gz` containers are re-encoded
/// with the port's own gzip (`level 6`, `mtime 0`) rather than the committed
/// fixtures' normalisation (`level 9`, `mtime 0`), so their raw bytes differ
/// even for an identical payload — they are compared on the decompressed
/// payload, not the raw container.
const NUCLEI_WRITTEN: [&str; 4] = [
    "nuclei/nucl_num",
    "nuclei/nucl_label.gz",
    "nuclei/nucl_charge.gz",
    "nuclei/nucl_coord.gz",
];

/// Load the generated configuration.
fn engine() -> Engine {
    let directory = common::config_dir();
    match qp_input::loader::load(&directory) {
        Ok(config) => Engine::new(config),
        Err(error) => panic!("loading {}: {error:#}", directory.display()),
    }
}

/// `qp_edit -c ezfio`, the gate's command.
fn check(root: &Path) -> Output {
    qp_edit(root, &["-c", "ezfio"])
}

/// Every regular file under `root`, raw bytes, keyed by `/`-separated path.
///
/// The root-level `work/` directory is skipped: `qp_edit` creates it (and
/// clears `qpstop`/`qpkill` there), and the harness excludes it from the EZFIO
/// state it compares (`SNAPSHOT_SKIP_DIRS` in `golden_harness.py`).
fn snapshot(root: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut files = BTreeMap::new();
    collect(root, Path::new(""), &mut files);
    files
}

/// Recursive worker of [`snapshot`]; `prefix` is relative to `root`.
fn collect(root: &Path, prefix: &Path, files: &mut BTreeMap<String, Vec<u8>>) {
    let directory = root.join(prefix);
    let mut entries: Vec<_> = fs::read_dir(&directory)
        .unwrap_or_else(|error| panic!("reading {}: {error}", directory.display()))
        .map(|entry| entry.expect("read one directory entry"))
        .collect();
    entries.sort_by_key(|entry| entry.file_name());
    for entry in entries {
        let relative = prefix.join(entry.file_name());
        if entry.file_type().expect("entry type").is_dir() {
            if prefix.as_os_str().is_empty() && relative == Path::new("work") {
                continue;
            }
            collect(root, &relative, files);
        } else {
            let key = relative.to_string_lossy().replace('\\', "/");
            let bytes = fs::read(entry.path())
                .unwrap_or_else(|error| panic!("reading {}: {error}", entry.path().display()));
            files.insert(key, bytes);
        }
    }
}

/// The decompressed payload of one file: what
/// `transcripts/<fixture>/added/` stores (suffix stripped, gzip undone).
fn payload(path: &Path) -> Vec<u8> {
    let raw = fs::read(path).unwrap_or_else(|error| panic!("reading {}: {error}", path.display()));
    decompress(&raw, &path.display().to_string())
}

/// Undo gzip when `raw` carries the gzip magic; pass everything else through.
/// `payload` is the on-disk form of this.  The w12 nuclei write-back rewrites
/// the `.gz` containers with the port's own gzip even when the payload is
/// unchanged, so the payload — not the raw bytes — is the comparison unit.
fn decompress(raw: &[u8], label: &str) -> Vec<u8> {
    if raw.starts_with(&[0x1f, 0x8b]) {
        let mut decoded = Vec::new();
        GzDecoder::new(raw)
            .read_to_end(&mut decoded)
            .unwrap_or_else(|error| panic!("decompressing {label}: {error}"));
        decoded
    } else {
        raw.to_vec()
    }
}

/// The paths present after a run but not before it, sorted.
fn created<'a>(
    before: &BTreeMap<String, Vec<u8>>,
    after: &'a BTreeMap<String, Vec<u8>>,
) -> Vec<&'a str> {
    after
        .keys()
        .filter(|path| !before.contains_key(*path))
        .map(String::as_str)
        .collect()
}

/// [`snapshot`] with gzip undone: every regular file under `root`, keyed by
/// `/`-separated path, mapped to its decompressed payload.  The comparison
/// unit for the w12 nuclei write-back, whose `.gz` containers are re-encoded
/// with the port's own gzip even when the payload is unchanged.
fn payload_snapshot(root: &Path) -> BTreeMap<String, Vec<u8>> {
    snapshot(root)
        .into_iter()
        .map(|(path, bytes)| (path.clone(), decompress(&bytes, &path)))
        .collect()
}

/// Assert that `after` has exactly the keys of `before`, each with the same
/// bytes: nothing was created, deleted or rewritten.
///
/// The first offending path is reported, not the two maps: a whole-map
/// `assert_eq!` would print every payload on failure.
fn assert_unchanged(
    before: &BTreeMap<String, Vec<u8>>,
    after: &BTreeMap<String, Vec<u8>>,
    what: &str,
) {
    let added = created(before, after);
    assert!(added.is_empty(), "{what}: created {added:?}");
    for (path, bytes) in before {
        match after.get(path) {
            None => panic!("{what}: {path} was deleted"),
            Some(now) => assert_eq!(now, bytes, "{what}: {path} was rewritten"),
        }
    }
}

#[test]
fn c1_check_exits_zero_with_empty_streams() {
    for name in CORPUS {
        let directory = fixture(name);
        let output = check(directory.path());
        assert_silent_success(&output, name);
    }
}

#[test]
fn c2_no_pre_existing_file_is_rewritten() {
    for name in CORPUS {
        let directory = fixture(name);
        let root = directory.path().join("ezfio");
        // Every file except the w12 nuclei write-back set is compared
        // raw-byte strict: the post-fix tool must not touch it at all, so even
        // a content-identical rewrite (the pre-fix MD5 stamps, D4) fails here.
        //
        // The four nuclei files are the w12 exception: check mode now applies
        // the nuclei section (the per-section `applies_in_check` flag), and the
        // handler rewrites all four on every `-c` even when unchanged.  Their
        // `.gz` containers are re-encoded with the port's own gzip, so the raw
        // bytes differ from the committed fixtures even for an identical
        // payload — they are compared on the decompressed payload instead.
        //
        // Provenance of the `c2h2`/`b2` pins (measured 2026-11-07, review
        // follow-up): the corpus was generated by the reference pre-15e725aa
        // build, whose angstrom→bohr conversion multiplied by the reciprocal
        // (`fl(a * (1 / K))`); the port divides (`fl(a / K)`) since w9
        // (15e725aa parity).  For 1.6644 Å (c2h2) and 3.0 Å (b2) the two
        // pipelines differ by 1 ulp at the 16th significant digit, exposed
        // by the 16-sig-digit `%24.15e` `nucl_coord` text; every other
        // corpus coordinate renders identically in both.  The w12 nuclei
        // write-back runs the divide pipeline, so the first `-c` on a
        // committed fixture normalises the stored payloads of `c2h2`/`b2` by
        // that 1 ulp: display-invisible (the 8-decimal Å strings are
        // unchanged), one-time, idempotent.  The pin is
        // `contracts/check-mutated.json` (`post-fix-changed`) plus the
        // decompressed payloads in `transcripts/<name>/changed/`, which this
        // test reads as the single source of truth.  The fixtures are
        // deliberately NOT regenerated: the CI golden step runs the
        // reference build against the corpus, and its multiply pipeline is
        // the reference's fixed point.
        let before = snapshot(&root);
        assert_silent_success(&check(directory.path()), name);
        let after = snapshot(&root);
        for (path, bytes) in &before {
            match after.get(path) {
                None => panic!("{name}: {path} was deleted"),
                Some(now) => {
                    if path == "nuclei/nucl_coord.gz" && (name == "c2h2" || name == "b2") {
                        // The pinned one-time normalisation: the stored
                        // payload becomes exactly the pinned expected
                        // payload, and it is a real change (a pin that
                        // already matches the fixture is stale, not a
                        // normalisation).
                        let pinned = Path::new(env!("CARGO_MANIFEST_DIR"))
                            .join("../../tests/golden/transcripts")
                            .join(name)
                            .join("changed/nuclei/nucl_coord");
                        let want = fs::read(&pinned).unwrap_or_else(|error| {
                            panic!("reading the pinned payload {}: {error}", pinned.display())
                        });
                        assert_ne!(
                            decompress(bytes, &format!("{name} before {path}")),
                            want,
                            "{name}: {path}: the fixture already carries the pinned value; the pin is stale"
                        );
                        assert_eq!(
                            decompress(now, &format!("{name} after {path}")),
                            want,
                            "{name}: {path} payload is not the pinned normalised payload"
                        );
                    } else if NUCLEI_WRITTEN.contains(&path.as_str()) {
                        assert_eq!(
                            decompress(bytes, &format!("{name} before {path}")),
                            decompress(now, &format!("{name} after {path}")),
                            "{name}: {path} payload was changed by check mode"
                        );
                    } else {
                        assert_eq!(now, bytes, "{name}: {path} was rewritten");
                    }
                }
            }
        }
    }
}

#[test]
fn c3_missing_defaults_are_materialised() {
    // The five fixtures without MOs get nothing: their determinants section
    // returns `None` as soon as the two keys `qp_create_ezfio` left behind are
    // present.
    for name in &CORPUS[..5] {
        let directory = fixture(name);
        let root = directory.path().join("ezfio");
        let before = snapshot(&root);
        assert_silent_success(&check(directory.path()), name);
        let after = snapshot(&root);
        let created = created(&before, &after);
        assert!(created.is_empty(), "{name}: check mode created {created:?}");
    }

    // `hcn` has converged MOs, so the determinants getter materialises its
    // whole record.
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    let before = snapshot(&root);
    assert_silent_success(&check(directory.path()), "hcn");
    let after = snapshot(&root);
    assert_eq!(created(&before, &after), HCN_CREATED.to_vec());

    let golden =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/transcripts/hcn/added");
    for path in HCN_CREATED {
        let name = path.strip_suffix(".gz").unwrap_or(path);
        let expected = fs::read(golden.join(name))
            .unwrap_or_else(|error| panic!("reading the golden payload {name}: {error}"));
        let actual = payload(&root.join(path));
        assert_eq!(
            actual, expected,
            "{path} differs from the reference payload"
        );
    }
}

/// Contract C3, audit item: a materialised default comes from the *installed*
/// `$QP_ROOT/data/ezfio_defaults/` tree, not from the configuration's copied
/// `default:` value.
///
/// The generated readers call `Qpackage.get_ezfio_default`, so the value a
/// build writes depends on the data tree it runs against.  The child gets a
/// private root whose file disagrees with the config (15) for a field the
/// fixture has: removing the key makes the getter materialise it.
#[test]
fn c3_b_an_absent_key_takes_the_runtime_default() {
    const SECTION: &str = "scf_utils";
    const FIELD: &str = "max_dim_diis";
    const CONFIG_DEFAULT: i64 = 15;
    const RUNTIME_DEFAULT: i64 = 42;

    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    fs::remove_file(root.join(SECTION).join(FIELD)).expect("remove the key");

    let qp_root = tempfile::tempdir().expect("create a fake QP_ROOT");
    let defaults = qp_root.path().join(qp_input::defaults::DEFAULTS_DIR);
    fs::create_dir_all(&defaults).expect("create the defaults directory");
    // One file per section is what the generator writes; the reference scans
    // every file, so the name is not load-bearing.
    fs::write(
        defaults.join("scf_utils.ezfio_interface_default"),
        format!("{SECTION}\n  {FIELD}   {RUNTIME_DEFAULT}\n"),
    )
    .expect("write the runtime default");

    let output = command(directory.path())
        .env("QP_ROOT", qp_root.path())
        .args(["-c", "ezfio"])
        .output()
        .expect("run qp_edit");
    assert_silent_success(&output, "runtime default");

    let ezfio = EzFio::open(&root).expect("reopen the fixture copy");
    let value = ezfio
        .get::<i64>(SECTION, FIELD)
        .expect("read the materialised key");
    assert_eq!(
        value, RUNTIME_DEFAULT,
        "the config default ({CONFIG_DEFAULT}) was used instead of the installed one"
    );

    // The materialised key now exists, so a second run must not touch it
    // (contract C2); this is what makes the back-fill a no-op on a converged
    // tree.
    let bytes = fs::read(root.join(SECTION).join(FIELD)).expect("read the key file");
    assert_silent_success(&check(directory.path()), "second -c");
    assert_eq!(
        fs::read(root.join(SECTION).join(FIELD)).expect("read the key file"),
        bytes,
        "the second -c rewrote the materialised key"
    );
}

/// Contract C3, deviation B8: an installed `$QP_ROOT/data/ezfio_defaults/`
/// tree that is missing a section's file still yields the built-in default.
///
/// `Qpackage.get_ezfio_default` raises `Not_found` after scanning the whole
/// directory (printing `<section>/<field> not found`), and `get` only catches
/// `Sys_error`; the port's `defaults.rs` returns `None` and
/// `materialise_default` falls back to the TOML `default:`.  This is the
/// Ampere shape (54 of 55 files, `symmetry` absent): the section is
/// materialised instead of lost.
#[test]
fn c3_c_a_missing_defaults_file_falls_back_to_the_built_in_default() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    let key = root.join("symmetry/auto_symmetry");
    let original = fs::read(&key).expect("read the fixture's auto_symmetry");
    fs::remove_file(&key).expect("remove auto_symmetry");

    // A QP root whose defaults tree exists but holds no symmetry file: the
    // gate's Ampere shape (54 of 55 files), in its limit.
    let qp_root = tempfile::tempdir().expect("create a fake QP_ROOT");
    fs::create_dir_all(qp_root.path().join(qp_input::defaults::DEFAULTS_DIR))
        .expect("create the defaults directory");

    let output = command(directory.path())
        .env("QP_ROOT", qp_root.path())
        .args(["-c", "ezfio"])
        .output()
        .expect("run qp_edit");
    assert_silent_success(&output, "missing defaults file");

    let ezfio = EzFio::open(&root).expect("reopen the fixture copy");
    assert!(
        ezfio
            .get::<bool>("symmetry", "auto_symmetry")
            .expect("read the materialised auto_symmetry"),
        "the config default of auto_symmetry is true"
    );
    // The fixture's own value came from the reference client: the fallback
    // must produce the same `T\n\n`, not a different spelling of true.
    assert_eq!(
        fs::read(&key).expect("read the materialised auto_symmetry"),
        original,
        "the fallback value does not match the built-in OCaml-layout T"
    );
}

#[test]
fn c4_the_qp_edit_scratch_pair_is_not_written() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    assert_silent_success(&check(directory.path()), "hcn");
    for name in [
        "determinants/psi_det_qp_edit.gz",
        "determinants/psi_coef_qp_edit.gz",
    ] {
        assert!(
            !root.join(name).exists(),
            "{name} belongs to the section writer that check mode skips"
        );
    }
}

/// Contract C4, second clause: a pair that already exists survives `-c`
/// byte-identically.
///
/// The reference build re-runs the section writer in check mode, and the
/// determinants writer refreshes the pair from `psi_det`/`psi_coef` whenever
/// `n_det <= n_det_qp_edit`, so its pair changes.  The port writes nothing.
/// The pair is forged from the working arrays the first `-c` materialises,
/// which is the state an interactive `qp_edit -w` leaves behind.
#[test]
fn c4_b_a_consistent_qp_edit_pair_survives_check() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    assert_silent_success(&check(directory.path()), "setup -c");
    for key in ["psi_det", "psi_coef"] {
        let pair = root.join("determinants").join(format!("{key}_qp_edit.gz"));
        fs::copy(root.join("determinants").join(format!("{key}.gz")), &pair)
            .unwrap_or_else(|error| panic!("forging {}: {error}", pair.display()));
    }

    let before = snapshot(&root);
    assert_silent_success(&check(directory.path()), "hcn with a pair");
    let after = snapshot(&root);
    assert_unchanged(&before, &after, "hcn with a consistent pair");
}

/// Contract C4, third clause: a pair that disagrees with the working arrays
/// is left stale, not repaired and not deleted.
///
/// The reference build re-copies `psi_det`/`psi_coef` into the pair when the
/// edited count equals the stored count; the port writes nothing, so the
/// forged values survive.  This is the accepted residual of DEVIATIONS.md B9.
#[test]
fn c4_c_a_stale_qp_edit_pair_is_not_rewritten() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    assert_silent_success(&check(directory.path()), "setup -c");

    // Same shapes as the working arrays (`n_int = 1`, `bit_kind = 8`,
    // `n_det_qp_edit = 1`, `n_states = 1`), different values.
    {
        let ezfio = EzFio::open(&root).expect("open the fixture copy");
        let det = Array::new(vec![1, 2, 1], vec![3i64, 4]).expect("psi_det_qp_edit shape");
        ezfio
            .set_array("determinants", "psi_det_qp_edit", &det)
            .expect("forge a stale psi_det_qp_edit");
        let coef = Array::new(vec![1, 1], vec![0.5f64]).expect("psi_coef_qp_edit shape");
        ezfio
            .set_array("determinants", "psi_coef_qp_edit", &coef)
            .expect("forge a stale psi_coef_qp_edit");
    }

    let before = snapshot(&root);
    assert_silent_success(&check(directory.path()), "hcn with a stale pair");
    let after = snapshot(&root);
    assert_unchanged(&before, &after, "hcn with a stale pair");

    // The files still carry the forged values, and `psi_det` still carries the
    // HF determinant: `-c` copied neither into the other.
    let ezfio = EzFio::open(&root).expect("reopen the fixture copy");
    let det = ezfio
        .get_array::<i64>("determinants", "psi_det_qp_edit")
        .expect("read the stale pair");
    assert_eq!(det.data(), &[3, 4]);
    let coef = ezfio
        .get_array::<f64>("determinants", "psi_coef_qp_edit")
        .expect("read the stale pair");
    assert_eq!(coef.data(), &[0.5]);
    let working = ezfio
        .get_array::<i64>("determinants", "psi_det")
        .expect("read psi_det");
    assert_eq!(working.data(), &[127, 127]);
}

/// Contract C4, read side: when the edited count differs from the stored one
/// (`n_det_qp_edit != n_det`, the truncated-wave-function state), the section
/// getter reads the pair *instead of* `psi_det`/`psi_coef`.  Check mode must
/// still write nothing: the pair, the working arrays and the counters keep
/// their bytes.
#[test]
fn c4_d_the_scratch_pair_is_read_without_being_touched() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    assert_silent_success(&check(directory.path()), "setup -c");

    // Two stored determinants, one edited: `read` selects `read_only = true`.
    {
        let ezfio = EzFio::open(&root).expect("open the fixture copy");
        ezfio
            .set::<i64>("determinants", "n_det", &2)
            .expect("set n_det");
        let det = Array::new(vec![1, 2, 2], vec![127i64, 127, 191, 127]).expect("psi_det shape");
        ezfio
            .set_array("determinants", "psi_det", &det)
            .expect("write psi_det");
        let coef = Array::new(vec![2, 1], vec![0.9f64, 0.1]).expect("psi_coef shape");
        ezfio
            .set_array("determinants", "psi_coef", &coef)
            .expect("write psi_coef");
        let det = Array::new(vec![1, 2, 1], vec![223i64, 127]).expect("pair shape");
        ezfio
            .set_array("determinants", "psi_det_qp_edit", &det)
            .expect("write psi_det_qp_edit");
        let coef = Array::new(vec![1, 1], vec![0.7f64]).expect("pair shape");
        ezfio
            .set_array("determinants", "psi_coef_qp_edit", &coef)
            .expect("write psi_coef_qp_edit");
    }
    // The getter's branch is chosen by these two counters: `2 != 1` is the
    // truncated state, where `read` switches to the scratch pair.
    {
        let ezfio = EzFio::open(&root).expect("reopen the forged tree");
        assert_eq!(
            ezfio
                .get::<i64>("determinants", "n_det")
                .expect("read n_det"),
            2
        );
        assert_eq!(
            ezfio
                .get::<i64>("determinants", "n_det_qp_edit")
                .expect("read n_det_qp_edit"),
            1
        );
    }

    let before = snapshot(&root);
    assert_silent_success(&check(directory.path()), "hcn with an edited subset");
    let after = snapshot(&root);
    assert_unchanged(&before, &after, "hcn with an edited subset");
}

#[test]
fn c5_a_perturbed_coordinate_is_restored_to_the_grid() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    let coordinate = root.join("nuclei/nucl_coord.gz");

    // The first component of the hcn fixture is exactly 0.0, on the
    // master-display grid.  Snapshot its payload before the perturbation: the
    // grid value is what `-c` must restore.
    let before = payload(&coordinate);

    // Perturb one coordinate by +5e-9 bohr through the same EZFIO client the
    // tool uses.
    {
        let ezfio = EzFio::open(&root).expect("open the fixture copy");
        let mut coord = ezfio
            .get_array::<f64>("nuclei", "nucl_coord")
            .expect("read nucl_coord");
        assert_eq!(coord.data()[0], 0.0, "fixture precondition");
        coord.data_mut()[0] = 5e-9;
        ezfio
            .set_array("nuclei", "nucl_coord", &coord)
            .expect("write the perturbed nucl_coord");
    }

    // The nudge (5e-9 bohr, below half the ≈1.89e-8 bohr grid step) rounds
    // back to the grid value `0.00000000` in the master `%16.8f` display, so
    // the check-mode write-back restores it.
    assert_silent_success(&check(directory.path()), "hcn");

    // Stronger assertion: the file's decompressed payload is back to the
    // pre-perturbation payload — byte equality of the decoded values, covering
    // all nine components, not just the perturbed one.  The raw container
    // bytes are not comparable (the port re-encodes the `.gz` at its own gzip
    // level), so the payload is the right unit.
    assert_eq!(
        payload(&coordinate),
        before,
        "-c did not restore nucl_coord to its pre-perturbation payload"
    );
    let ezfio = EzFio::open(&root).expect("reopen after -c");
    let coord = ezfio
        .get_array::<f64>("nuclei", "nucl_coord")
        .expect("read nucl_coord after -c");
    assert_eq!(
        coord.data()[0],
        0.0,
        "the 5e-9 bohr nudge was not restored to the grid value"
    );
}

/// Contract C5, second clause (w12): a coordinate that is *off* the
/// master-display grid is changed by `-c` to its grid value, and the change
/// touches exactly the nuclei coordinate keys — no other file's payload
/// differs.
#[test]
fn c5_b_an_off_grid_coordinate_is_quantised_to_the_grid() {
    // A fresh, full-precision value: not a multiple of the 1e-8 Å grid step.
    const OFF_GRID: f64 = 0.12345678912345678; // bohr
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");

    // Materialise the determinants defaults (contract C3) so the `-c` under
    // test creates no files: its change must be nuclei-only.
    assert_silent_success(&check(directory.path()), "setup -c");

    // Write the off-grid value into the first coordinate (C at the origin).
    {
        let ezfio = EzFio::open(&root).expect("open the fixture copy");
        let mut coord = ezfio
            .get_array::<f64>("nuclei", "nucl_coord")
            .expect("read nucl_coord");
        assert_eq!(coord.data()[0], 0.0, "fixture precondition");
        coord.data_mut()[0] = OFF_GRID;
        ezfio
            .set_array("nuclei", "nucl_coord", &coord)
            .expect("write the off-grid nucl_coord");
    }

    let before = payload_snapshot(&root);
    assert_silent_success(&check(directory.path()), "hcn");
    let after = payload_snapshot(&root);

    // No file was created or deleted.
    assert_eq!(
        before.keys().collect::<BTreeSet<_>>(),
        after.keys().collect::<BTreeSet<_>>(),
        "the file set changed"
    );
    // Every file except the coordinate's payload is untouched.  The other
    // three nuclei files are re-emitted with identical content, so only
    // `nucl_coord` (the perturbed component) changes.
    let coordinate = "nuclei/nucl_coord.gz";
    for (path, bytes) in &before {
        if path.as_str() == coordinate {
            continue;
        }
        assert_eq!(
            after.get(path),
            Some(bytes),
            "{path} payload changed (should be untouched)"
        );
    }
    // ... and the coordinate's payload actually changed.
    assert_ne!(
        before.get(coordinate),
        after.get(coordinate),
        "nucl_coord was not quantised"
    );

    // The coordinate is now on the master-display grid: a second `-c` is a
    // payload no-op (the reproducibility fixed point).  The assertion is on
    // the stored payload, the level the 16-significant-digit `%24.15e`
    // storage preserves — the in-memory grid value of an arbitrary off-grid
    // input need not round-trip that storage bit-for-bit (a 1-ulp difference
    // that re-quantises to the same 8-decimal display).
    assert_silent_success(&check(directory.path()), "second -c");
    let after_second = payload_snapshot(&root);
    assert_eq!(
        after_second, after,
        "the quantised coordinate is not a fixed point of -c"
    );
}

/// Contract C5, third clause (w12): a coordinate file already at the grid is
/// a fixed point of `-c`.  The first `-c` leaves every payload unchanged
/// (the content fixed point); the second `-c` is a byte-for-byte no-op (the
/// byte fixed point, now that the files carry the port's own gzip encoding).
#[test]
fn c5_c_an_on_grid_coordinate_file_is_a_fixed_point() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");

    // Materialise the determinants defaults (contract C3) so the runs under
    // test create no files.
    assert_silent_success(&check(directory.path()), "setup -c");

    // First -c: hcn is at the grid, so no payload changes.  The `.gz`
    // containers are re-encoded with the port's own gzip, so their raw bytes
    // do change — the payload is the content fixed point.
    let before = payload_snapshot(&root);
    assert_silent_success(&check(directory.path()), "first -c");
    let after_first = payload_snapshot(&root);
    assert_eq!(
        after_first, before,
        "the first -c changed a payload (hcn is at the grid)"
    );

    // Second -c: the files now carry the port's own (deterministic) encoding,
    // so the re-encode is a byte-for-byte no-op.
    let raw_first = snapshot(&root);
    assert_silent_success(&check(directory.path()), "second -c");
    let raw_second = snapshot(&root);
    assert_unchanged(&raw_first, &raw_second, "the second -c");
}

#[test]
fn c6_a_second_check_changes_nothing() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    assert_silent_success(&check(directory.path()), "first -c");
    let first = snapshot(&root);
    assert_silent_success(&check(directory.path()), "second -c");
    let second = snapshot(&root);
    assert_eq!(
        second,
        first,
        "a second -c is not idempotent (created: {:?})",
        created(&first, &second)
    );
}

#[test]
fn c7_a_sexplib_failure_is_reported_and_check_still_succeeds() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");

    // Render the real Electrons section, then drop one of the two fields the
    // generated `of_rst` requires.  Only this section is present in the file:
    // an absent header is a silent no-op, so the other tasks do nothing.
    let engine = engine();
    let text = {
        let ezfio = EzFio::open(&root).expect("open the fixture copy");
        engine.get(&ezfio, "electrons").expect("render electrons")
    };
    let broken = delete_line(&text, "elec_alpha_num =");
    fs::write(directory.path().join("broken.rst"), broken).expect("write broken.rst");

    let output = qp_edit(directory.path(), &["-c", "-r", "broken.rst", "ezfio"]);
    assert_eq!(output.status.code(), Some(0), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    // Wave S2.4: one plain line on stderr, the section is skipped, the check
    // still exits 0 (the sexplib reporter and its snippet are gone).
    assert_eq!(
        String::from_utf8(output.stderr).expect("utf-8 stderr"),
        "qp_edit: parse error in Electrons: record_of_sexp: missing field: elec_alpha_num\n"
    );
}

#[test]
fn c7_b_an_escaped_error_still_fails_the_check() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");

    // `Atom.of_string` rejects the element symbol, and `Nuclei_by_hand.of_rst`
    // lets that exception escape: the driver retries the parse once, silently
    // (the `Info: Read error in <title>` line is gone), and the second escape
    // is wrapped with the section title (wave S2.4: the `=-` banner is gone,
    // the driver prints `qp_edit: {error:#}`).
    let engine = engine();
    let title = engine
        .section("nuclei_by_hand")
        .expect("nuclei_by_hand is a section")
        .title
        .clone();
    let text = {
        let ezfio = EzFio::open(&root).expect("open the fixture copy");
        engine
            .get(&ezfio, "nuclei_by_hand")
            .expect("render Molecule")
    };
    let broken = break_first_atom(&text);
    fs::write(directory.path().join("broken.rst"), broken).expect("write broken.rst");

    let output = qp_edit(directory.path(), &["-c", "-r", "broken.rst", "ezfio"]);
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert_eq!(
        String::from_utf8(output.stderr).expect("utf-8 stderr"),
        format!("qp_edit: {title}: Atom.of_string: Element Xx unknown\n")
    );
}

/// B2-shape regression: a declared shape that no longer matches the file must
/// not fail the read.
///
/// The reference getters (`read_float_array`, `read_int64_array`, …) take rank
/// and extents from the payload and never consult `EZFIO.cfg`'s `size:`; the
/// port used to reject the read with `Error::Shape`, so `qp_edit -c` rendered
/// an empty `Determinants_by_hand` section and the `state_average_weight`
/// repair below never ran.  Growing `n_states` behind the arrays' back is the
/// `cis` state: the getter must read the stale `(1,1)` coefficients and repair
/// the weights to `n_states` uniform values, and check mode must leave every
/// other file alone.
#[test]
fn a_stale_declared_shape_is_read_and_repaired() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    assert_silent_success(&check(directory.path()), "setup -c");

    // Scalars are one 20-column padded integer per file, as the reference
    // `Ezfio.set_int` writes them.
    fs::write(
        root.join("determinants").join("n_states"),
        format!("{:>20}\n", 3),
    )
    .expect("grow n_states");

    let before = snapshot(&root);
    assert_silent_success(&check(directory.path()), "-c on a stale declared shape");
    let after = snapshot(&root);

    // Only the getter's repair may write: `read_state_average_weight` sees one
    // stored weight for three states and re-materialises a uniform vector.
    let created = created(&before, &after);
    assert!(
        created.is_empty(),
        "the repair created {created:?} instead of rewriting only the weights"
    );
    for (path, bytes) in &before {
        if path == "determinants/state_average_weight.gz" {
            continue;
        }
        match after.get(path) {
            None => panic!("{path} was deleted"),
            Some(now) => assert_eq!(now, bytes, "{path} was rewritten"),
        }
    }
    assert_ne!(
        after["determinants/state_average_weight.gz"],
        before["determinants/state_average_weight.gz"],
        "the stale weights were not repaired"
    );

    let ezfio = EzFio::open(&root).expect("reopen the fixture copy");
    let weights = ezfio
        .get_array::<f64>("determinants", "state_average_weight")
        .expect("read the repaired weights");
    assert_eq!(weights.dims(), &[3]);
    // `1.0 / 3.0` survives the `%24.15e` payload exactly (16 significant
    // digits round-trip this double).
    assert_eq!(weights.data(), &[1.0 / 3.0; 3]);

    // Check mode runs no writer, so the working array keeps the shape it was
    // materialised with by the first `-c`, and the scratch pair stays absent.
    let coef = ezfio
        .get_array::<f64>("determinants", "psi_coef")
        .expect("read psi_coef");
    assert_eq!(coef.dims(), &[1, 1], "check mode resized psi_coef");
    for name in ["psi_det_qp_edit.gz", "psi_coef_qp_edit.gz"] {
        let pair = root.join("determinants").join(name);
        assert!(!pair.exists(), "-c wrote {name}");
    }

    // The render proves the section parsed rather than falling back to the
    // empty body: three states, one determinant, and `to_rst` pads the
    // coefficients it does not have with `"0."`.
    assert_silent_success(
        &qp_edit(directory.path(), &["-w", "doc.rst", "ezfio"]),
        "-w",
    );
    let document = fs::read_to_string(directory.path().join("doc.rst")).expect("read doc.rst");
    assert!(
        document
            .contains("state_average_weight = (0.333333333333\t0.333333333333\t0.333333333333)"),
        "the repaired weights are missing from the render"
    );
    assert!(
        document.contains("\nDeterminants ::\n\n  1.\t0.\t0.\n"),
        "the determinants section is empty or its coefficient line regressed"
    );
}

/// B2-shape regression, apply path: the three-state section that
/// [`a_stale_declared_shape_is_read_and_repaired`] renders must be writable.
///
/// `-r` parses the padded coefficients back into a `(1, n_states)` array, so
/// the scratch pair is refreshed from `psi_coef`.  The weights do not take the
/// exact `1.0 / 3.0` value on this path: they come back through `to_rst`'s
/// `string_of_float` (`%.12g`), i.e. the twelve-significant-digit render.
#[test]
fn a_three_state_section_round_trips_through_rst() {
    let directory = fixture("hcn");
    let root = directory.path().join("ezfio");
    assert_silent_success(&check(directory.path()), "setup -c");
    fs::write(
        root.join("determinants").join("n_states"),
        format!("{:>20}\n", 3),
    )
    .expect("grow n_states");
    assert_silent_success(&check(directory.path()), "-c on a stale declared shape");

    assert_silent_success(
        &qp_edit(directory.path(), &["-w", "doc.rst", "ezfio"]),
        "-w",
    );
    assert_silent_success(
        &qp_edit(directory.path(), &["-r", "doc.rst", "ezfio"]),
        "-r",
    );

    let ezfio = EzFio::open(&root).expect("reopen after -r");
    let coef = ezfio
        .get_array::<f64>("determinants", "psi_coef")
        .expect("read psi_coef after -r");
    assert_eq!(coef.dims(), &[1, 3]);
    assert_eq!(coef.data(), &[1.0, 0.0, 0.0]);

    let working = root.join("determinants/psi_coef.gz");
    let scratch = root.join("determinants/psi_coef_qp_edit.gz");
    assert_eq!(
        fs::read(&scratch).unwrap_or_else(|error| panic!("reading {}: {error}", scratch.display())),
        fs::read(&working).unwrap_or_else(|error| panic!("reading {}: {error}", working.display())),
        "the writer must refresh psi_coef_qp_edit from psi_coef"
    );

    let weights = ezfio
        .get_array::<f64>("determinants", "state_average_weight")
        .expect("read state_average_weight after -r");
    assert_eq!(weights.dims(), &[3]);
    assert_eq!(
        weights.data(),
        &[0.333333333333f64; 3],
        "the RST round trip must carry the %.12g render, not 1.0 / 3.0"
    );
}

/// Delete the single line of `text` that contains `needle`.
fn delete_line(text: &str, needle: &str) -> String {
    let mut lines: Vec<&str> = text.split('\n').collect();
    let index = lines
        .iter()
        .position(|line| line.contains(needle))
        .unwrap_or_else(|| panic!("no line contains {needle:?} in {text:?}"));
    lines.remove(index);
    lines.join("\n")
}

/// Replace the element symbol of the first atom in a rendered `Molecule`
/// section with an unknown one, keeping the five-token `Atom.of_string` shape.
fn break_first_atom(section: &str) -> String {
    const MARKER: &str = "\nNuclear coordinates in xyz format (Angstroms) ::\n\n";
    let start = section
        .find(MARKER)
        .unwrap_or_else(|| panic!("no xyz marker in {section:?}"))
        + MARKER.len();
    let (head, rest) = section.split_at(start);
    let mut lines: Vec<&str> = rest.split('\n').collect();
    // `to_rst` emits: atom count, the "  " separator line, then the atoms.
    assert!(lines.len() > 2, "short xyz block: {rest:?}");
    let atom = lines[2].trim_start();
    assert!(atom.starts_with("C "), "unexpected first atom: {atom:?}");
    lines[2] = "  Xx     6 0 0 0";
    format!("{head}{}", lines.join("\n"))
}
