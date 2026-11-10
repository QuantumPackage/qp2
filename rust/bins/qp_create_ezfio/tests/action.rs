//! End-to-end tests for the action's *input* half: geometry, basis table,
//! pseudopotential table and the per-nucleus basis reads (wave S2: the
//! shared `ToolError`).
//!
//! These run the real binary, so they cover what the module unit tests
//! cannot: the phase order (a geometry failure must win over a basis
//! failure), the exit code (2 for every failure), and the fact that a
//! valid run reaches the EZFIO write. Every failure case keeps the
//! reference's scenario (same argv, same fixture mutation), asserts the
//! exit code, and pins the one stable `qp_create_ezfio: <message>` line on
//! standard error — the message is named in the test. No exact-stderr
//! multi-line pins.
//!
//! `QP_EDIT=true` is set even though the writer's trailing editor is the
//! `true` seam: the run must not depend on a host `qp_edit` (see
//! `NOTES.md`).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use tempfile::TempDir;

/// The hydrogen basis fixture: one `S` shell with one primitive.
const H_BASIS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/h.basis");
/// Two hydrogens 0.7 angstrom apart, XYZ form.
const H2_XYZ: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/h2.xyz");
/// The same molecule as a Gaussian Z-matrix.
const H2_ZMT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/h2.zmt");
/// A one-element pseudopotential library.
const H_PSEUDO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/h.pseudo");

/// The tool name of the one-line standard-error error.
const TOOL: &str = "qp_create_ezfio";

/// An OCaml-layout int scalar: `%20d\n` (20 columns, right-aligned).
fn int_scalar(value: i64) -> String {
    format!("{value:>20}\n")
}

/// The text of one scalar file under the created EZFIO directory.
fn scalar(ezfio: &Path, group: &str, key: &str) -> String {
    fs::read_to_string(ezfio.join(group).join(key))
        .unwrap_or_else(|error| panic!("read {group}/{key}: {error}"))
}

/// Run the binary with `QP_ROOT = root`, from a scratch working directory.
fn run(work: &Path, root: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_qp_create_ezfio"))
        .env("QP_ROOT", root)
        .env("QP_EDIT", "true")
        .current_dir(work)
        .args(args)
        .output()
        .expect("the compiled test binary must be runnable")
}

/// Decode stderr, which is ASCII in every case exercised here.
fn stderr(output: &Output) -> String {
    String::from_utf8(output.stderr.clone()).expect("stderr is valid UTF-8")
}

/// A scratch directory and its path as a string, for `&[&str]` arguments.
fn workdir() -> (TempDir, String) {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir
        .path()
        .to_str()
        .expect("tempdir paths are UTF-8")
        .to_owned();
    (dir, path)
}

/// Write `content` under `name` in `work` and return the absolute path.
fn write(work: &Path, name: &str, content: &str) -> PathBuf {
    let path = work.join(name);
    fs::write(&path, content).expect("write fixture");
    path
}

/// `path` as a `&str`, for `&[&str]` arguments.
fn text(path: &Path) -> &str {
    path.to_str().expect("scratch paths are UTF-8")
}

/// Every failure: exit 2, empty stdout, and the single line
/// `qp_create_ezfio: <message>` on standard error.
fn assert_fatal(output: &Output, message: &str) {
    assert_eq!(output.status.code(), Some(2), "stderr: {}", stderr(output));
    assert!(
        output.stdout.is_empty(),
        "stdout: {}",
        String::from_utf8_lossy(&output.stdout)
    );
    assert_eq!(
        stderr(output),
        format!("{TOOL}: {message}\n"),
        "the one line"
    );
}

#[test]
fn a_short_xyz_fails_on_the_missing_atom_lines() {
    // Three announced atoms, two printed: the XYZ reader's missing-lines
    // failure, which the reference did not catch.
    let (work, cwd) = workdir();
    let input = write(work.path(), "short.xyz", "3\nshort\nH 0 0 0\nH 0 0 0.7\n");
    let output = run(Path::new(&cwd), work.path(), &["-b", H_BASIS, text(&input)]);
    assert_fatal(
        &output,
        &format!(
            "too few atom lines in {}: 2 lines for 3 announced atoms",
            input.display()
        ),
    );
}

#[test]
fn an_unknown_element_is_reported() {
    // `zzz` is not an atom count, so the file is re-read as a Z-matrix and
    // the first line is then an element name.
    let (work, cwd) = workdir();
    let input = write(work.path(), "zzz.xyz", "zzz\nH 0 0 0\n");
    let output = run(Path::new(&cwd), work.path(), &["-b", H_BASIS, text(&input)]);
    assert_fatal(
        &output,
        &format!(
            "cannot parse geometry {}: Element Zzz unknown",
            input.display()
        ),
    );
}

#[test]
fn a_zero_nucleus_index_is_rejected() {
    // The auxiliary key is parsed while the basis table is built, after the
    // geometry and before any file is read.
    let (work, cwd) = workdir();
    let spec = format!("0,H:{H_BASIS}");
    let output = run(Path::new(&cwd), work.path(), &["-b", &spec, H2_XYZ]);
    assert_fatal(&output, "nucleus index 0 must be at least 1");
}

#[test]
fn a_nucleus_index_above_ten_thousand_warns_on_stderr() {
    // The `More than 10_000 nuclei` warning moved to standard error (wave
    // S2) and is followed by the out-of-range failure.
    let (work, cwd) = workdir();
    let spec = format!("10001,H:{H_BASIS}");
    let output = run(Path::new(&cwd), work.path(), &["-b", &spec, H2_XYZ]);
    assert_eq!(output.status.code(), Some(2));
    assert!(
        output.stdout.is_empty(),
        "the warning no longer goes to stdout"
    );
    assert_eq!(
        stderr(&output),
        format!(
            "More than 10_000 nuclei\n{TOOL}: nucleus index 10001 is above the maximum of 10000\n"
        )
    );
}

#[test]
fn a_missing_basis_is_reported() {
    // `$QP_ROOT` is the empty scratch directory, so neither the local nor
    // the database lookup succeeds.
    let (work, cwd) = workdir();
    let output = run(Path::new(&cwd), work.path(), &["-b", "nosuch", H2_XYZ]);
    assert_fatal(&output, "Basis nosuch not found");
}

#[test]
fn an_empty_basis_name_opens_the_database_directory() {
    // `--basis=` gives the name `""`. `$QP_ROOT/data/basis/` exists, so the
    // directory is opened; the first read fails with `Is a directory`,
    // which attempt 2 swallows.
    let (work, cwd) = workdir();
    fs::create_dir_all(work.path().join("data").join("basis")).expect("create data dir");
    // A local copy of the input keeps the default output name inside the
    // scratch directory.
    let input = write(
        work.path(),
        "h2.xyz",
        fs::read_to_string(H2_XYZ).as_deref().unwrap(),
    );
    let output = run(Path::new(&cwd), work.path(), &["--basis=", text(&input)]);
    assert_fatal(&output, "Basis not found for atom 1 (H)");
}

#[test]
fn a_missing_pseudo_file_is_reported() {
    // The pseudo table is built after the basis table, so a valid basis is
    // required to reach this message (note the trailing dot).
    let (work, cwd) = workdir();
    let output = run(
        Path::new(&cwd),
        work.path(),
        &["-b", H_BASIS, "-p", "nosuch", H2_XYZ],
    );
    assert_fatal(&output, "Pseudo file nosuch not found.");
}

#[test]
fn a_valid_xyz_run_writes_the_ezfio() {
    let (work, cwd) = workdir();
    let out = work.path().join("h2.ezfio");
    let output = run(
        Path::new(&cwd),
        work.path(),
        &["-b", H_BASIS, "-o", text(&out), H2_XYZ],
    );
    assert_eq!(output.status.code(), Some(0), "stderr: {}", stderr(&output));
    assert_eq!(
        output.stdout,
        format!("{}\n", out.display()).into_bytes(),
        "the directory goes to stdout"
    );
    assert!(output.stderr.is_empty());
    // h2 with the one-S-shell fixture: 2 nuclei, 2 electrons (one per
    // spin), 2 AOs. The `%20d` int layout and the 3-byte logical scalars
    // are the reference's Ocaml writer forms.
    assert_eq!(scalar(&out, "nuclei", "nucl_num"), int_scalar(2));
    assert_eq!(scalar(&out, "electrons", "elec_alpha_num"), int_scalar(1));
    assert_eq!(scalar(&out, "electrons", "elec_beta_num"), int_scalar(1));
    assert_eq!(scalar(&out, "ao_basis", "ao_num"), int_scalar(2));
    assert_eq!(scalar(&out, "basis", "ao_normalized"), "T\n\n");
    assert_eq!(scalar(&out, "ao_basis", "ao_cartesian"), "F\n\n");
    assert_eq!(scalar(&out, "basis", "typ"), "Gaussian\n");
    // The md5 stamp is 32 lowercase hex digits; its exact value is pinned
    // by the 27-row `basis-md5.jsonl` cross-check, not here.
    let md5 = scalar(&out, "ao_basis", "ao_md5");
    assert_eq!(md5.len(), 33, "ao_md5 is 32 hex digits + newline: {md5:?}");
    let hex = md5.trim_end();
    assert!(hex.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f')));
    // The `set_file` stamps (ezfio.ml:431-434): VOLATILE in the golden
    // harness, but `creation` must exist — the Fortran consumers of a
    // created tree (`scf`, ezfio.irp.f:30) probe for it and stop 1
    // without it.
    let creation = scalar(&out, "ezfio", "creation");
    assert!(
        creation.len() > 1 && creation.ends_with('\n'),
        "creation: {creation:?}"
    );
    // `echo $USER`: empty line when the variable is unset.
    assert_eq!(
        scalar(&out, "ezfio", "user"),
        format!("{}\n", std::env::var("USER").unwrap_or_default())
    );
    assert_eq!(scalar(&out, "ezfio", "library"), "qp-ezfio\n");
}

#[test]
fn a_valid_zmatrix_run_writes_the_ezfio() {
    let (work, cwd) = workdir();
    let out = work.path().join("h2.ezfio");
    let output = run(
        Path::new(&cwd),
        work.path(),
        &["-b", H_BASIS, "-o", text(&out), H2_ZMT],
    );
    assert_eq!(output.status.code(), Some(0), "stderr: {}", stderr(&output));
    assert_eq!(output.stdout, format!("{}\n", out.display()).into_bytes());
    // The Z-matrix and the XYZ give the same molecule.
    assert_eq!(scalar(&out, "nuclei", "nucl_num"), int_scalar(2));
    assert_eq!(scalar(&out, "ao_basis", "ao_num"), int_scalar(2));
}

#[test]
fn a_valid_pseudo_run_adjusts_the_molecule() {
    // h.pseudo removes one electron per H: the electron counts drop by 2
    // (alpha first), the charges drop to 0, do_pseudo is set, and
    // nucl_charge_remove is one value per nucleus.
    let (work, cwd) = workdir();
    let out = work.path().join("h2.ezfio");
    let output = run(
        Path::new(&cwd),
        work.path(),
        &["-b", H_BASIS, "-p", H_PSEUDO, "-o", text(&out), H2_XYZ],
    );
    assert_eq!(output.status.code(), Some(0), "stderr: {}", stderr(&output));
    assert_eq!(output.stdout, format!("{}\n", out.display()).into_bytes());
    assert_eq!(scalar(&out, "electrons", "elec_alpha_num"), int_scalar(0));
    assert_eq!(scalar(&out, "electrons", "elec_beta_num"), int_scalar(0));
    assert_eq!(scalar(&out, "pseudo", "do_pseudo"), "T\n\n");
    assert_eq!(scalar(&out, "pseudo", "pseudo_klocmax"), int_scalar(0));
    assert_eq!(scalar(&out, "pseudo", "pseudo_lmax"), int_scalar(0));
}

#[test]
fn an_existing_output_directory_wins_over_the_basis_read() {
    // The `Sys.file_exists` check sits between the table builds and the
    // per-nucleus reads: an existing `-o` wins over a per-nucleus basis
    // read failure (h.basis has no helium section, so `He` fails the
    // read — the table itself is fine).
    let (work, cwd) = workdir();
    let out = work.path().join("h2.ezfio");
    fs::create_dir_all(&out).expect("create the colliding directory");
    let input = write(work.path(), "he.xyz", "1\nhe\nHe 0 0 0\n");
    let output = run(
        Path::new(&cwd),
        work.path(),
        &["-b", H_BASIS, "-o", text(&out), text(&input)],
    );
    assert_fatal(&output, &format!("{} already exists", out.display()));
}

#[test]
fn an_unreadable_basis_fails_after_the_output_check() {
    // Same He case, but `-o` absent: the per-nucleus read runs and fails.
    let (work, cwd) = workdir();
    let input = write(work.path(), "he.xyz", "1\nhe\nHe 0 0 0\n");
    let output = run(
        Path::new(&cwd),
        work.path(),
        &[
            "-b",
            H_BASIS,
            "-o",
            text(&work.path().join("he.ezfio")),
            text(&input),
        ],
    );
    assert_fatal(&output, "Basis not found for atom 1 (He)");
}

#[test]
fn a_missing_parent_directory_is_reported() {
    // `-o` passes the exists check (the path is absent), then the single-
    // level `mkdir` of the directory creation fails.
    let (work, cwd) = workdir();
    let out = work.path().join("no").join("parent").join("h2.ezfio");
    let output = run(
        Path::new(&cwd),
        work.path(),
        &["-b", H_BASIS, "-o", text(&out), H2_XYZ],
    );
    assert_fatal(
        &output,
        &format!(
            "cannot create {}: No such file or directory (os error 2)",
            out.display()
        ),
    );
    assert!(!out.exists(), "no partial directory is left behind");
}
