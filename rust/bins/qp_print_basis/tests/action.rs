//! End-to-end tests for `qp_print_basis`.
//!
//! The success expectations are captures of the reference build
//! (`~/qp2/ocaml/qp_print_basis.native`, `QP_ROOT=~/qp2`), byte for byte:
//!
//! ```text
//! ~/qp2/ocaml/qp_print_basis.native                      # no argument
//! ~/qp2/ocaml/qp_print_basis.native /nonexistent
//! ~/qp2/ocaml/qp_print_basis.native <copy of tests/golden/fixtures/h2>
//! ~/qp2/ocaml/qp_print_basis.native <copy of tests/golden/fixtures/hcn>
//! ~/qp2/ocaml/qp_print_basis.native <empty directory>
//! ```
//!
//! The golden fixtures are **copied** into a temporary directory before they
//! are printed: the reference read path stamps `ao_basis/ao_md5` and
//! materialises `ao_cartesian` in the tree it opens, and the checked-in
//! fixtures must stay pristine.
//!
//! `tests/data/h2_basis.rst` is the complete 1145-byte stdout of the `h2`
//! run (its basis block plus the newline `print_endline` adds);
//! `tests/data/hcn_basis_prefix.rst` is the first 200 bytes of the `hcn`
//! stdout. The full `hcn` document is byte-compared by the parent's dual-run
//! gate, not here.
//!
//! Wave S2: the error text is Rust-native, so the failure cases assert the
//! exit code and at most one stable substring (named in the test name),
//! never the exact line; the stdout of the frozen data path stays
//! byte-pinned.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use tempfile::TempDir;

/// The `h2` stdout, complete.
const H2_BASIS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/h2_basis.rst");
/// The first 200 bytes of the `hcn` stdout.
const HCN_PREFIX: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/hcn_basis_prefix.rst"
);
/// Length of the complete `hcn` stdout (measured, 62782 bytes).
const HCN_LEN: usize = 62_782;

/// Copy `rust/tests/golden/fixtures/<name>` under `<temporary>/ezfio`.
///
/// The fixture root *is* an EZFIIO directory; the returned handle owns the copy.
fn fixture(name: &str) -> TempDir {
    let source = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/fixtures")
        .join(name);
    let directory = tempfile::tempdir().expect("create a temporary directory");
    copy_tree(&source, &directory.path().join("ezfio"));
    directory
}

/// Recursively copy `from` into `to`.
fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).expect("create the destination directory");
    for entry in fs::read_dir(from).expect("read the fixture directory") {
        let entry = entry.expect("read one fixture entry");
        let target = to.join(entry.file_name());
        if entry.file_type().expect("fixture entry type").is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), &target).expect("copy a fixture file");
        }
    }
}

/// Run the binary with `args`, from a scratch working directory.
fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_qp_print_basis"))
        .args(args)
        .output()
        .expect("the compiled test binary must be runnable")
}

/// Decode stderr, which is ASCII in every case exercised here.
fn stderr(output: &Output) -> String {
    String::from_utf8(output.stderr.clone()).expect("stderr is valid UTF-8")
}

/// Read the file at `path` as bytes.
fn bytes(path: &str) -> Vec<u8> {
    fs::read(path).unwrap_or_else(|error| panic!("read {path}: {error}"))
}

/// The reference dies on `Sys.argv.(1)`; the port says the argument is
/// missing. Stable substring: `argument is missing`.
#[test]
fn no_argument_is_the_missing_argument_error() {
    let output = run(&[]);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert!(
        stderr(&output).contains("argument is missing"),
        "{}",
        stderr(&output)
    );
}

/// The one stable substring is `EZFIIO` (the message names the file form).
#[test]
fn nonexistent_path_is_the_ezfio_file_failure() {
    let output = run(&["/nonexistent"]);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert!(stderr(&output).contains("EZFIIO"), "{}", stderr(&output));
}

/// `Path::exists("")` is false, like the reference's `Sys.file_exists ""`
/// (measured). Stable substring: `EZFIIO`.
#[test]
fn empty_argument_is_the_ezfio_file_failure() {
    let output = run(&[""]);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert!(stderr(&output).contains("EZFIIO"), "{}", stderr(&output));
}

/// An existing path that is not an EZFIIO tree (measured with both a
/// directory and a regular file). Stable substring: `not an EZFIIO
/// directory`.
#[test]
fn existing_directory_without_marker_is_not_an_ezfio() {
    let directory = TempDir::new().expect("tempdir");
    let path = directory.path().to_str().expect("UTF-8 tempdir").to_owned();
    let output = run(&[&path]);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert!(
        stderr(&output).contains("not an EZFIIO directory"),
        "{}",
        stderr(&output)
    );
}

#[test]
fn h2_prints_the_basis_then_fails_on_the_missing_mos() {
    let fixture = fixture("h2");
    let path = fixture.path().join("ezfio");
    let output = run(&[path.to_str().expect("UTF-8 tempdir")]);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    // The basis block is already on stdout when the MO read fails.
    assert_eq!(output.stdout, bytes(H2_BASIS), "h2 basis block");
    // Stable substring: `mo set`.
    assert!(stderr(&output).contains("mo set"), "{}", stderr(&output));
}

#[test]
fn hcn_prints_the_basis_then_the_mos() {
    let fixture = fixture("hcn");
    let path = fixture.path().join("ezfio");
    let output = run(&[path.to_str().expect("UTF-8 tempdir")]);
    assert_eq!(output.status.code(), Some(0), "{output:?}");
    assert!(output.stderr.is_empty(), "{output:?}");
    let prefix = bytes(HCN_PREFIX);
    assert_eq!(output.stdout.len(), HCN_LEN, "hcn stdout length");
    assert!(
        output.stdout.starts_with(&prefix),
        "hcn stdout does not start with the reference prefix"
    );
}

/// `Sys.argv.(1)` is the only argument read; the rest are ignored.
#[test]
fn extra_arguments_are_ignored() {
    let fixture = fixture("hcn");
    let path = fixture.path().join("ezfio");
    let path = path.to_str().expect("UTF-8 tempdir");
    let alone = run(&[path]);
    let with_extras = run(&[path, "extra", "-h", "--nope"]);
    assert_eq!(with_extras.status.code(), Some(0), "{with_extras:?}");
    assert_eq!(with_extras.stdout, alone.stdout);
    assert_eq!(with_extras.stderr, alone.stderr);
}

/// The fixture copy is the one that gets printed, and its AO digest is left
/// alone by this port (decision D4: read paths do not restamp).
#[test]
fn the_printed_tree_is_the_copy_and_stays_unchanged() {
    let fixture = fixture("h2");
    let path: PathBuf = fixture.path().join("ezfio");
    let before = fs::read(path.join("ao_basis").join("ao_md5")).expect("read ao_md5");
    let output = run(&[path.to_str().expect("UTF-8 tempdir")]);
    assert!(!output.stdout.is_empty(), "{output:?}");
    let after = fs::read(path.join("ao_basis").join("ao_md5")).expect("read ao_md5");
    assert_eq!(before, after, "ao_md5 must not be restamped");
}
