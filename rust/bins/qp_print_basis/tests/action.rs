//! End-to-end tests for `qp_print_basis`.
//!
//! Every expectation below is a capture of the reference build
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
/// The fixture root *is* an EZFIO directory; the returned handle owns the copy.
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

/// The reference's `Fatal error: exception Failure("…")` line.
fn failure(message: &str) -> String {
    format!("Fatal error: exception Failure(\"{message}\")\n")
}

/// Read the file at `path` as bytes.
fn bytes(path: &str) -> Vec<u8> {
    fs::read(path).unwrap_or_else(|error| panic!("read {path}: {error}"))
}

#[test]
fn no_argument_is_the_ocaml_index_error() {
    let output = run(&[]);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert_eq!(
        stderr(&output),
        "Fatal error: exception Invalid_argument(\"index out of bounds\")\n"
    );
}

#[test]
fn nonexistent_path_is_the_ezfio_file_failure() {
    let output = run(&["/nonexistent"]);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert_eq!(stderr(&output), failure("Error reading EZFIO file"));
}

/// `Sys.file_exists ""` is false, like the reference's (measured).
#[test]
fn empty_argument_is_the_ezfio_file_failure() {
    let output = run(&[""]);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert_eq!(stderr(&output), failure("Error reading EZFIO file"));
}

/// An existing path that is not an EZFIO tree: `Ezfio.set_file` calls
/// `mkdir`, which returns `EEXIST` (measured with both a directory and a
/// regular file).
#[test]
fn existing_directory_without_marker_is_mkdir_eexist() {
    let directory = TempDir::new().expect("tempdir");
    let path = directory.path().to_str().expect("UTF-8 tempdir").to_owned();
    let output = run(&[&path]);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert_eq!(
        stderr(&output),
        format!("Fatal error: exception Unix.Unix_error(Unix.EEXIST, \"mkdir\", \"{path}\")\n")
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
    assert_eq!(stderr(&output), failure("Error reading the mo set"));
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
