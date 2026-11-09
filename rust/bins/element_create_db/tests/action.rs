//! End-to-end tests for `element_create_db`.
//!
//! The reference writes its table into the live tree, so the port is run with
//! a fake `QP_ROOT` (a temporary directory) in every case:
//!
//! ```text
//! mkdir -p $FAKE/data
//! QP_ROOT=$FAKE ~/qp2/ocaml/element_create_db.byte      # exit 0, file below
//! env -u QP_ROOT ~/qp2/ocaml/element_create_db.byte     # the QP_ROOT line
//! QP_ROOT=$EMPTY ~/qp2/ocaml/element_create_db.byte     # the write error
//! ```
//!
//! `tests/data/list_element.txt` is that first capture, complete (2987 bytes,
//! 110 lines); the dual-run gate diffs it against a fresh reference run.
//!
//! Wave S2: the error text is Rust-native; the failure tests assert the exit
//! code and at most one stable substring (named in the test name), never the
//! exact line.

use std::fs;
use std::path::Path;
use std::process::{Command, Output};

use tempfile::TempDir;

/// The complete reference table.
const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/list_element.txt");

/// Run the binary with `QP_ROOT` set to `root` (`None` = unset).
fn run(root: Option<&Path>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_element_create_db"));
    match root {
        Some(root) => command.env("QP_ROOT", root),
        None => command.env_remove("QP_ROOT"),
    };
    command
        .output()
        .expect("the compiled test binary must be runnable")
}

/// Decode stderr, which is ASCII in every case exercised here.
fn stderr(output: &Output) -> String {
    String::from_utf8(output.stderr.clone()).expect("stderr is valid UTF-8")
}

#[test]
fn writes_the_reference_table() {
    let root = TempDir::new().expect("tempdir");
    fs::create_dir(root.path().join("data")).expect("create data");
    let output = run(Some(root.path()));

    assert_eq!(output.status.code(), Some(0), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert!(output.stderr.is_empty(), "{output:?}");

    let path = root.path().join("data").join("list_element.txt");
    let written = fs::read(&path).expect("the table must exist");
    assert_eq!(
        written,
        fs::read(GOLDEN).expect("read the golden table"),
        "the generated table differs from the reference"
    );

    let text = String::from_utf8(written).expect("the table is ASCII");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 110, "one line per charge 0..=109");
    assert_eq!(lines[0], "  0   X Ghost 0.000000");
    assert_eq!(lines[109], "109  Mt Meitnerium 278.000000");
}

/// The QP_ROOT guard fails before anything is opened: exit 2, no file, no
/// stdout; the one stable substring is the variable name.
#[test]
fn unset_qp_root_is_the_guard_failure() {
    let output = run(None);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert!(stderr(&output).contains("QP_ROOT"), "{}", stderr(&output));
}

/// Writing to `$QP_ROOT/data/…` with no `data/` directory: exit 2, and the
/// one stable substring is the standard `No such file or directory` text.
#[test]
fn missing_data_directory_is_a_write_error() {
    let root = TempDir::new().expect("tempdir");
    let output = run(Some(root.path()));
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    assert!(
        stderr(&output).contains("No such file or directory"),
        "{}",
        stderr(&output)
    );
}

/// The reference's `open_out` truncates: a longer pre-existing file is
/// replaced, not overwritten in place.
#[test]
fn rewriting_truncates_the_target() {
    let root = TempDir::new().expect("tempdir");
    fs::create_dir(root.path().join("data")).expect("create data");
    let path = root.path().join("data").join("list_element.txt");
    fs::write(&path, "x".repeat(10_000)).expect("seed the target");
    let output = run(Some(root.path()));
    assert_eq!(output.status.code(), Some(0), "{output:?}");
    assert_eq!(
        fs::metadata(&path).expect("stat").len(),
        fs::metadata(GOLDEN).expect("stat").len()
    );
}
