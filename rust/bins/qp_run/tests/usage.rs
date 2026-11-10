//! Command-line handling (`ocaml/qp_run.ml:232-270`).
//!
//! Wave S2.1: `clap` owns the command line. `-h`/`--help` prints its help
//! once and exits **0** (the reference's print-and-continue, with its
//! double-print on a failing argument count, is gone); a wrong number of
//! positionals is a `clap` usage error with its own plain diagnostic,
//! exit 2.

mod common;

use std::os::unix::fs::PermissionsExt;
use std::process::{Command, Stdio};
use std::time::Duration;

use common::{Fixture, Running};

/// No argument at all: the `clap` missing-positionals error, exit 2, no
/// help block and no program list on stdout.
#[test]
fn no_arguments_is_a_usage_error() {
    let fixture = Fixture::new();
    fixture.program("scf", "exit 0");
    let outcome = fixture.run(&[]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert!(outcome.stdout.is_empty(), "{outcome}");
    assert!(
        outcome
            .stderr
            .contains("the following required arguments were not provided"),
        "{outcome}"
    );
    assert!(outcome.stderr.contains("<PROGRAM>"), "{outcome}");
    assert!(outcome.stderr.contains("<EZFIO_DIR>"), "{outcome}");
}

/// `-h` prints the `clap` help once, on standard output, and exits 0.
#[test]
fn help_flag_prints_the_help_once_and_exits_zero() {
    let fixture = Fixture::new();
    fixture.program("scf", "exit 0");
    let outcome = fixture.run(&["-h"]);
    assert_eq!(outcome.code(), 0, "{outcome}");
    assert!(outcome.stderr.is_empty(), "{outcome}");
    assert_eq!(outcome.stdout.matches("Usage:").count(), 1, "{outcome}");
    // `clap`'s option rendering lists the reference's options.
    assert!(outcome.stdout.contains("--slave"), "{outcome}");
    assert!(outcome.stdout.contains("--prefix"), "{outcome}");
    assert!(outcome.stdout.contains("PROGRAM"), "{outcome}");
    assert!(outcome.stdout.contains("EZFIO_DIR"), "{outcome}");
}

/// One positional is not enough: the missing `EZFIO_DIR`, exit 2.
#[test]
fn a_single_argument_is_a_usage_error() {
    let fixture = Fixture::new();
    let outcome = fixture.run(&["scf"]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert!(outcome.stdout.is_empty(), "{outcome}");
    assert!(outcome.stderr.contains("<EZFIO_DIR>"), "{outcome}");
}

/// Three positionals: the extra one, exit 2.
#[test]
fn a_third_argument_is_a_usage_error() {
    let fixture = Fixture::new();
    fixture.program("scf", "exit 0");
    let outcome = fixture.run(&["scf", fixture.ezfio_str(), "extra"]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert!(outcome.stdout.is_empty(), "{outcome}");
    assert!(
        outcome.stderr.contains("unexpected argument 'extra' found"),
        "{outcome}"
    );
}

/// `$QP_ROOT` unset and no `data/executables` above the executable: the
/// reference's `Qpackage.root` message, exit 2. The root is resolved before
/// the command line is parsed (wave S2.1), so a valid command line still
/// shows the root error.
///
/// The binary is copied into a fresh directory because the Cargo target tree
/// sits inside a QP2 worktree, whose `data/executables` would be found by the
/// derivation from `current_exe`.
#[test]
fn missing_qp_root_is_reported() {
    let dir = tempfile::tempdir().expect("tempdir");
    let copy = dir.path().join("qp_run");
    std::fs::copy(common::BIN, &copy).expect("copy the binary");
    let mut permissions = std::fs::metadata(&copy).expect("metadata").permissions();
    permissions.set_mode(0o755);
    std::fs::set_permissions(&copy, permissions).expect("chmod the copy");

    let mut command = Command::new(&copy);
    command
        .args(["scf", "some.ezfio"])
        .env_remove("QP_ROOT")
        .stdin(Stdio::null());
    let outcome = Running::spawn(command, &dir.path().join("logs")).wait(Duration::from_secs(10));
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert!(
        outcome
            .stderr
            .starts_with("qp_run: QP_ROOT environment variable is not set."),
        "{outcome}"
    );
    assert!(
        outcome
            .stderr
            .contains("Please source the quantum_package.rc file."),
        "{outcome}"
    );
}
