//! Command-line handling (`ocaml/qp_run.ml:232-270`).
//!
//! `set_specs` prints the help when `-h` is present and then *continues*; the
//! anonymous-argument check prints the help a second time before failing. So
//! `qp_run` alone shows the help once and exits 2, `qp_run -h` shows it twice
//! and exits 2, and `qp_run -h scf <dir>` shows it once and runs the job.

mod common;

use std::os::unix::fs::PermissionsExt;
use std::process::{Command, Stdio};
use std::time::Duration;

use common::{Fixture, HELP_MARKER, Running};

/// No argument at all: one help block, then `Inconsistent command line`.
#[test]
fn no_arguments_prints_the_help_once() {
    let fixture = Fixture::new();
    fixture.program("scf", "exit 0");
    let outcome = fixture.run(&[]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert_eq!(outcome.stdout.matches(HELP_MARKER).count(), 1, "{outcome}");
    assert!(
        outcome.stdout.contains(" * scf"),
        "the description lists the programs:\n{outcome}"
    );
    assert!(
        outcome.stderr.contains("Inconsistent command line"),
        "{outcome}"
    );
}

/// `-h` with no anonymous argument: the help is printed by `-h` and again by
/// the failing anonymous check, then the tool exits 2.
#[test]
fn help_flag_prints_the_help_twice() {
    let fixture = Fixture::new();
    fixture.program("scf", "exit 0");
    let outcome = fixture.run(&["-h"]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert_eq!(outcome.stdout.matches(HELP_MARKER).count(), 2, "{outcome}");
    assert!(
        outcome.stderr.contains("Inconsistent command line"),
        "{outcome}"
    );
    // `clap`'s option rendering follows the reference description.
    assert!(
        outcome.stdout.contains("--slave") && outcome.stdout.contains("--prefix"),
        "{outcome}"
    );
}

/// `-h` with a complete command line: the help is printed once and the run
/// goes on, which is `set_specs`' documented behaviour.
#[test]
fn help_flag_continues_to_the_run() {
    let fixture = Fixture::new();
    fixture.program("scf", "exit 0");
    let outcome = fixture.run(&["-h", "scf", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 0, "{outcome}");
    assert_eq!(outcome.stdout.matches(HELP_MARKER).count(), 1, "{outcome}");
    assert!(outcome.stdout.contains("Wall time:"), "{outcome}");
    assert!(!outcome.stdout.contains("Inconsistent command line"));
}

/// One anonymous argument is not enough: same path as no argument.
#[test]
fn a_single_argument_is_inconsistent() {
    let fixture = Fixture::new();
    let outcome = fixture.run(&["scf"]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert_eq!(outcome.stdout.matches(HELP_MARKER).count(), 1, "{outcome}");
    assert!(
        outcome.stderr.contains("Inconsistent command line"),
        "{outcome}"
    );
}

/// `$QP_ROOT` unset and no `data/executables` above the executable: the
/// reference's `Qpackage.root` message, exit 2.
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
