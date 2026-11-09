//! Interactive-driver tests: `$EDITOR`, the temporary document, the backup.
//!
//! These lie outside the C1–C7 check-mode contract
//! (`tests/golden/contracts/check-mode.md`): they pin the path that renders
//! the document, runs the external editor and reads it back with `apply`,
//! which no contract rule exercises (check mode skips the editor and every
//! write).
//!
//! The editor is a shell script written by each test.  It records the path it
//! was given (to check the temporary file's lifetime) and applies an `awk`
//! program passed through `$QP_EDIT_TEST_AWK`, so the three tests share one
//! script and differ only in the edit.

#![cfg(unix)]

mod common;

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Output;

use qp_ezfio::EzFio;
use qp_types::point3d::units::BOHR_TO_ANGSTROM;

/// The fake `$EDITOR`.
const EDITOR: &str = r#"#!/bin/sh
# Record the document path, then rewrite it with the test's awk program.
printf '%s' "$1" > "$QP_EDIT_TEST_RECORD"
awk "$QP_EDIT_TEST_AWK" "$1" > "$1.tmp" && mv "$1.tmp" "$1"
"#;

/// Leave the document untouched.
const NO_OP: &str = r#"{ print }"#;

/// Turn on periodic boundary conditions and move the first atom by 0.01
/// Angstrom along x (the `hcn` fixture's carbon sits at the origin).
///
/// The match is the master-format all-zero line (w12): one or more spaces
/// between every token, the coordinate fields rendered as `%16.8f`
/// Angstrom.  The replacement line stays ragged on purpose: it pins that a
/// ragged user edit parses inside a master-format document.
const MOVE_AND_SET_PERIODIC: &str = r#"
/is_periodic = false/ { print "     is_periodic = true"; next }
/^  C +6 +0\.00000000 +0\.00000000 +0\.00000000$/ { print "  C      6 0.01 0 0"; next }
{ print }
"#;

/// Move the first atom, then make the generic `Nuclei` section unparsable.
/// `nuclei_by_hand` is the first task and `nuclei` a much later one, so the
/// moved coordinate is written before the failure and has to be rolled back
/// by the backup.
///
/// The match is the master-format all-zero line (w12); the replacement
/// line stays ragged on purpose (see `MOVE_AND_SET_PERIODIC`).
const MOVE_AND_BREAK_NUCLEI: &str = r#"
/^  C +6 +0\.00000000 +0\.00000000 +0\.00000000$/ { print "  C      6 0.01 0 0"; next }
/io_nuclear_repulsion = None/ { print "     io_nuclear_repulsion = Bogus"; next }
{ print }
"#;

/// The eight entries of the `nuclei` group in the `hcn` fixture.
const NUCLEI_ENTRIES: [&str; 8] = [
    ".version",
    "io_nuclear_repulsion",
    "is_periodic",
    "nucl_charge.gz",
    "nucl_coord.gz",
    "nucl_label.gz",
    "nucl_num",
    "point_charges",
];

/// Write [`EDITOR`] under `directory` and return its path.
fn install_editor(directory: &Path) -> PathBuf {
    let path = directory.join("editor.sh");
    fs::write(&path, EDITOR).expect("write the fake editor");
    fs::set_permissions(&path, fs::Permissions::from_mode(0o755))
        .expect("make the fake editor executable");
    path
}

/// Run `qp_edit` interactively on `ezfio` with `program` as the edit.
///
/// Returns the process output and the path the editor was given.
fn run(directory: &Path, program: &str) -> (Output, PathBuf) {
    let editor = install_editor(directory);
    let record = directory.join("editor-record");
    let output = common::command(directory)
        .env("EDITOR", &editor)
        .env("QP_EDIT_TEST_RECORD", &record)
        .env("QP_EDIT_TEST_AWK", program)
        .arg("ezfio")
        .output()
        .expect("run qp_edit");
    let path = fs::read_to_string(&record).expect("the editor did not run");
    (output, PathBuf::from(path))
}

/// A fresh `hcn` fixture, its EZFIO root.
fn fixture() -> (tempfile::TempDir, PathBuf) {
    let directory = common::fixture("hcn");
    let root = directory.path().join("ezfio");
    (directory, root)
}

#[test]
fn an_unchanged_document_keeps_every_nuclei_entry() {
    let (directory, root) = fixture();
    let (output, document) = run(directory.path(), NO_OP);
    common::assert_silent_success(&output, "no-op editor");

    // The document is a real temporary file outside the EZFIO tree, named
    // like the reference's `Filename.temp_file "qp_edit_" ".rst"`, and it is
    // removed when the run ends.
    let name = document
        .file_name()
        .expect("the document has a file name")
        .to_string_lossy();
    assert!(
        name.starts_with("qp_edit_") && name.ends_with(".rst"),
        "unexpected document name {name}"
    );
    assert!(!document.exists(), "the temporary document was not removed");

    // Regression: reading the whole document back with `apply` must not drop
    // the four `nucl_*` entries of the handler-owned group.
    for entry in NUCLEI_ENTRIES {
        assert!(
            root.join("nuclei").join(entry).exists(),
            "nuclei/{entry} disappeared"
        );
    }
    assert_eq!(
        fs::read_to_string(root.join("nuclei/nucl_num")).expect("read nucl_num"),
        "                   3\n"
    );

    // The backup is removed on success; `work/` itself stays behind.
    let work = root.join("work");
    assert!(work.is_dir(), "create_backup did not create work/");
    assert_eq!(
        fs::read_dir(&work).expect("read work/").count(),
        0,
        "work/ is not empty: {:?}",
        fs::read_dir(&work)
            .expect("read work/")
            .map(|entry| entry.expect("entry").file_name())
            .collect::<Vec<_>>()
    );
}

#[test]
fn the_editor_edits_reach_the_ezfio() {
    let (directory, root) = fixture();
    let (output, document) = run(directory.path(), MOVE_AND_SET_PERIODIC);
    common::assert_silent_success(&output, "editor");
    assert!(!document.exists(), "the temporary document was not removed");

    // The generic `Nuclei` section went through `set`. The OCaml-layer
    // logical scalar is 3 bytes (`T\n\n` — CONTRACTS §6, writer fix
    // `d5439a70`); the 2-byte `T\n` is the Fortran form.
    assert_eq!(
        fs::read_to_string(root.join("nuclei/is_periodic")).expect("read is_periodic"),
        "T\n\n"
    );

    // The handler parsed the xyz block and converted Angstrom back to bohr by
    // dividing by `bohr_to_angstrom`.  The EZFIO `%24.15e` layout keeps 15
    // digits after the decimal point, hence the tolerance.
    let ezfio = EzFio::open(&root).expect("open the edited EZFIO");
    let coord = ezfio
        .get_array::<f64>("nuclei", "nucl_coord")
        .expect("read nucl_coord");
    let expected = 0.01 / BOHR_TO_ANGSTROM;
    assert!(
        (coord.data()[0] - expected).abs() < 1e-15,
        "first coordinate is {:?}, expected {expected:?}",
        coord.data()[0]
    );
}

#[test]
fn a_failed_read_back_restores_the_backup() {
    let (directory, root) = fixture();
    let coordinate = root.join("nuclei/nucl_coord.gz");
    let before = fs::read(&coordinate).expect("read the fixture nucl_coord");

    let (output, _) = run(directory.path(), MOVE_AND_BREAK_NUCLEI);
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    assert!(output.stdout.is_empty(), "{output:?}");
    // Wave S2.4: the domain violation (`Bogus` is not a `Disk_access` choice)
    // escapes the parse, the silent retry escapes again, and the driver
    // prints the error wrapped with the section title — no `Info:` line, no
    // `=-` banner.
    let stderr = String::from_utf8(output.stderr).expect("utf-8 stderr");
    assert_eq!(
        stderr,
        "qp_edit: Nuclei: value \"Bogus\" is outside the domain of type \"Disk_access\" (one of Read|Write|None)\n"
    );

    // `nuclei_by_hand` ran before `nuclei` and wrote the moved coordinate;
    // the failed run must bring the pre-run bytes back.
    assert_eq!(
        fs::read(&coordinate).expect("read nucl_coord after the failed run"),
        before,
        "the backup was not restored"
    );
    for entry in NUCLEI_ENTRIES {
        assert!(
            root.join("nuclei").join(entry).exists(),
            "nuclei/{entry} disappeared"
        );
    }
    assert!(
        !root.join("work/backup.tar").exists(),
        "restore_backup left the archive behind"
    );
}
