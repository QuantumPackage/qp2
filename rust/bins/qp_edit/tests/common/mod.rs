//! Helpers shared by the `qp_edit` integration tests.
//!
//! Every test copies a golden fixture into a temporary directory and runs the
//! real binary with `QP_EDIT_CONFIG` pointing at the generated artifact.  The
//! checked-in fixtures are never opened for writing: the read paths
//! materialise missing defaults (contract C3).
//!
//! Each integration test is its own crate, and each compiles its own copy of
//! this module with a different subset of the helpers in use, so unused-item
//! warnings are silenced here rather than sprinkled over the test files.

#![allow(dead_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The generated configuration directory.
///
/// `$QP_EDIT_CONFIG` wins, exactly like the binary's loader.  Otherwise the
/// worktree layout is assumed: the artifact is gitignored and lives at the QP
/// root, one level above the Cargo workspace (`rust/`).
pub fn config_dir() -> PathBuf {
    if let Some(directory) = std::env::var_os("QP_EDIT_CONFIG") {
        return PathBuf::from(directory);
    }
    let directory = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../data/qp_edit_config");
    assert!(
        directory.join(qp_input::REGISTRY_FILE).is_file(),
        "missing {}: the configuration is a gitignored build artifact, regenerate it \
         with scripts/ezfio_interface/ei_handler.py (see rust/DESIGN-qp_edit-config.md)",
        directory.display()
    );
    directory
}

/// Copy `tests/golden/fixtures/<name>` under `<temporary>/ezfio`.
///
/// The fixture root *is* an EZFIO directory; the returned handle owns the
/// copy and removes it on drop.
pub fn fixture(name: &str) -> tempfile::TempDir {
    let source = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/fixtures")
        .join(name);
    let directory = tempfile::tempdir().expect("create a temporary directory");
    copy_tree(&source, &directory.path().join("ezfio"));
    directory
}

/// Recursively copy `from` into `to`.
pub fn copy_tree(from: &Path, to: &Path) {
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

/// A command running the binary from `root` with the generated configuration.
///
/// Callers add arguments and environment variables (an `$EDITOR`, for the
/// interactive tests) before [`Command::output`].
///
/// `$QP_ROOT` is removed: the binary reads the runtime defaults from
/// `$QP_ROOT/data/ezfio_defaults/`, and the tests must not depend on the
/// developer's installation.  A test that exercises that path sets its own
/// fake root with [`Command::env`].
pub fn command(root: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_qp_edit"));
    command
        .current_dir(root)
        .env_remove("QP_ROOT")
        .env("QP_EDIT_CONFIG", config_dir());
    command
}

/// Run the binary from `root` with the generated configuration.
pub fn qp_edit(root: &Path, args: &[&str]) -> Output {
    command(root).args(args).output().expect("run qp_edit")
}

/// Assert exit 0 and both streams empty (contract C1).
pub fn assert_silent_success(output: &Output, what: &str) {
    assert_eq!(output.status.code(), Some(0), "{what}: {output:?}");
    assert!(
        output.stdout.is_empty(),
        "{what}: stdout is not empty: {}",
        String::from_utf8_lossy(&output.stdout)
    );
    assert!(
        output.stderr.is_empty(),
        "{what}: stderr is not empty: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}
