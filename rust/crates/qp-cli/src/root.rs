//! `$QP_ROOT` resolution (`Qpackage.root`, `ocaml/Qpackage.ml:6`).
//!
//! The reference reads the `QP_ROOT` environment variable and fails when it is
//! unset. The Rust tools keep that variable as the first source, but fall back
//! to deriving the root from [`std::env::current_exe`] so that a copied binary
//! works without a wrapper script and without a path baked at compile time
//! (issue #224). The derived root is the first ancestor of the executable's
//! directory that contains `data/executables`, which covers both the OCaml
//! layout (`$QP_ROOT/ocaml/_build/qp_run.native`, `$QP_ROOT/ocaml/qp_run`) and
//! the `$QP_ROOT/bin` layout.
//!
//! No path from the build machine is ever embedded in the binary: the only
//! absolute path in this module is the one passed to
//! [`std::env::current_exe`].

use std::path::{Path, PathBuf};

use crate::Result;

/// Directory levels searched above the executable, inclusive of its own
/// directory.
///
/// Four covers `$QP_ROOT/bin`, `$QP_ROOT/ocaml` and
/// `$QP_ROOT/ocaml/_build` with one level to spare.
pub const SEARCH_DEPTH: usize = 4;

/// The directory containing `data/executables`, if any ancestor of `exe` is a
/// QP2 root.
///
/// `exe` may be a symlink target or the path the loader saw; the search only
/// uses its parent chain.
#[must_use]
pub fn derive_root(exe: &Path) -> Option<PathBuf> {
    let directory = exe.parent()?;
    directory
        .ancestors()
        .take(SEARCH_DEPTH)
        .find(|candidate| candidate.join("data").join("executables").is_file())
        .map(Path::to_path_buf)
}

/// Resolve the root from an environment value and an executable path.
///
/// An empty `QP_ROOT` is treated as unset (the reference would use the empty
/// string and then fail to open `data/executables`).
///
/// # Errors
///
/// [`crate::Error::Root`] when the variable is unset and no ancestor of `exe`
/// is a QP2 root.
pub fn resolve(env: Option<&str>, exe: &Path) -> Result<PathBuf> {
    match env.map(str::trim).filter(|value| !value.is_empty()) {
        Some(value) => Ok(PathBuf::from(value)),
        None => derive_root(exe).ok_or(crate::Error::Root),
    }
}

/// The QP2 root of this process: `$QP_ROOT`, else the ancestor of the running
/// executable that contains `data/executables`.
///
/// # Errors
///
/// [`crate::Error::Root`] when neither source works; its message is the
/// reference's, which names the environment variable because setting it is
/// the fix in that case. [`crate::Error::Io`] when
/// [`std::env::current_exe`] fails.
pub fn root() -> Result<PathBuf> {
    let env = std::env::var("QP_ROOT").ok();
    let env = env.as_deref();
    if env.is_some_and(|value| !value.trim().is_empty()) {
        return resolve(env, Path::new(""));
    }
    let exe = std::env::current_exe()?;
    resolve(env, &exe)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build `<tmp>/root/{data/executables,ocaml/_build/qp_run.native}`.
    fn fake_root() -> (tempfile::TempDir, PathBuf, PathBuf) {
        let dir = tempfile::tempdir().expect("tempdir");
        let root = dir.path().join("root");
        std::fs::create_dir_all(root.join("data")).expect("data");
        std::fs::create_dir_all(root.join("ocaml").join("_build")).expect("ocaml");
        std::fs::write(root.join("data").join("executables"), "scf /x/scf\n").expect("write");
        let exe = root.join("ocaml").join("_build").join("qp_run.native");
        (dir, root, exe)
    }

    /// Both reference layouts resolve: `ocaml/_build/` (3 levels up) and
    /// `ocaml/` (2 levels up).
    #[test]
    fn derive_from_both_layouts() {
        let (_dir, root, exe) = fake_root();
        assert_eq!(derive_root(&exe), Some(root.clone()));
        assert_eq!(
            derive_root(&root.join("ocaml").join("qp_run")),
            Some(root.clone())
        );
        assert_eq!(
            derive_root(&root.join("bin").join("qp_run")),
            Some(root.clone())
        );
    }

    /// A tree without `data/executables` has no root, and the error text is
    /// the reference's.
    #[test]
    fn missing_tree_is_an_error() {
        let dir = tempfile::tempdir().expect("tempdir");
        let exe = dir.path().join("a").join("b").join("qp_run");
        assert_eq!(derive_root(&exe), None);
        let error = resolve(None, &exe).expect_err("no root");
        assert!(error.to_string().contains("QP_ROOT environment variable"));
    }

    /// The environment wins over the executable path, and an empty value is
    /// treated as unset.
    #[test]
    fn environment_has_priority() {
        let (_dir, root, exe) = fake_root();
        let other = PathBuf::from("/some/other/root");
        assert_eq!(resolve(Some("/some/other/root"), &exe).expect("env"), other);
        assert_eq!(resolve(Some(""), &exe).expect("derived"), root);
        assert_eq!(resolve(Some("   "), &exe).expect("derived"), root);
    }
}
