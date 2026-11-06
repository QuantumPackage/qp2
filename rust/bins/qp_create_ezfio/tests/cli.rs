//! End-to-end tests for the `qp_create_ezfio` front end.
//!
//! These run the real binary, so they cover what the unit tests cannot:
//! `QP_ROOT` handling, stdout/stderr routing, the exit status, the
//! `argv[0]` substitution in the help block, and the non-UTF-8 argument
//! guard. Every expected string was measured against the reference build
//! (`~/qp2/ocaml/qp_create_ezfio.native`), except the non-UTF-8 message,
//! which has no reference counterpart (documented in `main.rs`).

use std::process::{Command, Output};

/// The reference help capture, produced with `argv[0] = "qp_create_ezfio"`.
const HELP: &str = include_str!("data/help_bare.txt");

/// A hydrogen basis fixture; absolute so the action finds it regardless of
/// the working directory.
const H_BASIS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/h.basis");
/// Two hydrogens 0.7 angstrom apart.
const H2_XYZ: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/h2.xyz");

/// A command whose `QP_ROOT` is set (to a value that need not exist) and
/// whose trailing `qp_edit` is the `true` seam (the tests must not depend
/// on a host `qp_edit`; the reference's `Sys.command` runs it through
/// `/bin/sh` and would print the shell's `command not found` on stderr).
fn command(args: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_qp_create_ezfio"));
    command
        .env("QP_ROOT", "/nonexistent")
        .env("QP_EDIT", "true")
        .args(args);
    command
}

/// Run the binary and capture both streams.
fn run(args: &[&str]) -> Output {
    command(args)
        .output()
        .expect("the compiled test binary must be runnable")
}

/// Decode stderr, which is ASCII in every case exercised here.
fn stderr(output: &Output) -> String {
    String::from_utf8(output.stderr.clone()).expect("stderr is valid UTF-8")
}

#[cfg(unix)]
#[test]
fn help_matches_the_reference_byte_for_byte() {
    use std::os::unix::process::CommandExt;

    let output = command(&["-h"])
        .arg0("qp_create_ezfio")
        .output()
        .expect("the compiled test binary must be runnable");
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(output.stdout, HELP.as_bytes());
    // `-h` prints the help and parsing continues; with no FILE it then fails.
    assert_eq!(stderr(&output), "qp_create_ezfio: input file is missing\n");
}

#[test]
fn missing_qp_root_fails_before_parsing() {
    let output = Command::new(env!("CARGO_BIN_EXE_qp_create_ezfio"))
        .env_remove("QP_ROOT")
        .args(["-h", "-b", "x", "f.xyz"])
        .output()
        .expect("the compiled test binary must be runnable");
    assert_eq!(output.status.code(), Some(2));
    assert!(
        output.stdout.is_empty(),
        "the QP_ROOT failure must precede the help block"
    );
    assert_eq!(
        stderr(&output),
        "Fatal error: exception Failure(\"QP_ROOT environment variable is not set.\\nPlease source the quantum_package.rc file.\")\n"
    );
}

#[test]
fn empty_qp_root_is_accepted() {
    // An empty root is a valid presence check result. Both inputs are
    // absolute local files, so the empty root never has to resolve a
    // database name.
    let work = tempfile::tempdir().expect("tempdir");
    let out = work.path().join("h2.ezfio");
    let output = Command::new(env!("CARGO_BIN_EXE_qp_create_ezfio"))
        .env("QP_ROOT", "")
        .env("QP_EDIT", "true")
        .args(["-b", H_BASIS, "-o", &out.to_string_lossy(), H2_XYZ])
        .output()
        .expect("the compiled test binary must be runnable");
    assert_eq!(output.status.code(), Some(0), "stderr: {}", stderr(&output));
    assert_eq!(
        String::from_utf8(output.stdout).expect("stdout is valid UTF-8"),
        format!("{}\n", out.display())
    );
    assert!(output.stderr.is_empty());
}

#[test]
fn front_end_failures_match_the_reference() {
    let cases: &[(&[&str], &str)] = &[
        (&[], "qp_create_ezfio: --basis option is missing.\n"),
        (&["-b", "x"], "qp_create_ezfio: input file is missing\n"),
        (
            &["-Z"],
            "Fatal error: exception Getopt.Error(\"Unknown option -Z\")\n",
        ),
        (
            &["--basis"],
            "Fatal error: exception Getopt.Error(\"Option --basis must have an argument\")\n",
        ),
        (
            &["--au=x"],
            "Fatal error: exception Getopt.Error(\"Option --au does not accept argument\")\n",
        ),
        (
            &["-b", ""],
            "Fatal error: exception Invalid_argument(\"index out of bounds\")\n",
        ),
        (
            &["-b", "x", "-c", "abc", "f.xyz"],
            "qp_create_ezfio: int_of_string\n",
        ),
        (
            &["-b", "x", "--ghost=", "f.xyz"],
            "qp_create_ezfio: float_of_string\n",
        ),
    ];
    for (args, expected) in cases {
        let output = run(args);
        assert_eq!(output.status.code(), Some(2), "args {args:?}");
        assert_eq!(stderr(&output), *expected, "args {args:?}");
        assert!(output.stdout.is_empty(), "args {args:?}");
    }
}

#[test]
fn valid_invocation_writes_the_default_named_ezfio() {
    // Without `-o`, the directory name comes from the input: the last
    // `.xyz` extension is replaced by `.ezfio`.
    let work = tempfile::tempdir().expect("tempdir");
    std::fs::copy(H_BASIS, work.path().join("h.basis")).expect("copy basis");
    std::fs::copy(H2_XYZ, work.path().join("h2.xyz")).expect("copy xyz");
    let output = command(&["-b", "h.basis", "h2.xyz"])
        .current_dir(work.path())
        .output()
        .expect("the compiled test binary must be runnable");
    assert_eq!(output.status.code(), Some(0), "stderr: {}", stderr(&output));
    // The name goes to stdout as given (relative, like the reference).
    assert_eq!(
        String::from_utf8(output.stdout).expect("stdout is valid UTF-8"),
        "h2.ezfio\n"
    );
    assert!(work.path().join("h2.ezfio").is_dir());
    assert!(output.stderr.is_empty());
}

#[test]
fn terminator_makes_a_later_dash_h_a_file_name() {
    // A real file named `-h` in the working directory: after `--` it must be
    // opened as the input, not consumed as the help flag.
    let work = tempfile::tempdir().expect("tempdir");
    std::fs::write(work.path().join("-h"), "2\nh2\nH 0 0 0\nH 0 0 0.7\n").expect("write fixture");
    std::fs::copy(H_BASIS, work.path().join("h.basis")).expect("copy fixture");
    let output = command(&["-b", "h.basis", "--", "-h"])
        .current_dir(work.path())
        .output()
        .expect("the compiled test binary must be runnable");
    assert_eq!(output.status.code(), Some(0), "stderr: {}", stderr(&output));
    // The file `-h` is the input: the default name has no dot, so the
    // whole name plus `.ezfio`. And no help was printed.
    assert_eq!(
        String::from_utf8(output.stdout).expect("stdout is valid UTF-8"),
        "-h.ezfio\n"
    );
    assert!(work.path().join("-h.ezfio").is_dir());
    assert!(output.stderr.is_empty());
}

#[cfg(unix)]
#[test]
fn help_is_printed_before_the_action() {
    use std::os::unix::process::CommandExt;

    // `-h` prints the help and parsing continues: the action still runs,
    // and the directory line follows the help block on stdout.
    let work = tempfile::tempdir().expect("tempdir");
    let out = work.path().join("h2.ezfio");
    let output = command(&["-h", "-b", H_BASIS, "-o", &out.to_string_lossy(), H2_XYZ])
        .arg0("qp_create_ezfio")
        .output()
        .expect("the compiled test binary must be runnable");
    assert_eq!(output.status.code(), Some(0), "stderr: {}", stderr(&output));
    assert_eq!(
        output.stdout,
        format!("{HELP}{}\n", out.display()).into_bytes()
    );
    assert!(output.stderr.is_empty());
}

#[cfg(unix)]
#[test]
fn non_utf8_argument_is_rejected_without_panicking() {
    use std::ffi::OsString;
    use std::os::unix::ffi::OsStringExt;

    let output = command(&["-b", "x"])
        .arg(OsString::from_vec(vec![0xff, b'.', b'x', b'y', b'z']))
        .output()
        .expect("the compiled test binary must be runnable");
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert_eq!(
        stderr(&output),
        "qp_create_ezfio: command line argument 3 is not valid UTF-8\n"
    );
}
