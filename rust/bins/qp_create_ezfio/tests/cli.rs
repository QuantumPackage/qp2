//! End-to-end tests for the `qp_create_ezfio` front end (wave S2: clap +
//! the shared `ToolError`).
//!
//! These run the real binary, so they cover what the unit tests cannot:
//! `QP_ROOT` handling, stdout/stderr routing, the exit status, the clap
//! grammar, and the manual check order. Every failure case keeps the
//! reference's scenario (same argv, same fixture mutation), asserts the
//! exit code (2), and pins at most one stable substring — named in the
//! test — on the streams. No exact-stderr pins, no help-byte pins.

#[cfg(unix)]
use std::os::unix::process::CommandExt;
use std::process::{Command, Output};

/// The hydrogen basis fixture; absolute so the action finds it regardless
/// of the working directory.
const H_BASIS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/h.basis");
/// Two hydrogens 0.7 angstrom apart.
const H2_XYZ: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/h2.xyz");

/// The tool name of the one-line standard-error error.
const TOOL: &str = "qp_create_ezfio";

/// A command whose `QP_ROOT` is set (to a value that need not exist) and
/// whose trailing `qp_edit` is the `true` seam (the tests must not depend
/// on a host `qp_edit`).
fn command(args: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_qp_create_ezfio"));
    command
        .env("QP_ROOT", "/nonexistent")
        .env("QP_EDIT", "true")
        .args(args);
    // The usage line of the help is checked through `argv[0]`; `arg0`
    // keeps it the bare program name.
    #[cfg(unix)]
    command.arg0(TOOL);
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
fn the_help_flag_prints_help_and_exits_zero() {
    for flag in ["-h", "--help"] {
        let output = run(&[flag]);
        assert_eq!(output.status.code(), Some(0), "flag {flag}");
        let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
        assert!(stdout.contains("Usage:"), "flag {flag}; stdout: {stdout}");
        // Every option of the reference's table is documented.
        for name in [
            "-o",
            "--output",
            "-b",
            "--basis",
            "-a",
            "--au",
            "-c",
            "--charge",
            "-g",
            "--ghost",
            "-m",
            "--multiplicity",
            "-p",
            "--pseudo",
            "-x",
            "--cartesian",
            "<FILE>",
        ] {
            assert!(stdout.contains(name), "flag {flag} missing {name}");
        }
        assert_eq!(stderr(&output), "", "flag {flag}");
    }
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
        format!("{TOOL}: QP_ROOT is not set (source quantum_package.rc)\n")
    );
}

#[test]
fn an_empty_qp_root_is_accepted() {
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
fn front_end_failures_exit_two() {
    // (args, the stable substring named by the test, the stream it is on).
    let cases: &[(&[&str], &str)] = &[
        // clap: the `FILE` positional is required.
        (&[], "the following required arguments were not provided"),
        (
            &["-b", "x"],
            "the following required arguments were not provided",
        ),
        (
            &["-b", "x", "a.xyz", "b.xyz"],
            "unexpected argument 'b.xyz' found",
        ),
        // clap: unknown options.
        (&["-Z"], "unexpected argument '-Z' found"),
        (&["--nosuch"], "unexpected argument '--nosuch' found"),
        // clap: a value-taking option without a value.
        (&["-b"], "a value is required for '--basis"),
        // clap: a flag given a value.
        (&["--au=x"], "unexpected value 'x' for '--au'"),
        // manual: the mandatory `--basis` check.
        (&["-c", "0", "f.xyz"], "--basis option is missing"),
        // manual: the conversions, charge → ghost → multiplicity,
        // decimal-only (wave S2.3).
        (
            &["-b", "x", "-c", "abc", "f.xyz"],
            "invalid integer 'abc' in charge",
        ),
        (
            &["-b", "x", "--charge=", "f.xyz"],
            "invalid integer '' in charge",
        ),
        (
            &["-b", "x", "-c", "0x10", "f.xyz"],
            "invalid integer '0x10' in charge",
        ),
        (
            &["-b", "x", "-g", "x", "f.xyz"],
            "invalid float 'x' in ghost threshold",
        ),
        (
            &["-b", "x", "--ghost=", "f.xyz"],
            "invalid float '' in ghost threshold",
        ),
        (
            &["-b", "x", "-m", "y", "f.xyz"],
            "invalid integer 'y' in multiplicity",
        ),
        (
            &["-b", "x", "-m", "0x2", "f.xyz"],
            "invalid integer '0x2' in multiplicity",
        ),
    ];
    for (args, needle) in cases {
        let output = run(args);
        assert_eq!(output.status.code(), Some(2), "args {args:?}");
        assert!(
            stderr(&output).contains(needle),
            "args {args:?}; stderr: {}",
            stderr(&output)
        );
        assert!(output.stdout.is_empty(), "args {args:?}");
    }
}

#[test]
fn the_mandatory_basis_check_comes_before_the_conversions() {
    // No `--basis`: the mandatory check wins over the bad charge.
    let output = run(&["-c", "abc", "f.xyz"]);
    assert_eq!(output.status.code(), Some(2));
    assert!(
        stderr(&output).contains("--basis option is missing"),
        "stderr: {}",
        stderr(&output)
    );
    // With `--basis`: the charge conversion fails before the multiplicity
    // one.
    let output = run(&["-b", "x", "-c", "abc", "-m", "yy", "f.xyz"]);
    assert_eq!(output.status.code(), Some(2));
    assert!(
        stderr(&output).contains("invalid integer 'abc' in charge"),
        "stderr: {}",
        stderr(&output)
    );
}

#[test]
fn a_valid_invocation_writes_the_default_named_ezfio() {
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
fn the_terminator_makes_a_later_dash_h_a_file_name() {
    // A real file named `-h` in the working directory: after `--` it must
    // be opened as the input, not consumed as the help flag.
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
fn a_non_utf8_argument_is_rejected_without_panicking() {
    use std::ffi::OsString;
    use std::os::unix::ffi::OsStringExt;

    // clap rejects a non-UTF-8 `argv` byte sequence with its own plain
    // diagnostic (wave S2 deviation; the reference walked raw bytes and
    // rejected it with its own message).
    let output = command(&["-b", "x"])
        .arg(OsString::from_vec(vec![0xff, b'.', b'x', b'y', b'z']))
        .output()
        .expect("the compiled test binary must be runnable");
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert!(
        stderr(&output).contains("invalid UTF-8"),
        "stderr: {}",
        stderr(&output)
    );
}
