//! Byte-level checks of the `qp_tunnel` command line against the captured
//! reference outputs.
//!
//! Every case supplies `argv[0] = "qp_tunnel"` (through
//! [`std::os::unix::process::CommandExt::arg0`]) so the help block is the
//! 481-byte capture; the one long-`argv[0]` case re-invokes the binary the
//! way a renamed reference copy was, and the help block re-wraps with it
//! (measured `Format` behaviour, reproduced by the bin-local parser).
//!
//! The dead-remote case cannot be byte-identical: the reference prints its
//! four `Thread … killed` lines and then hangs in `zmq_ctx_term` with
//! `Unix.EINTR`; this port reaps the threads and exits 2 with the lowest-id
//! failure. The four thread lines are compared as a set because their order
//! is a thread-scheduling race in both implementations.

mod support;

use std::fs;
use std::path::Path;
use std::process::ExitStatus;
use std::sync::Mutex;

use support::{HELP, Running, parse_local_port};

/// Serialises the cases that bind the fixed local window.
///
/// `qp_tunnel` picks its window by binding and closing each candidate
/// (`discover_local_port`), so two test children running at once can both
/// choose `42379`; the loser then fails its *local* bind, and the fatal
/// line names `tcp://*:42379` where these cases assert the *remote*
/// endpoint. One case at a time is deterministic, and every assertion is
/// unchanged. Poisoning is ignored: a failed case has already reported.
static LOCAL_WINDOW: Mutex<()> = Mutex::new(());

/// `Failure("EZFIO_FILE or ADDRESS is missing")`, as the runtime prints it.
const MISSING_ARGUMENT: &str =
    "Fatal error: exception Failure(\"EZFIO_FILE or ADDRESS is missing\")\n";

/// The uncaught `Failure` of `Qpackage.root`, escaped by `Printexc`.
const QP_ROOT_MISSING: &str = "Fatal error: exception Failure(\"QP_ROOT environment variable is not set.\\nPlease source the quantum_package.rc file.\")\n";

/// Run one command line to completion in `directory`.
fn run(
    directory: &Path,
    arguments: &[&str],
    qp_root: Option<&str>,
) -> (ExitStatus, String, String) {
    let mut child = Running::spawn(directory, arguments, qp_root);
    let status = child.wait();
    (status, child.stdout(), child.stderr())
}

/// A fresh directory per case, removed when the test ends.
fn tempdir() -> tempfile::TempDir {
    tempfile::tempdir().expect("tempdir")
}

/// Every failure exits 2 and leaves stdout empty.
fn assert_fatal(status: ExitStatus, stdout: &str, stderr: &str, expected: &str) {
    assert_eq!(status.code(), Some(2), "stderr: {stderr}");
    assert_eq!(stdout, "", "a pre-banner failure prints nothing on stdout");
    assert_eq!(stderr, expected);
}

#[test]
fn no_argument_prints_one_help_block_and_fails() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &[], Some("/nonexistent"));
    assert_eq!(status.code(), Some(2));
    assert_eq!(stdout, HELP, "the reference prints the help exactly once");
    assert_eq!(stderr, MISSING_ARGUMENT);
}

#[test]
fn help_alone_prints_the_help_twice() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &["-h"], Some("/nonexistent"));
    assert_eq!(status.code(), Some(2));
    // Once from `-h`, once from the anonymous-argument check, as measured.
    assert_eq!(stdout, format!("{HELP}{HELP}"));
    assert_eq!(stdout.len(), 962);
    assert_eq!(stderr, MISSING_ARGUMENT);
}

#[test]
fn help_before_a_malformed_address_prints_once() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(
        directory.path(),
        &["-h", "notanaddress"],
        Some("/nonexistent"),
    );
    assert_eq!(status.code(), Some(2));
    assert_eq!(stdout, HELP);
    assert_eq!(
        stderr,
        "Fatal error: exception Failure(\"notanaddress : Malformed address\")\n"
    );
}

#[test]
fn malformed_addresses_match_the_reference() {
    for address in ["tcp:1.2.3.4", "notanaddress"] {
        let directory = tempdir();
        let (status, stdout, stderr) = run(directory.path(), &[address], Some("/nonexistent"));
        assert_fatal(
            status,
            &stdout,
            &stderr,
            &format!("Fatal error: exception Failure(\"{address} : Malformed address\")\n"),
        );
    }
}

#[test]
fn a_non_numeric_port_is_int_of_string() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(
        directory.path(),
        &["tcp:1.2.3.4:notaport"],
        Some("/nonexistent"),
    );
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "Fatal error: exception Failure(\"int_of_string\")\n",
    );
}

#[test]
fn a_missing_address_file_names_the_path() {
    let directory = tempdir();
    let ezfio = directory.path().join("h2.ezfio");
    fs::create_dir_all(&ezfio).expect("create the EZFIO directory");
    let (status, stdout, stderr) = run(
        directory.path(),
        &[ezfio.to_str().unwrap()],
        Some("/nonexistent"),
    );
    let address = ezfio.join("work").join("qp_run_address");
    assert_fatal(
        status,
        &stdout,
        &stderr,
        &format!(
            "Fatal error: exception Sys_error(\"{}: No such file or directory\")\n",
            address.display()
        ),
    );
    // The reference's `Qpackage.ezfio_work` created `work/` before failing.
    assert!(ezfio.join("work").is_dir());
}

#[test]
fn an_empty_address_file_is_end_of_file() {
    let directory = tempdir();
    let ezfio = directory.path().join("h3.ezfio");
    fs::create_dir_all(ezfio.join("work")).expect("create work/");
    fs::write(ezfio.join("work").join("qp_run_address"), "").expect("empty address");
    let (status, stdout, stderr) = run(
        directory.path(),
        &[ezfio.to_str().unwrap()],
        Some("/nonexistent"),
    );
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "Fatal error: exception End_of_file\n",
    );
}

#[test]
fn qp_root_is_checked_before_parsing() {
    // Even `-h` dies: the reference reads Qpackage.root from a module
    // initialiser, before `Getopt.parse_cmdline`.
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &["-h"], None);
    assert_fatal(status, &stdout, &stderr, QP_ROOT_MISSING);
}

#[test]
fn an_empty_qp_root_is_accepted() {
    // Presence is what is checked, not validity: the run reaches the parser
    // and fails on the address.
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &["notanaddress"], Some(""));
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "Fatal error: exception Failure(\"notanaddress : Malformed address\")\n",
    );
}

#[test]
fn a_dead_remote_reports_all_four_threads_and_exits() {
    let _window = LOCAL_WINDOW
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let directory = tempdir();
    let (status, stdout, stderr) =
        run(directory.path(), &["tcp:127.0.0.1:1"], Some("/nonexistent"));
    assert_eq!(status.code(), Some(2));

    // Thread 1 starts before the banner and fails before the others; the
    // banner is the only stdout.
    let local_port = parse_local_port(&stdout);
    assert_eq!(
        stdout,
        format!("Connect to:\ntcp://127.0.0.1:{local_port}\n")
    );
    assert!(!stdout.contains("Ready"));

    // The reference's four lines, in any order: each thread's `create_socket`
    // fails on the single-slash endpoint, and no stop request can suppress
    // the setup of the other threads.
    let expected = [
        "Thread 1 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:1.\")",
        "Thread 2 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:3.\")",
        "Thread 3 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:2.\")",
        "Thread 4 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:10.\")",
        "Fatal error: exception Failure(\"Unable to establish connection to tcp:127.0.0.1:1.\")",
    ];
    assert_dead_remote(&stdout, &stderr, &expected);
}

/// The dead-remote shape: the banner on stdout, then the five stderr lines
/// (four thread lines as a set, the lowest-id failure last).
fn assert_dead_remote(stdout: &str, stderr: &str, expected: &[&str]) {
    let local_port = parse_local_port(stdout);
    assert_eq!(
        stdout,
        format!("Connect to:\ntcp://127.0.0.1:{local_port}\n")
    );
    assert!(!stdout.contains("Ready"));

    let mut lines: Vec<&str> = stderr.lines().collect();
    assert_eq!(lines.len(), expected.len(), "stderr:\n{stderr}");
    // The fatal line is last; the four thread lines are a set.
    assert_eq!(lines.pop(), Some(expected[4]));
    lines.sort_unstable();
    let mut wanted: Vec<&str> = expected[..4].to_vec();
    wanted.sort_unstable();
    assert_eq!(lines, wanted);
}

#[test]
fn a_hex_port_is_an_ocaml_int() {
    // `int_of_string "0x10"` is 16; the endpoints render it in decimal
    // (measured against the reference).
    let _window = LOCAL_WINDOW
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let directory = tempdir();
    let (status, stdout, stderr) = run(
        directory.path(),
        &["tcp:127.0.0.1:0x10"],
        Some("/nonexistent"),
    );
    assert_eq!(status.code(), Some(2));
    let expected = [
        "Thread 1 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:16.\")",
        "Thread 2 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:18.\")",
        "Thread 3 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:17.\")",
        "Thread 4 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:25.\")",
        "Fatal error: exception Failure(\"Unable to establish connection to tcp:127.0.0.1:16.\")",
    ];
    assert_dead_remote(&stdout, &stderr, &expected);
}

#[test]
fn a_negative_port_keeps_its_sign() {
    // `int_of_string "-5"` is -5; `port + 9` renders `4` (measured).
    let _window = LOCAL_WINDOW
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let directory = tempdir();
    let (status, stdout, stderr) = run(
        directory.path(),
        &["tcp:127.0.0.1:-5"],
        Some("/nonexistent"),
    );
    assert_eq!(status.code(), Some(2));
    let expected = [
        "Thread 1 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:-5.\")",
        "Thread 2 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:-3.\")",
        "Thread 3 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:-4.\")",
        "Thread 4 killed on uncaught exception Failure(\"Unable to establish connection to tcp:127.0.0.1:4.\")",
        "Fatal error: exception Failure(\"Unable to establish connection to tcp:127.0.0.1:-5.\")",
    ];
    assert_dead_remote(&stdout, &stderr, &expected);
}

#[test]
fn a_port_above_the_ocaml_int_bound_is_int_of_string() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(
        directory.path(),
        &["tcp:127.0.0.1:4611686018427387904"],
        Some("/nonexistent"),
    );
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "Fatal error: exception Failure(\"int_of_string\")\n",
    );
}

#[test]
fn unknown_options_are_getopt_errors() {
    // Measured against the reference: the reference's `Getopt` walk, not a
    // clap diagnostic.
    let cases = [
        (&["-z", "ezfio"], "Getopt.Error(\"Unknown option -z\")"),
        (
            &["--nosuch", "ezfio"],
            "Getopt.Error(\"Unknown option --nosuch\")",
        ),
        (
            &["--nosuch=1", "ezfio"],
            "Getopt.Error(\"Unknown option --nosuch\")",
        ),
        (&["--h", "ezfio"], "Getopt.Error(\"Unknown option --h\")"),
        (&["-g1", "ezfio"], "Getopt.Error(\"Unknown option -1\")"),
        (&["-hx", "ezfio"], "Getopt.Error(\"Unknown option -x\")"),
        (&["-h=x", "ezfio"], "Getopt.Error(\"Unknown option -=\")"),
    ];
    for (arguments, message) in cases {
        let directory = tempdir();
        let (status, stdout, stderr) = run(directory.path(), arguments, Some("/nonexistent"));
        assert_fatal(
            status,
            &stdout,
            &stderr,
            &format!("Fatal error: exception {message}\n"),
        );
    }
}

#[test]
fn a_value_on_a_no_argument_option() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(
        directory.path(),
        &["--get-input=1", "ezfio"],
        Some("/nonexistent"),
    );
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "Fatal error: exception Getopt.Error(\"Option --get-input does not accept argument\")\n",
    );
}

#[test]
fn concatenated_no_argument_flags() {
    // Both flags are set, parsing continues, the address fails.
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &["-gh", "ezfio"], Some("/nonexistent"));
    assert_eq!(status.code(), Some(2));
    assert_eq!(
        stdout, HELP,
        "-h (inside the cluster) printed the help once"
    );
    assert_eq!(
        stderr,
        "Fatal error: exception Failure(\"ezfio : Malformed address\")\n"
    );
}

#[test]
fn double_dash_is_consumed_and_switches_to_anonymous() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &["--", "ezfio"], Some("/nonexistent"));
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "Fatal error: exception Failure(\"ezfio : Malformed address\")\n",
    );

    // Two anonymous arguments after `--`: the second help block, as with
    // any doubled argument.
    let directory = tempdir();
    let (status, stdout, stderr) = run(
        directory.path(),
        &["--", "ezfio", "extra"],
        Some("/nonexistent"),
    );
    assert_eq!(status.code(), Some(2));
    assert_eq!(stdout, HELP);
    assert_eq!(stderr, MISSING_ARGUMENT);
}

#[test]
fn the_help_rewraps_with_a_long_argv0() {
    // The capture was taken from a renamed copy of the reference binary
    // invoked as `./qp_tunnel_really_really_really_long_name_for_wrapping_tests`;
    // invoking the port the same way must reproduce it byte for byte.
    let arg0 = "./qp_tunnel_really_really_really_long_name_for_wrapping_tests";
    let directory = tempdir();
    let mut child =
        support::Running::spawn_with_arg0(directory.path(), &[], Some("/nonexistent"), arg0);
    let status = child.wait();
    assert_eq!(status.code(), Some(2));
    assert_eq!(child.stdout(), include_str!("data/help_long.txt"));
    assert_eq!(child.stderr(), MISSING_ARGUMENT);
}
