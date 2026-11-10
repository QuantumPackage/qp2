//! The `qp_tunnel` command line (wave S2: clap).
//!
//! Every case checks the exit code of the reference's failure class (2 for
//! every setup failure, 0 for `-h`) and at most one stable substring of the
//! new Rust wording; no test pins a whole help block or an OCaml `Fatal
//! error: exception` line.
//!
//! The dead-remote case cannot be a single line: the four proxy threads
//! each fail their `connect` and print their own line before the fatal line
//! names the lowest-id failure. The four thread lines are compared as a set
//! because their order is a thread-scheduling race in both implementations.

mod support;

use std::fs;
use std::path::Path;
use std::process::ExitStatus;
use std::sync::Mutex;

use support::{Running, parse_local_port};

/// Serialises the cases that bind the fixed local window.
///
/// `qp_tunnel` picks its window by binding and closing each candidate
/// (`discover_local_port`), so two test children running at once can both
/// choose `42379`; the loser then fails its *local* bind, and the fatal
/// line names `tcp://*:42379` where these cases assert the *remote*
/// endpoint. One case at a time is deterministic, and every assertion is
/// unchanged. Poisoning is ignored: a failed case has already reported.
static LOCAL_WINDOW: Mutex<()> = Mutex::new(());

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

/// Exit 2, nothing on stdout, and one stable substring on stderr.
fn assert_fatal(status: ExitStatus, stdout: &str, stderr: &str, needle: &str) {
    assert_eq!(status.code(), Some(2), "stderr: {stderr}");
    assert_eq!(stdout, "", "a pre-banner failure prints nothing on stdout");
    assert!(
        stderr.contains(needle),
        "stderr:\n{stderr}\nwas expected to contain {needle:?}"
    );
}

#[test]
fn no_argument_exits_two_with_usage() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &[], Some("/nonexistent"));
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "the following required arguments were not provided",
    );
}

#[test]
fn the_help_flag_prints_help_and_exits_zero() {
    for flag in ["-h", "--help"] {
        let directory = tempdir();
        let (status, stdout, stderr) = run(directory.path(), &[flag], Some("/nonexistent"));
        assert_eq!(status.code(), Some(0), "flag {flag}; stderr: {stderr}");
        assert!(stdout.contains("Usage:"), "flag {flag}; stdout:\n{stdout}");
        assert_eq!(stderr, "", "flag {flag}");
    }
}

#[test]
fn an_unknown_option_exits_two() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &["-z", "ezfio"], Some("/nonexistent"));
    assert_fatal(status, &stdout, &stderr, "unexpected argument '-z' found");
}

#[test]
fn a_value_on_the_flag_is_rejected() {
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
        "unexpected value '1' for '--get-input'",
    );
}

#[test]
fn a_second_positional_is_rejected() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &["a", "b"], Some("/nonexistent"));
    assert_fatal(status, &stdout, &stderr, "unexpected argument 'b' found");
}

#[test]
fn a_malformed_address_names_it() {
    for address in ["tcp:1.2.3.4", "notanaddress"] {
        let directory = tempdir();
        let (status, stdout, stderr) = run(directory.path(), &[address], Some("/nonexistent"));
        assert_fatal(
            status,
            &stdout,
            &stderr,
            &format!("malformed address {address} (expected proto:host:port)"),
        );
    }
}

#[test]
fn a_non_numeric_port_names_it() {
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
        "invalid integer 'notaport' in the port of tcp:1.2.3.4:notaport",
    );
}

#[test]
fn a_hex_port_is_rejected() {
    // Wave S2: the port is a decimal `i64`, so the reference's
    // `int_of_string` hex grammar is gone.
    let directory = tempdir();
    let (status, stdout, stderr) = run(
        directory.path(),
        &["tcp:127.0.0.1:0x10"],
        Some("/nonexistent"),
    );
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "invalid integer '0x10' in the port of tcp:127.0.0.1:0x10",
    );
}

#[test]
fn a_port_above_the_i64_bound_is_rejected() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(
        directory.path(),
        &["tcp:127.0.0.1:9223372036854775808"],
        Some("/nonexistent"),
    );
    assert_fatal(status, &stdout, &stderr, "invalid integer");
}

#[test]
fn double_dash_passes_the_positional() {
    let directory = tempdir();
    let (status, stdout, stderr) = run(
        directory.path(),
        &["--", "notanaddress"],
        Some("/nonexistent"),
    );
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "malformed address notanaddress (expected proto:host:port)",
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
    assert_fatal(status, &stdout, &stderr, &format!("{}", address.display()));
    // `ensure_work` created `work/` before the read failed.
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
    assert_fatal(status, &stdout, &stderr, "end of file");
}

#[test]
fn qp_root_is_checked_before_help() {
    // Even `-h` dies: the reference reads Qpackage.root from a module
    // initialiser, before the parser.
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &["-h"], None);
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "QP_ROOT is not set (source quantum_package.rc)",
    );
}

#[test]
fn an_empty_qp_root_is_accepted() {
    // Presence is what is checked, not validity: the run reaches the
    // address and fails there.
    let directory = tempdir();
    let (status, stdout, stderr) = run(directory.path(), &["notanaddress"], Some(""));
    assert_fatal(
        status,
        &stdout,
        &stderr,
        "malformed address notanaddress (expected proto:host:port)",
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
        format!("Connect to the tunnel with: tcp://127.0.0.1:{local_port}\n")
    );
    assert!(!stdout.contains("ready"));

    // Each thread's `connect` fails on the single-slash endpoint, and no
    // stop request can suppress the setup of the other threads.
    let expected = [
        "proxy thread 1 stopped: cannot connect to tcp:127.0.0.1:1",
        "proxy thread 2 stopped: cannot connect to tcp:127.0.0.1:3",
        "proxy thread 3 stopped: cannot connect to tcp:127.0.0.1:2",
        "proxy thread 4 stopped: cannot connect to tcp:127.0.0.1:10",
        "qp_tunnel: cannot connect to tcp:127.0.0.1:1",
    ];
    assert_dead_remote(&stdout, &stderr, &expected);
}

#[test]
fn a_negative_port_keeps_its_sign() {
    // `port + 9` renders `4` for a port of -5.
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
        "proxy thread 1 stopped: cannot connect to tcp:127.0.0.1:-5",
        "proxy thread 2 stopped: cannot connect to tcp:127.0.0.1:-3",
        "proxy thread 3 stopped: cannot connect to tcp:127.0.0.1:-4",
        "proxy thread 4 stopped: cannot connect to tcp:127.0.0.1:4",
        "qp_tunnel: cannot connect to tcp:127.0.0.1:-5",
    ];
    assert_dead_remote(&stdout, &stderr, &expected);
}

/// The dead-remote shape: the banner on stdout, then the five stderr lines
/// (four thread lines as a set, the fatal line last).
fn assert_dead_remote(stdout: &str, stderr: &str, expected: &[&str]) {
    let local_port = parse_local_port(stdout);
    assert_eq!(
        stdout,
        format!("Connect to the tunnel with: tcp://127.0.0.1:{local_port}\n")
    );
    assert!(!stdout.contains("ready"));

    let mut lines: Vec<&str> = stderr.lines().collect();
    assert_eq!(lines.len(), expected.len(), "stderr:\n{stderr}");
    // The fatal line is last; the four thread lines are a set.
    assert_eq!(lines.pop(), Some(expected[4]));
    lines.sort_unstable();
    let mut wanted: Vec<&str> = expected[..4].to_vec();
    wanted.sort_unstable();
    assert_eq!(lines, wanted);
}
