//! The master-mode body of `qp_run` (`ocaml/qp_run.ml:19-227`): the address
//! the job inherits, the address file's lifetime, the job's exit status, the
//! `qp_edit -c` pre-check and the two early failures (missing EZFIO directory,
//! unknown program).

mod common;

use common::{Fixture, Running, TIMEOUT, observe_script, observed};

/// The port the task server advertises, from its startup line.
fn advertised_port(stdout: &str) -> u32 {
    let line = stdout
        .lines()
        .find(|line| line.starts_with("Task server running : "))
        .unwrap_or_else(|| panic!("no task-server line in:\n{stdout}"));
    line.rsplit(':')
        .next()
        .expect("the address ends with a port")
        .parse()
        .expect("numeric port")
}

/// The job sees `QP_RUN_ADDRESS`, the address file exists and matches it
/// during the run, and both are gone afterwards; the banner and the wall time
/// are printed.
#[test]
fn the_job_inherits_the_advertised_address() {
    let fixture = Fixture::new();
    fixture.program("scf", &observe_script("exit 0"));
    let outcome = fixture.run(&["scf", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 0, "{outcome}");

    let address = format!("tcp://127.0.0.1:{}", advertised_port(&outcome.stdout));
    let seen = observed(&fixture.work());
    assert_eq!(
        seen.get("address").map(String::as_str),
        Some(address.as_str()),
        "{outcome}"
    );
    assert_eq!(
        seen.get("arg1").map(String::as_str),
        Some(fixture.ezfio_str()),
        "the job's first argument is the EZFIO directory"
    );
    assert_eq!(
        seen.get("master").map(String::as_str),
        Some(""),
        "master mode sets no QP_RUN_ADDRESS_MASTER"
    );
    assert_eq!(
        seen.get("file").map(String::as_str),
        Some(address.as_str()),
        "work/qp_run_address holds the advertised address while the job runs"
    );
    assert!(
        !fixture.work().join("qp_run_address").exists(),
        "the address file is removed when the run ends"
    );
    assert!(
        fixture.work().join("empty").exists(),
        "ezfio_work created the work/ marker"
    );

    assert!(outcome.stdout.contains("Date: "), "{outcome}");
    assert!(outcome.stdout.contains("Quantum Package"), "{outcome}");
    assert!(
        outcome
            .stdout
            .contains(&format!("EZFIO Dir : {}", fixture.ezfio_str())),
        "{outcome}"
    );
    assert!(outcome.stdout.contains("Wall time: "), "{outcome}");
}

/// A non-zero job status is reported and becomes the tool's status, after the
/// server stopped and the wall time was printed.
#[test]
fn the_job_exit_code_is_propagated() {
    let fixture = Fixture::new();
    fixture.program("failing", "exit 3");
    let outcome = fixture.run(&["failing", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 3, "{outcome}");
    assert!(
        outcome.stdout.contains("Program exited with code 3."),
        "{outcome}"
    );
    assert!(
        outcome.stdout.contains("Wall time: "),
        "the reference prints the wall time before exiting:\n{outcome}"
    );
    assert!(!fixture.work().join("qp_run_address").exists());
}

/// A failing `qp_edit -c` is `Error: Input inconsistent`, exit 2, and nothing
/// else runs: no server, no job, no address file.
#[test]
fn a_failing_input_check_stops_the_run() {
    let fixture = Fixture::new();
    fixture.program("scf", &observe_script("exit 0"));
    let record = fixture.dir().join("edit_arguments");
    fixture.shim(
        "qp_edit_fail",
        &format!("printf '%s\\n' \"$*\" > {}\nexit 1", record.display()),
    );
    let mut command = fixture.command(&["scf", fixture.ezfio_str()]);
    command.env(
        "QP_RUN_QP_EDIT",
        fixture.root().join("bin").join("qp_edit_fail"),
    );
    let outcome = Running::spawn(command, &fixture.logs()).wait(TIMEOUT);

    assert_eq!(outcome.code(), 2, "{outcome}");
    assert!(
        outcome
            .stderr
            .starts_with("qp_run: Error: Input inconsistent"),
        "{outcome}"
    );
    // The pre-check is `"<editor> -c <ezfio>"`, like `Sys.command`.
    assert_eq!(
        std::fs::read_to_string(&record)
            .expect("the editor ran")
            .trim(),
        format!("-c {}", fixture.ezfio_str())
    );
    assert!(
        !fixture.work().join("observed").exists(),
        "the job must not run\n{outcome}"
    );
    assert!(!fixture.work().join("qp_run_address").exists());
    assert!(!outcome.stdout.contains("Task server running"), "{outcome}");
}

/// A directory that does not exist fails before any other check.
#[test]
fn a_missing_ezfio_directory_is_reported() {
    let fixture = Fixture::new();
    let missing = fixture.dir().join("not-an-ezfio");
    let outcome = fixture.run(&["scf", missing.to_str().expect("UTF-8 path")]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert_eq!(
        outcome.stderr,
        format!("qp_run: EZFIO directory {} not found\n", missing.display())
    );
}

/// An unknown program prints the visible choices on standard output and fails
/// with the reference message, exit 2, before the banner.
#[test]
fn an_unknown_program_lists_the_choices() {
    let fixture = Fixture::new();
    fixture.program("scf", "exit 0");
    let outcome = fixture.run(&["nope", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert!(outcome.stdout.contains("Possible choices:"), "{outcome}");
    assert!(outcome.stdout.contains("* scf"), "{outcome}");
    assert_eq!(outcome.stderr, "qp_run: Executable nope not found\n");
    assert!(!outcome.stdout.contains("Quantum Package"), "{outcome}");
}
