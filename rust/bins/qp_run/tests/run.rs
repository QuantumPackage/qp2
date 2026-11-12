//! The launcher body of `qp_run`: program resolution, the banner, the
//! `qp_edit -c` pre-check, the job's exit status, and the two early failures
//! (missing EZFIO directory, unknown program). The job inherits the
//! environment: no task-server address variables, no address file, no
//! launcher-created `work/`.

mod common;

use common::{Fixture, observe_script, observed};

/// The job sees the EZFIO directory as its first argument and none of the
/// task-server environment; the launcher creates no `work/` and no address
/// file; the banner and the wall time are printed.
#[test]
fn the_job_inherits_no_task_server_environment() {
    let fixture = Fixture::new();
    fixture.program("scf", &observe_script("exit 0"));
    let outcome = fixture.run(&["scf", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 0, "{outcome}");

    let seen = observed(&fixture.work());
    assert_eq!(
        seen.get("arg1").map(String::as_str),
        Some(fixture.ezfio_str()),
        "the job's first argument is the EZFIO directory"
    );
    assert_eq!(
        seen.get("address").map(String::as_str),
        Some("UNSET"),
        "QP_RUN_ADDRESS is not set for the job\n{outcome}"
    );
    assert_eq!(
        seen.get("master").map(String::as_str),
        Some("UNSET"),
        "QP_RUN_ADDRESS_MASTER is not set for the job\n{outcome}"
    );
    assert_eq!(
        seen.get("file").map(String::as_str),
        Some("absent"),
        "work/qp_run_address is never created\n{outcome}"
    );
    assert_eq!(
        seen.get("workdir").map(String::as_str),
        Some("absent"),
        "qp_run creates no work/; the Fortran executable does\n{outcome}"
    );
    assert!(
        !fixture.work().join("qp_run_address").exists(),
        "the address file is absent after the run"
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

/// A non-zero job status is reported and becomes the tool's status, after
/// the wall time was printed.
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

/// A failing `qp_edit -c` is `the input is inconsistent (qp_edit -c failed)`,
/// exit 2, and nothing else runs: no job, no `work/`, no address file.
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
    let outcome = common::Running::spawn(command, &fixture.logs()).wait(common::TIMEOUT);

    assert_eq!(outcome.code(), 2, "{outcome}");
    assert_eq!(
        outcome.stderr, "qp_run: the input is inconsistent (qp_edit -c failed)\n",
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
        format!(
            "qp_run: the EZFIO directory does not exist: {}\n",
            missing.display()
        )
    );
}

/// An unknown program prints the visible choices on standard output and
/// fails with the message `program not found in data/executables`, exit 2,
/// before the banner.
#[test]
fn an_unknown_program_lists_the_choices() {
    let fixture = Fixture::new();
    fixture.program("scf", "exit 0");
    let outcome = fixture.run(&["nope", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert!(outcome.stdout.contains("Possible choices:"), "{outcome}");
    assert!(outcome.stdout.contains("* scf"), "{outcome}");
    assert_eq!(
        outcome.stderr,
        "qp_run: program not found in data/executables: nope\n"
    );
    assert!(!outcome.stdout.contains("Quantum Package"), "{outcome}");
}
