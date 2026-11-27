//! The launcher body of `qp_run`: program resolution, the banner and the
//! run-configuration block, the `qp_edit -c` pre-check, the job's exit
//! status, and the two early failures (missing EZFIO directory, unknown
//! program). The job inherits the environment, minus the task-server
//! address variables and plus the one buffering switch: no address file,
//! no launcher-created `work/`.

mod common;

use common::{Fixture, observe_script, observed};

/// The job sees the EZFIO directory as its first argument and none of the
/// task-server environment; the launcher creates no `work/` and no address
/// file; the banner is printed.
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

    assert!(outcome.stdout.contains("Quantum Package"), "{outcome}");
    assert!(
        outcome
            .stdout
            .contains(&format!("EZFIO Dir : {}", fixture.ezfio_str())),
        "the banner names the EZFIO directory\n{outcome}"
    );
}

/// The job's Fortran runtime does not buffer a redirected log.
///
/// libgfortran keeps stdout in memory until `exit` when it is a regular file,
/// so `qp_run fci h2o.ezfio > fci.out` left the file empty for the whole run.
/// The launcher sets the runtime's own switch for every job; the value is
/// asserted here because the shim is a shell script and cannot buffer.
#[test]
fn the_job_output_is_not_buffered() {
    let fixture = Fixture::new();
    fixture.program("scf", &observe_script("exit 0"));
    let outcome = fixture.run(&["scf", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 0, "{outcome}");
    assert_eq!(
        observed(&fixture.work())
            .get("unbuffered")
            .map(String::as_str),
        Some("1"),
        "GFORTRAN_UNBUFFERED_PRECONNECTED must reach the job\n{outcome}"
    );
}

/// A non-zero job status is reported and becomes the tool's status.
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

/// `-p/--prefix` runs the job through the prefix, end to end: the prefix
/// command receives exactly the program path and the EZFIO directory, and
/// the job still runs with the EZFIO directory as its first argument.
/// (The `gdb` special case of `qp_cli::launch::job_command` — `gdb --args `
/// — has no real `gdb` on every test machine; the mapping itself is unit
/// tested there.)
#[test]
fn the_prefix_is_run_before_the_job() {
    let fixture = Fixture::new();
    let seen = fixture.dir().join("prefix_seen");
    fixture.shim(
        "prefixer",
        &format!("printf '%s\\n' \"$@\" > {}\nexec \"$@\"", seen.display()),
    );
    let prefix = fixture.root().join("bin").join("prefixer");
    let program = fixture.program("scf", &observe_script("exit 0"));
    let outcome = fixture.run(&[
        "-p",
        prefix.to_str().expect("UTF-8 path"),
        "scf",
        fixture.ezfio_str(),
    ]);
    assert_eq!(outcome.code(), 0, "{outcome}");

    let prefix_args = std::fs::read_to_string(&seen).expect("the prefix ran");
    let mut lines = prefix_args.lines();
    assert_eq!(
        lines.next(),
        Some(program.to_str().expect("UTF-8 path")),
        "the prefix receives the program path"
    );
    assert_eq!(
        lines.next(),
        Some(fixture.ezfio_str()),
        "the prefix receives the EZFIO directory"
    );
    assert_eq!(
        lines.next(),
        None,
        "the prefix receives exactly two arguments"
    );

    let observed = observed(&fixture.work());
    assert_eq!(
        observed.get("arg1").map(String::as_str),
        Some(fixture.ezfio_str()),
        "the job ran through the prefix with the EZFIO directory as $1"
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
    assert!(!outcome.stdout.contains("EZFIO Dir : "), "{outcome}");
}

/// The run-configuration block prints one line per group of readable keys,
/// with the fixed labels and field order; a group with no readable key
/// (here `AO integrals` / `MO integrals`) prints no line at all. The
/// fixture's EZFIO is created bare, so the test writes the keys of the
/// groups it asserts.
#[test]
fn the_run_configuration_block_prints_the_readable_groups() {
    let fixture = Fixture::new();
    {
        let ezfio = qp_ezfio::EzFio::open(fixture.ezfio()).expect("open the fixture EZFIO");
        ezfio
            .set("determinants", "n_states", &2i64)
            .expect("n_states");
        ezfio.set("determinants", "s2_eig", &false).expect("s2_eig");
        ezfio
            .set("determinants", "n_det_max", &100_000i64)
            .expect("n_det_max");
        ezfio
            .set("determinants", "weight_selection", &0i64)
            .expect("weight_selection");
        ezfio.set("perturbation", "do_pt2", &false).expect("do_pt2");
        ezfio
            .set("perturbation", "pt2_max", &1e-6f64)
            .expect("pt2_max");
        ezfio
            .set("scf_utils", "scf_algorithm", &"DF".to_string())
            .expect("scf_algorithm");
        ezfio
            .set("scf_utils", "n_it_scf_max", &100i64)
            .expect("n_it_scf_max");
        ezfio
            .set("scf_utils", "thresh_scf", &1e-13f64)
            .expect("thresh_scf");
        ezfio
            .set("davidson_keywords", "davidson_sze_max", &50i64)
            .expect("davidson_sze_max");
    }
    fixture.program("scf", "exit 0");
    let outcome = fixture.run(&["scf", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 0, "{outcome}");
    assert!(
        outcome.stdout.contains(
            "  Run: n_states=2  s2_eig=F  n_det_max=100000  \
             weight_selection=0  do_pt2=F  pt2_max=1e-6",
        ),
        "the Run: line names the readable determinants and perturbation keys\n{outcome}"
    );
    assert!(
        outcome
            .stdout
            .contains("  SCF: algorithm=DF  n_it_scf_max=100  thresh_scf=1e-13"),
        "the SCF: line names the readable scf_utils keys\n{outcome}"
    );
    assert!(
        outcome.stdout.contains("  Davidson: sze_max=50"),
        "the Davidson: line names davidson_sze_max\n{outcome}"
    );
    assert!(
        !outcome.stdout.contains("AO integrals") && !outcome.stdout.contains("MO integrals"),
        "a group with no readable key prints no line at all\n{outcome}"
    );
}
