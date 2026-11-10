//! Slave mode (`qp_run -s`, `ocaml/qp_run.ml:129-141`).
//!
//! A slave still starts its own task server; it differs by skipping the
//! `qp_edit -c` pre-check, reading the master's address from
//! `<work>/qp_run_address` instead of writing it, and leaving that file in
//! place when it exits.

mod common;

use common::{Fixture, Running, TIMEOUT, observe_script, observed};

/// The master's address, as a master `qp_run` would have written it.
const MASTER: &str = "tcp://10.1.2.3:45678";

/// Write the address file a master leaves behind.
fn write_master(fixture: &Fixture) {
    std::fs::create_dir_all(fixture.work()).expect("create work/");
    std::fs::write(fixture.work().join("qp_run_address"), format!("{MASTER}\n"))
        .expect("write the master address");
}

/// The job sees both variables: `QP_RUN_ADDRESS` is the slave's own server,
/// `QP_RUN_ADDRESS_MASTER` is the file's address. The file is not removed and
/// the `qp_edit` pre-check is skipped.
#[test]
fn a_slave_reads_the_master_and_keeps_the_file() {
    let fixture = Fixture::new();
    fixture.program("report", &observe_script("exit 0"));
    write_master(&fixture);
    // If the slave ran the editor, this failing shim would abort the run.
    let mut command = fixture.command(&["-s", "report", fixture.ezfio_str()]);
    command.env(
        "QP_RUN_QP_EDIT",
        fixture.root().join("bin").join("qp_edit_fail"),
    );
    let outcome = Running::spawn(command, &fixture.logs()).wait(TIMEOUT);
    assert_eq!(outcome.code(), 0, "{outcome}");

    let seen = observed(&fixture.work());
    assert_eq!(
        seen.get("master").map(String::as_str),
        Some(MASTER),
        "{outcome}"
    );
    let address = seen.get("address").expect("QP_RUN_ADDRESS is set");
    assert!(
        address.starts_with("tcp://127.0.0.1:"),
        "the slave advertises its own server, not the master's: {address:?}"
    );
    assert_ne!(address, MASTER);
    assert_eq!(
        seen.get("file").map(String::as_str),
        Some(MASTER),
        "the job still sees the master's file while it runs"
    );
    assert!(outcome.stdout.contains("Wall time: "), "{outcome}");

    assert_eq!(
        std::fs::read_to_string(fixture.work().join("qp_run_address")).expect("the file is kept"),
        format!("{MASTER}\n"),
        "a slave does not remove the master's file"
    );
}

/// No address file: the missing-master failure, exit 2, before any
/// server or job.
#[test]
fn a_slave_without_an_address_file_fails() {
    let fixture = Fixture::new();
    fixture.program("report", &observe_script("exit 0"));
    std::fs::create_dir_all(fixture.work()).expect("create work/");
    let outcome = fixture.run(&["-s", "report", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert_eq!(
        outcome.stderr,
        "qp_run: no master is running (work/qp_run_address missing or empty)\n"
    );
    assert!(
        !fixture.work().join("observed").exists(),
        "the job must not run\n{outcome}"
    );
    assert!(!outcome.stdout.contains("Task server running"), "{outcome}");
}

/// An empty address file is treated like a missing one: there is no usable
/// master. (The reference raises `End_of_file` there — also exit 2, with a
/// different message; the deviation is documented on `read_master_address`.)
#[test]
fn an_empty_address_file_is_no_master() {
    let fixture = Fixture::new();
    fixture.program("report", &observe_script("exit 0"));
    std::fs::create_dir_all(fixture.work()).expect("create work/");
    std::fs::write(fixture.work().join("qp_run_address"), "").expect("write an empty file");
    let outcome = fixture.run(&["-s", "report", fixture.ezfio_str()]);
    assert_eq!(outcome.code(), 2, "{outcome}");
    assert_eq!(
        outcome.stderr,
        "qp_run: no master is running (work/qp_run_address missing or empty)\n"
    );
}
