//! CLI parity tests for the two `Command_line` specifications that the
//! executables parse with `clap`: `qp_edit` and `qp_run`.
//!
//! The expected messages are transcribed from the reference sources:
//! `scripts/ezfio_interface/qp_edit_template` (the generated `qp_edit`) and
//! `ocaml/qp_run.ml:232–270`. `qp_create_ezfio`, `qp_set_mo_class` and
//! `qp_tunnel` are not here: their binaries parse with a
//! statement-by-statement `getopt` port (P4) and pin that parser against
//! renamed reference copies instead.
//!
//! The parsing rules were checked against the vendored `getopt`
//! (`~/qp2/external/opampack/.../lib/getopt/getopt.ml`) with a standalone
//! OCaml program (`$HOME/OpenCode/tmp/w2b/getopt_oracle/`). It established
//! the `clap` forms asserted below: `-n5` (value attached to a short
//! option), `--ndet=5` (long option through `=`), `--` (end of options),
//! `""` (anonymous argument), and `-n -x` (a next argv element starting
//! with `-` is never consumed as a value).
//!
//! Wave S2.1 (text release): the hand-written help, the `EZFIO_DIR is
//! missing`/`Inconsistent command line` messages and the print-then-
//! continue `-h` behavior are gone. `clap` now owns the positional count
//! (its own diagnostics, exit 2) and the help flag (`DisplayHelp`, exit
//! 0, no continuation into the argument checks). The value-conversion
//! messages (`[-n|--ndet] expects an integer`, `[-s|--state] expects a
//! range`) are unchanged.

use clap::Parser;
use clap::error::ErrorKind;
use qp_types::command_line::{
    EXIT_FAILURE, EXIT_USAGE, EditArgs, EditConfig, RunArgs, RunConfig, UsageError,
    parse_state_range,
};

fn edit(argv: &[&str]) -> Result<EditConfig, UsageError> {
    EditArgs::try_parse_from(argv)
        .map_err(|error| UsageError::new(error.to_string()))?
        .validate()
}

fn run(argv: &[&str]) -> Result<RunConfig, UsageError> {
    RunArgs::try_parse_from(argv)
        .map(RunArgs::validate)
        .map_err(|error| UsageError::new(error.to_string()))
}

#[test]
fn edit_messages_and_option_spellings() {
    // A missing or doubled `EZFIO_DIR` is now `clap`'s own required-
    // argument error (exit 2), not `EZFIO_DIR is missing`.
    let missing = EditArgs::try_parse_from(["qp_edit"]).unwrap_err();
    assert_eq!(missing.kind(), ErrorKind::MissingRequiredArgument);
    assert_eq!(missing.exit_code(), EXIT_USAGE);
    assert!(
        missing.to_string().contains("<EZFIO_DIR>"),
        "{}",
        missing.to_string()
    );

    let extra = edit(&["qp_edit", "dir", "extra"]).expect_err("extra positional");
    assert!(
        extra
            .message()
            .contains("unexpected argument 'extra' found"),
        "{}",
        extra.message()
    );
    assert_eq!(
        EditArgs::try_parse_from(["qp_edit", "dir", "extra"])
            .unwrap_err()
            .kind(),
        ErrorKind::UnknownArgument
    );

    // The value-conversion messages are unchanged from the reference.
    assert_eq!(
        edit(&["qp_edit", "-n", "x", "dir"]),
        Err(UsageError::new("[-n|--ndet] expects an integer"))
    );
    assert_eq!(
        edit(&["qp_edit", "-s", "x", "dir"]),
        Err(UsageError::new("[-s|--state] expects a range"))
    );
    // `-s ""`: `Range.of_string` raises `Invalid_argument` on `s.[0]`, and
    // qp_edit's `with _` turns it into the range message.
    assert_eq!(
        edit(&["qp_edit", "-s", "", "dir"]),
        Err(UsageError::new("[-s|--state] expects a range"))
    );

    // All three spellings of `-n 5` and the `-s` list.
    for argv in [
        &["qp_edit", "-n", "5", "dir"][..],
        &["qp_edit", "-n5", "dir"][..],
        &["qp_edit", "--ndet=5", "dir"][..],
    ] {
        let config = edit(argv).expect("ndet spelling");
        assert_eq!(config.ndet, Some(5), "{argv:?}");
        assert_eq!(config.ezfio_dir, "dir");
    }
    let config = edit(&["qp_edit", "-c", "-s", "[1,3-5]", "dir"]).expect("state list");
    assert!(config.check);
    assert_eq!(config.state, Some(vec![1, 3, 4, 5]));
    assert_eq!(config.write, None);
    assert_eq!(config.read, None);
    // Options after the anonymous argument.
    assert_eq!(
        edit(&["qp_edit", "dir", "-n", "5"])
            .expect("option after dir")
            .ndet,
        Some(5)
    );
    // `--` protects a dash-leading directory name.
    assert_eq!(
        edit(&["qp_edit", "--", "-dir"]).expect("-- dir").ezfio_dir,
        "-dir"
    );
    // An empty string is an anonymous argument, not "missing".
    assert_eq!(edit(&["qp_edit", ""]).expect("empty dir").ezfio_dir, "");
    // Write and read are both accepted here; their incompatibility is a
    // `run`-time failure, not a command-line one.
    let config = edit(&["qp_edit", "-w", "a", "-r", "b", "dir"]).expect("read+write");
    assert_eq!(config.write.as_deref(), Some("a"));
    assert_eq!(config.read.as_deref(), Some("b"));
}

#[test]
fn run_matches_the_reference() {
    // A wrong positional count is now `clap`'s required-argument error
    // (exit 2), not `Inconsistent command line`.
    let missing_program = RunArgs::try_parse_from(["qp_run"]).unwrap_err();
    assert_eq!(missing_program.kind(), ErrorKind::MissingRequiredArgument);
    assert_eq!(missing_program.exit_code(), EXIT_USAGE);
    assert!(
        missing_program.to_string().contains("<PROGRAM>"),
        "{}",
        missing_program.to_string()
    );

    let missing_dir = run(&["qp_run", "scf"]).expect_err("missing EZFIO_DIR");
    assert!(
        missing_dir.message().contains("<EZFIO_DIR>"),
        "{}",
        missing_dir.message()
    );
    let extra = run(&["qp_run", "scf", "dir", "extra"]).expect_err("extra positional");
    assert!(
        extra
            .message()
            .contains("unexpected argument 'extra' found"),
        "{}",
        extra.message()
    );

    assert_eq!(
        run(&["qp_run", "scf", "dir"]).expect("two arguments"),
        RunConfig {
            slave: false,
            prefix: None,
            program: "scf".to_owned(),
            ezfio_dir: "dir".to_owned(),
        }
    );
    let config = run(&["qp_run", "-s", "-p", "gdb", "scf", "dir"]).expect("slave + prefix");
    assert!(config.slave);
    assert_eq!(config.prefix.as_deref(), Some("gdb"));
    assert_eq!(
        run(&["qp_run", "--prefix=valgrind", "scf", "dir"])
            .expect("long prefix")
            .prefix
            .as_deref(),
        Some("valgrind")
    );
}

#[test]
fn help_is_the_clap_built_in_and_stops_the_program() {
    // `-h`/`--help` is `clap`'s `DisplayHelp`: the binary prints the help
    // once and exits 0. The reference printed its own help and *continued*
    // into the argument checks (most tools then printed the help a second
    // time and exited 2); the continuation is a deliberate S2 break.
    for (label, error) in [
        (
            "edit -h",
            EditArgs::try_parse_from(["qp_edit", "-h"]).unwrap_err(),
        ),
        (
            "edit --help",
            EditArgs::try_parse_from(["qp_edit", "--help"]).unwrap_err(),
        ),
        (
            "run -h",
            RunArgs::try_parse_from(["qp_run", "-h"]).unwrap_err(),
        ),
        (
            "run --help",
            RunArgs::try_parse_from(["qp_run", "--help"]).unwrap_err(),
        ),
    ] {
        assert_eq!(error.kind(), ErrorKind::DisplayHelp, "{label}");
        assert_eq!(error.exit_code(), 0, "{label}");
    }
}

#[test]
fn usage_errors_exit_2_as_measured() {
    assert_eq!(EXIT_USAGE, 2);
    assert_eq!(EXIT_FAILURE, 1);

    // Unknown option: the reference raises an uncaught
    // `Getopt.Error("Unknown option -Z")`, exit 2.
    let error = EditArgs::try_parse_from(["qp_edit", "-Z", "dir"]).expect_err("-Z must fail");
    assert_eq!(error.exit_code(), EXIT_USAGE);
    assert_eq!(error.kind(), ErrorKind::UnknownArgument);
    // A value-taking option at the end of argv: the reference raises
    // `Getopt.Error("Option ndet must have an argument")`, exit 2.
    let error = EditArgs::try_parse_from(["qp_edit", "-n"]).expect_err("-n needs a value");
    assert_eq!(error.exit_code(), EXIT_USAGE);
    // `--` ends option parsing, so a dash-leading directory name parses.
    assert_eq!(
        run(&["qp_run", "scf", "--", "-dir"])
            .expect("-- dir")
            .ezfio_dir,
        "-dir"
    );
}

#[test]
fn state_range_matches_range_of_string() {
    assert_eq!(parse_state_range("7"), Some(vec![7]));
    assert_eq!(parse_state_range("007"), Some(vec![7]));
    assert_eq!(parse_state_range("[1-3,5]"), Some(vec![1, 2, 3, 5]));
    assert_eq!(
        parse_state_range("[72-73,36-37]"),
        Some(vec![36, 37, 72, 73])
    );
    assert_eq!(parse_state_range("[1-2,2-3]"), Some(vec![1, 2, 3]));
    assert_eq!(parse_state_range("[]"), Some(vec![]));
    assert_eq!(parse_state_range("[1,,2]"), Some(vec![1, 2]));
    assert_eq!(parse_state_range("[1-3]"), Some(vec![1, 2, 3]));
    // Failures: `parse_state_range` returns `None` and the caller emits the
    // single `[-s|--state] expects a range` message.
    assert_eq!(parse_state_range("[3-1]"), None);
    assert_eq!(parse_state_range("x"), None);
    assert_eq!(parse_state_range(""), None);
    assert_eq!(parse_state_range("[1"), None);
    assert_eq!(parse_state_range("1]"), None);
    assert_eq!(parse_state_range("1-3"), None);
    assert_eq!(parse_state_range("[-2]"), None);
    assert_eq!(parse_state_range("[1-3-5]"), None);
    // OCaml `int_of_string` also accepts `0x`/`0o`/`0b` and `_` separators;
    // this port's decimal-only subset is documented in the module notes.
    assert_eq!(parse_state_range("0x10"), None);
    assert_eq!(parse_state_range("1_0"), None);
}
