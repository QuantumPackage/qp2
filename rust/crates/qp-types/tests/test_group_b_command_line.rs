//! CLI parity tests for the three `Command_line` specifications that the
//! executables parse with `clap`.
//!
//! The expected messages are transcribed from the reference sources:
//! `ocaml/qp_create_ezfio.ml:700–825`, `scripts/ezfio_interface/qp_edit_template`
//! (the generated `qp_edit`), `ocaml/qp_run.ml:232–270`, and
//! `ocaml/Command_line.ml:161–199`. `qp_set_mo_class` and `qp_tunnel` are not
//! here: their binaries parse with a statement-by-statement `getopt` port
//! (P4) and pin that parser against renamed reference copies instead.
//!
//! The parsing rules were checked against the vendored `getopt`
//! (`~/qp2/external/opampack/.../lib/getopt/getopt.ml`) with a standalone
//! OCaml program (`$HOME/OpenCode/tmp/w2b/getopt_oracle/`). It established the
//! forms asserted below: `-bcc-pvdz` (value attached to a short option),
//! `-axq` (clustered argument-less shorts), `--basis=` (empty value through
//! `=`), `--` (end of options), `-` and `""` (anonymous arguments), and
//! `-c -1` (`Option -c must have an argument`: a next argv element starting
//! with `-` is never consumed as a value). Where `clap` cannot reproduce a
//! `getopt` behaviour, the deviation is asserted and named.

use clap::Parser;
use clap::error::ErrorKind;
use qp_types::command_line::{
    CreateArgs, CreateConfig, EXIT_FAILURE, EXIT_USAGE, EditArgs, EditConfig, RunArgs, RunConfig,
    UsageError, parse_state_range,
};

fn create(argv: &[&str]) -> Result<CreateConfig, UsageError> {
    CreateArgs::try_parse_from(argv)
        .map_err(|error| UsageError::new(error.to_string()))?
        .validate()
}

fn edit(argv: &[&str]) -> Result<EditConfig, UsageError> {
    EditArgs::try_parse_from(argv)
        .map_err(|error| UsageError::new(error.to_string()))?
        .validate()
}

fn run(argv: &[&str]) -> Result<RunConfig, UsageError> {
    RunArgs::try_parse_from(argv)
        .map_err(|error| UsageError::new(error.to_string()))?
        .validate()
}

/// `clap` exit code for a raw parse of `argv`, as the executable would get it.
fn clap_exit(argv: &[&str]) -> i32 {
    match CreateArgs::try_parse_from(argv) {
        Ok(_) => 0,
        Err(error) => error.exit_code(),
    }
}

#[test]
fn create_messages_and_values_match_the_reference() {
    assert_eq!(
        create(&["qp_create_ezfio"]),
        Err(UsageError::new("--basis option is missing."))
    );
    assert_eq!(
        create(&["qp_create_ezfio", "-b", "cc-pvdz"]),
        Err(UsageError::new("input file is missing"))
    );
    assert_eq!(
        create(&["qp_create_ezfio", "-b", "x", "-c", "x", "f"]),
        Err(UsageError::new("int_of_string"))
    );
    assert_eq!(
        create(&["qp_create_ezfio", "-b", "x", "-g", "x", "f"]),
        Err(UsageError::new("float_of_string"))
    );
    assert_eq!(
        create(&["qp_create_ezfio", "-b", "x", "-m", "x", "f"]),
        Err(UsageError::new("int_of_string"))
    );

    // Every option at once; the second anonymous argument is ignored
    // (`anon_args ()` matches `x :: _`).
    let config = create(&[
        "qp_create_ezfio",
        "-o",
        "out",
        "--basis",
        "cc-pvdz",
        "-a",
        "--charge",
        "m1",
        "-g",
        "1.5",
        "-m",
        "3",
        "-p",
        "bfd",
        "-q",
        "-x",
        "mol.xyz",
        "ignored.xyz",
    ])
    .expect("valid create command line");
    assert_eq!(
        config,
        CreateConfig {
            output: Some("out".to_owned()),
            basis: "cc-pvdz".to_owned(),
            au: true,
            charge: -1,
            ghost: 1.5,
            multiplicity: 3,
            pseudo: Some("bfd".to_owned()),
            qmc: true,
            cartesian: true,
            file: "mol.xyz".to_owned(),
        }
    );

    // Defaults: charge 0, ghost 0.0, multiplicity 1, no pseudo/output.
    let defaults = create(&["qp_create_ezfio", "-b", "x", "f"]).expect("defaults");
    assert_eq!(defaults.output, None);
    assert_eq!(defaults.charge, 0);
    assert_eq!(defaults.ghost, 0.0);
    assert_eq!(defaults.multiplicity, 1);
    assert_eq!(defaults.pseudo, None);
    assert!(!defaults.au && !defaults.qmc && !defaults.cartesian);
}

#[test]
fn create_accepts_the_getopt_option_spellings() {
    // `-bcc-pvdz`: value attached to the short option.
    assert_eq!(
        create(&["qp_create_ezfio", "-bcc-pvdz", "f"])
            .expect("attached value")
            .basis,
        "cc-pvdz"
    );
    // `--basis=cc-pvdz`: empty/attached value through `=`.
    assert_eq!(
        create(&["qp_create_ezfio", "--basis=cc-pvdz", "f"])
            .expect("long =value")
            .basis,
        "cc-pvdz"
    );
    // `--basis=`: the `=` form accepts an empty value (see the deviation test
    // for `-b ""`, which crashes the reference getopt).
    assert_eq!(
        create(&["qp_create_ezfio", "--basis=", "f"])
            .expect("empty basis")
            .basis,
        ""
    );
    // `-axq`: argument-less shorts concatenate.
    let flags = create(&["qp_create_ezfio", "-b", "x", "-axq", "f"]).expect("clustered");
    assert!(flags.au && flags.cartesian && flags.qmc);
    // Options may follow the anonymous argument.
    assert_eq!(
        create(&["qp_create_ezfio", "f", "-b", "x"])
            .expect("option after file")
            .file,
        "f"
    );
    // `--` ends option parsing, so a leading dash is a file name.
    assert_eq!(
        create(&["qp_create_ezfio", "-b", "x", "--", "-weird.xyz"])
            .expect("-- file")
            .file,
        "-weird.xyz"
    );
    // A lone `-` is an anonymous argument, not an option.
    assert_eq!(
        create(&["qp_create_ezfio", "-b", "x", "-"])
            .expect("dash file")
            .file,
        "-"
    );
    // Charge syntax: `m` is the minus sign, `m0` is zero.
    for (text, expected) in [("0", 0), ("7", 7), ("m1", -1), ("m0", 0), ("m13", -13)] {
        assert_eq!(
            create(&["qp_create_ezfio", "-b", "x", "-c", text, "f"])
                .expect("charge form")
                .charge,
            expected,
            "charge {text}"
        );
    }
}

#[test]
fn create_help_skips_the_mandatory_option_check() {
    // `qp_create_ezfio -h` prints the help, `set_specs` skips the mandatory
    // `-b` check, and the anonymous check then fails: CONTRACTS §2.4.
    let args = CreateArgs::try_parse_from(["qp_create_ezfio", "-h"]).expect("-h parses");
    assert!(args.help);
    assert_eq!(
        args.validate(),
        Err(UsageError::new("input file is missing"))
    );
    // With a file the tool would go on and run with an empty basis.
    let config = create(&["qp_create_ezfio", "-h", "mol.xyz"]).expect("help + file");
    assert_eq!(config.basis, "");
    assert_eq!(config.file, "mol.xyz");
    // `--help` is the same option.
    assert!(
        CreateArgs::try_parse_from(["qp_create_ezfio", "--help"])
            .expect("--help")
            .help
    );
}

#[test]
fn edit_messages_and_option_spellings() {
    assert_eq!(
        edit(&["qp_edit"]),
        Err(UsageError::new("EZFIO_DIR is missing"))
    );
    assert_eq!(
        edit(&["qp_edit", "dir", "extra"]),
        Err(UsageError::new("EZFIO_DIR is missing"))
    );
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
    assert_eq!(
        run(&["qp_run"]),
        Err(UsageError::new("Inconsistent command line"))
    );
    assert_eq!(
        run(&["qp_run", "scf"]),
        Err(UsageError::new("Inconsistent command line"))
    );
    assert_eq!(
        run(&["qp_run", "scf", "dir", "extra"]),
        Err(UsageError::new("Inconsistent command line"))
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
fn every_tool_recognises_help_and_still_checks_anonymous_arguments() {
    // `set_specs` prepends `-h` to every spec list; the anonymous check still
    // runs, which is why the reference prints the help twice and exits 2.
    assert!(
        EditArgs::try_parse_from(["qp_edit", "-h"])
            .expect("-h")
            .help
    );
    assert!(
        RunArgs::try_parse_from(["qp_run", "--help"])
            .expect("-h")
            .help
    );
    assert_eq!(
        edit(&["qp_edit", "-h"]),
        Err(UsageError::new("EZFIO_DIR is missing"))
    );
    assert_eq!(
        run(&["qp_run", "-h"]),
        Err(UsageError::new("Inconsistent command line"))
    );
}

#[test]
fn usage_errors_exit_2_as_measured() {
    assert_eq!(EXIT_USAGE, 2);
    assert_eq!(EXIT_FAILURE, 1);

    // Unknown long option: the reference raises an uncaught
    // `Getopt.Error("Unknown option --nope")`, exit 2.
    let error =
        CreateArgs::try_parse_from(["qp_create_ezfio", "--nope"]).expect_err("--nope must fail");
    assert_eq!(error.exit_code(), EXIT_USAGE);
    assert_eq!(error.kind(), ErrorKind::UnknownArgument);
    let error = EditArgs::try_parse_from(["qp_edit", "-Z", "dir"]).expect_err("-Z must fail");
    assert_eq!(error.exit_code(), EXIT_USAGE);
    // A value-taking option at the end of argv: the reference raises
    // `Getopt.Error("Option basis must have an argument")`, exit 2.
    let error = CreateArgs::try_parse_from(["qp_create_ezfio", "--basis"])
        .expect_err("--basis needs a value");
    assert_eq!(error.exit_code(), EXIT_USAGE);
    let error = EditArgs::try_parse_from(["qp_edit", "-n"]).expect_err("-n needs a value");
    assert_eq!(error.exit_code(), EXIT_USAGE);
    // `-c -1` is rejected by both parsers (getopt: "Option -c must have an
    // argument"); only the text differs.
    let error = CreateArgs::try_parse_from(["qp_create_ezfio", "-b", "x", "-c", "-1", "f"])
        .expect_err("-c -1 must fail");
    assert_eq!(error.exit_code(), EXIT_USAGE);
    assert_eq!(clap_exit(&["qp_create_ezfio", "--", "-b"]), 0);
}

#[test]
fn getopt_deviations_are_the_documented_ones() {
    // (1) The reference `getopt` evaluates `args.(no+1).[0]`, so an empty
    // *next* argv element after a value-taking option raises
    // `Invalid_argument("index out of bounds")` before any tool code runs.
    // `clap` accepts it as an empty value; this is the port's deviation.
    let config = create(&["qp_create_ezfio", "-b", "", "mol.xyz"]).expect("empty -b value");
    assert_eq!(config.basis, "");
    assert_eq!(config.file, "mol.xyz");
    // (2) Empty `--charge=`: `getopt` stores "", then qp_create_ezfio's
    // `x.[0]` raises `Invalid_argument` (uncaught). The port reports the
    // `int_of_string` message instead; same exit code 2.
    assert_eq!(
        create(&["qp_create_ezfio", "-b", "x", "--charge=", "f"]),
        Err(UsageError::new("int_of_string"))
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
