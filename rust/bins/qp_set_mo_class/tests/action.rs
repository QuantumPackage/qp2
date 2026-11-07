//! End-to-end matrix for `qp_set_mo_class`, pinned against the reference
//! binary's measured outputs (2026-11-05; scratch
//! `$HOME/OpenCode/tmp/w3-p4-tools/set-mo-class/runs*`).
//!
//! Every case asserts the exit code, the exact stdout bytes and the exact
//! stderr line. Write cases additionally check the decompressed
//! `mo_basis/mo_class` payload through `qp_ezfio` (the `.gz` bytes are not
//! comparable across gzip implementations; the decompressed payload is).
//!
//! Fixtures are copied from `rust/tests/golden/fixtures` into per-test
//! temp dirs and never mutated in place. `QP_ROOT` is a per-test fake root
//! (the `src/bitmask/bitmasks_module.f90` file is created where the `-q`
//! fallback path forces `Qpackage.bit_kind_size`).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

#[cfg(unix)]
use std::os::unix::process::CommandExt;

use qp_ezfio::EzFio;

const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/golden/fixtures");

/// The reference help capture for `argv[0] = "qp_set_mo_class"` (639 B).
const HELP: &str = include_str!("data/help_bare.txt");

/// The golden `-q` stdout on the `hcn` fixture with the `bit_kind_size = 64`
/// fallback (the fixture has no `determinants/n_int`): measured `runs/q-hcn`.
const GOLDEN_Q: &str = "Electrons: 7 7\n\
                         MO  : 55\n\
                         n_int: 1\n\
                         Core\n\
                         Inactive\n\
                         Active   : (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55)\n\
                         Virtual\n\
                         Deleted\n";

/// The reference's assertion renderings (reference build locations).
const RANGE_REVERSED: &str = "File \"Range.ml\", line 26, characters 8-14: Assertion failed";
const RANGE_FIRST: &str = "File \"Range.ml\", line 45, characters 4-10: Assertion failed";
const RANGE_LAST: &str = "File \"Range.ml\", line 46, characters 4-10: Assertion failed";
const MO_NUM_MIN: &str = "File \"Qptypes.ml\", line 471, characters 6-12: Assertion failed";
const ELEC_ALPHA_MIN: &str = "File \"Qptypes.ml\", line 83, characters 20-26: Assertion failed";
const ELEC_BETA_MIN: &str = "File \"Qptypes.ml\", line 70, characters 20-26: Assertion failed";

/// The `bit_kind_size` line of the reference's bitmask module (64/16).
const BK64: &str = "  integer, parameter             :: bit_kind_size = 64\n";
const BK16: &str = "  integer, parameter             :: bit_kind_size = 16\n";

// -------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------

fn fatal(message: &str) -> String {
    format!("Fatal error: exception {message}\n")
}

fn fixture(name: &str) -> PathBuf {
    Path::new(FIXTURES).join(name)
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), &target).unwrap();
        }
    }
}

/// A fresh temp dir that becomes the process cwd of the child.
fn workdir() -> tempfile::TempDir {
    tempfile::tempdir().unwrap()
}

/// Copy `fixtures/<name>` to `<workdir>/ezfio`.
fn ezfio(work: &Path, name: &str) {
    copy_dir(&fixture(name), &work.join("ezfio"));
}

/// The fake `QP_ROOT` under `<work>/qproot`; `bitmask` is written to
/// `src/bitmask/bitmasks_module.f90` when `Some`.
fn fake_root(work: &Path, bitmask: Option<&str>) -> PathBuf {
    let root = work.join("qproot");
    fs::create_dir_all(root.join("data")).unwrap();
    if let Some(body) = bitmask {
        let file = root.join("src/bitmask/bitmasks_module.f90");
        fs::create_dir_all(file.parent().unwrap()).unwrap();
        fs::write(&file, body).unwrap();
    }
    root
}

struct Outcome {
    exit: i32,
    stdout: Vec<u8>,
    stderr: String,
}

fn run(root: &Path, work: &Path, args: &[&str]) -> Outcome {
    let mut command = Command::new(env!("CARGO_BIN_EXE_qp_set_mo_class"));
    command.env("QP_ROOT", root).current_dir(work).args(args);
    // The reference captures were taken with the bare program name
    // (invocation through PATH); `arg0` reproduces that `argv[0]`.
    #[cfg(unix)]
    command.arg0("qp_set_mo_class");
    let output = command
        .output()
        .expect("the compiled binary must be runnable");
    Outcome {
        exit: output.status.code().unwrap_or(-1),
        stdout: output.stdout,
        stderr: String::from_utf8(output.stderr).expect("UTF-8 stderr"),
    }
}

/// Run with `QP_ROOT` removed from the environment.
fn run_no_root(work: &Path, args: &[&str]) -> Outcome {
    let mut command = Command::new(env!("CARGO_BIN_EXE_qp_set_mo_class"));
    command.env_remove("QP_ROOT").current_dir(work).args(args);
    #[cfg(unix)]
    command.arg0("qp_set_mo_class");
    let output = command
        .output()
        .expect("the compiled binary must be runnable");
    Outcome {
        exit: output.status.code().unwrap_or(-1),
        stdout: output.stdout,
        stderr: String::from_utf8(output.stderr).expect("UTF-8 stderr"),
    }
}

fn set_scalar(work: &Path, rel: &str, text: &str) {
    let path = work.join("ezfio").join(rel);
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    fs::write(path, text).unwrap();
}

fn remove(work: &Path, rel: &str) {
    let path = work.join("ezfio").join(rel);
    if path.is_dir() {
        fs::remove_dir_all(path).unwrap();
    } else {
        fs::remove_file(path).unwrap();
    }
}

/// The decompressed `mo_basis/mo_class` values, in order.
fn payload(work: &Path) -> Vec<String> {
    let ezfio = EzFio::open(work.join("ezfio")).expect("the run produced an EZFIO dir");
    let array = ezfio
        .get_array::<String>("mo_basis", "mo_class")
        .expect("mo_class was written");
    array.data().to_vec()
}

fn dims(work: &Path) -> Vec<usize> {
    let ezfio = EzFio::open(work.join("ezfio")).expect("the run produced an EZFIO dir");
    let array = ezfio
        .get_array::<String>("mo_basis", "mo_class")
        .expect("mo_class was written");
    array.dims().to_vec()
}

fn assert_fatal(outcome: &Outcome, exception: &str) {
    let expected = format!("Fatal error: exception {exception}\n");
    assert_eq!(outcome.exit, 2, "exit code (stderr: {})", outcome.stderr);
    assert_eq!(outcome.stderr, expected);
}

// -------------------------------------------------------------------
// The QP_ROOT guard (module-initialiser timing: before any parse output)
// -------------------------------------------------------------------

#[test]
fn qp_root_unset_fails_before_anything() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let outcome = run_no_root(work.path(), &[]);
    assert_eq!(outcome.exit, 2);
    assert!(outcome.stdout.is_empty(), "no help, no output");
    assert_eq!(
        outcome.stderr,
        fatal(
            "Failure(\"QP_ROOT environment variable is not set.\\nPlease source the quantum_package.rc file.\")"
        )
    );
}

// -------------------------------------------------------------------
// CLI: parse errors and the help / anonymous-argument flow
// -------------------------------------------------------------------

#[test]
fn no_argument_prints_help_then_fails() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &[]);
    assert_eq!(outcome.exit, 2);
    assert_eq!(outcome.stdout, HELP.as_bytes());
    assert_fatal(&outcome, "Failure(\"EZFIO_FILE is missing\")");
}

#[test]
fn h_flag_prints_help_twice() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-h"]);
    assert_eq!(outcome.exit, 2);
    assert_eq!(outcome.stdout, format!("{HELP}{HELP}").as_bytes());
    assert_fatal(&outcome, "Failure(\"EZFIO_FILE is missing\")");
}

#[test]
fn two_anonymous_arguments_fail_with_help() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["a", "b"]);
    assert_eq!(outcome.exit, 2);
    assert_eq!(outcome.stdout, HELP.as_bytes());
    assert_fatal(&outcome, "Failure(\"EZFIO_FILE is missing\")");
}

#[test]
fn h_with_a_directory_continues_the_run() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-h", "ezfio"]);
    assert_eq!(outcome.exit, 2);
    assert_eq!(outcome.stdout, HELP.as_bytes());
    assert_fatal(
        &outcome,
        "Failure(\"Orbital 1 is not specified (mo_num = 55)\")",
    );
}

#[test]
fn unknown_short_option_is_a_getopt_error() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-z", "ezfio"]);
    assert_eq!(outcome.stdout, b""[..]);
    assert_fatal(&outcome, "Getopt.Error(\"Unknown option -z\")");
}

#[test]
fn unknown_long_option_is_a_getopt_error() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["--nosuch", "ezfio"]);
    assert_fatal(&outcome, "Getopt.Error(\"Unknown option --nosuch\")");
}

#[test]
fn a_parse_error_suppresses_the_help() {
    // `-h` was seen, but the parse fails before the `set_specs` help
    // decision (measured: no help in stdout).
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-h", "-z", "ezfio"]);
    assert_eq!(outcome.stdout, b""[..]);
    assert_fatal(&outcome, "Getopt.Error(\"Unknown option -z\")");
}

#[test]
fn missing_value_messages_name_the_option_form() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    for (args, message) in [
        (vec!["-c"], "Option -c must have an argument"),
        (vec!["--core"], "Option --core must have an argument"),
        (vec!["-c", "-q", "ezfio"], "Option -c must have an argument"),
    ] {
        let outcome = run(&root, work.path(), &args);
        assert_fatal(&outcome, &format!("Getopt.Error(\"{message}\")"));
    }
}

#[test]
fn an_empty_next_argument_raises_index_out_of_bounds() {
    // `-c ""`: the reference's `args.(no+1).[0]` lookahead indexes the
    // empty token (measured; no `Number of MOs:` line — the failure is
    // before the `check` wrapper).
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "", "ezfio"]);
    assert_eq!(outcome.stdout, b""[..]);
    assert_fatal(&outcome, "Invalid_argument(\"index out of bounds\")");
}

#[test]
fn only_no_argument_options_can_be_concatenated() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-qc", "1", "ezfio"]);
    assert_fatal(
        &outcome,
        "Getopt.Error(\"Only non-argument short-options can be concatenated (error with option c in -qc)\")",
    );
    // An unknown character in the tail is the lookup's own error.
    let outcome = run(&root, work.path(), &["-qz", "ezfio"]);
    assert_fatal(&outcome, "Getopt.Error(\"Unknown option -z\")");
}

#[test]
fn long_options_reject_values_and_take_eq_values() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["--help=1"]);
    assert_fatal(
        &outcome,
        "Getopt.Error(\"Option --help does not accept argument\")",
    );
}

#[test]
fn double_dash_switches_everything_to_anonymous() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    // `-- ezfio`: one anonymous argument, the run proceeds with all
    // classes defaulting to `[]` (measured `runs/cli-dashdash`).
    let outcome = run(&root, work.path(), &["--", "ezfio"]);
    assert_eq!(outcome.stdout, b""[..]);
    assert_fatal(
        &outcome,
        "Failure(\"Orbital 1 is not specified (mo_num = 55)\")",
    );
}

#[test]
fn short_options_have_no_eq_syntax() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    // `-c=1`: the value is literally `=1` → not a digit, not `[` → the
    // first-character assertion (measured `runs2/c-eq`).
    let outcome = run(&root, work.path(), &["-c=1", "ezfio"]);
    assert_fatal(&outcome, RANGE_FIRST);
    // `-c3`: the value is `3` → orbital 3 classified, 1 is not.
    let outcome = run(&root, work.path(), &["-c3", "ezfio"]);
    assert_fatal(
        &outcome,
        "Failure(\"Orbital 1 is not specified (mo_num = 55)\")",
    );
    // `--core=5`: the value is `5` (the long form has `=` syntax).
    let outcome = run(&root, work.path(), &["--core=5", "ezfio"]);
    assert_fatal(
        &outcome,
        "Failure(\"Orbital 1 is not specified (mo_num = 55)\")",
    );
}

// -------------------------------------------------------------------
// Range grammar and `MO_number` validation (set mode, `hcn` = 55 MOs)
// -------------------------------------------------------------------

#[test]
fn range_error_classes() {
    let cases = [
        ("x", RANGE_FIRST),
        ("]", RANGE_FIRST),
        ("[1", RANGE_LAST),
        ("[3-1]", RANGE_REVERSED),
        ("0", MO_NUM_MIN),
        ("2a", "Failure(\"int_of_string\")"),
        ("[-1]", "Failure(\"int_of_string\")"),
        ("[1-2,a]", "Failure(\"int_of_string\")"),
        ("[a,1-2]", "Failure(\"int_of_string\")"),
    ];
    for (range_value, stderr) in cases {
        let work = workdir();
        ezfio(work.path(), "hcn");
        let root = fake_root(work.path(), None);
        let outcome = run(&root, work.path(), &["-c", range_value, "ezfio"]);
        assert_eq!(outcome.stdout, b""[..], "range value {range_value}");
        assert_fatal(&outcome, stderr);
    }
}

#[test]
fn an_unspecified_orbital_is_reported_in_order() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[]", "ezfio"]);
    assert_fatal(
        &outcome,
        "Failure(\"Orbital 1 is not specified (mo_num = 55)\")",
    );
    let outcome = run(&root, work.path(), &["-c", "[1,2,]", "ezfio"]);
    assert_fatal(
        &outcome,
        "Failure(\"Orbital 3 is not specified (mo_num = 55)\")",
    );
}

#[test]
fn the_maximum_check_prints_the_mo_count_first() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[56]", "ezfio"]);
    assert_eq!(outcome.stdout, b"Number of MOs: 55\n");
    assert_fatal(&outcome, "Invalid_argument(\"MO_number: 56\")");
}

#[test]
fn a_duplicate_orbital_names_both_spaces() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "-a", "[1]", "ezfio"]);
    assert_fatal(
        &outcome,
        "Failure(\"Orbital 1 is defined both in the Core and Active spaces\")",
    );
    let outcome = run(&root, work.path(), &["-c", "[1-55]", "-a", "[55]", "ezfio"]);
    assert_fatal(
        &outcome,
        "Failure(\"Orbital 55 is defined both in the Core and Active spaces\")",
    );
}

#[test]
fn the_ten_thousand_warning_is_raw_stdout() {
    // `mo_num = 10001`: `MO_number.of_int(10001)` prints the warning
    // (no newline) before the `not specified` failure (measured
    // `runs/e-warn`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", &format!("{:20}", 10001));
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[10001]", "ezfio"]);
    assert_eq!(outcome.stdout, b"More than 10_000 MOs");
    assert_fatal(
        &outcome,
        "Failure(\"Orbital 1 is not specified (mo_num = 10001)\")",
    );
}

#[test]
fn mo_num_one_skips_the_maximum_check() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", &format!("{:20}", 1));
    let root = fake_root(work.path(), None);
    // `max = 1` skips the check; `apply_class` then indexes a
    // one-element array (measured `runs/e-max1`).
    let outcome = run(&root, work.path(), &["-c", "[2]", "ezfio"]);
    assert_eq!(outcome.stdout, b""[..]);
    assert_fatal(&outcome, "Invalid_argument(\"index out of bounds\")");
    // The valid case classifies the single MO.
    let outcome = run(&root, work.path(), &["-c", "[1]", "ezfio"]);
    assert_eq!(outcome.exit, 0);
    assert_eq!(
        outcome.stdout,
        b"Core     : (1)\nInactive\nActive\nVirtual\nDeleted\n"
    );
    assert_eq!(outcome.stderr, "");
    assert_eq!(payload(work.path()), vec!["Core".to_owned()]);
}

#[test]
fn mo_num_zero_writes_an_empty_array() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", &format!("{:20}", 0));
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["ezfio"]);
    assert_eq!(outcome.exit, 0);
    assert_eq!(
        outcome.stdout,
        b"Core\nInactive\nActive\nVirtual\nDeleted\n"
    );
    assert_eq!(outcome.stderr, "");
    assert_eq!(payload(work.path()), Vec::<String>::new());
    assert_eq!(dims(work.path()), vec![0]);
}

#[test]
fn a_negative_mo_num_is_array_init() {
    // `Array.init mo_num` on a negative count (measured).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", "-3");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "ezfio"]);
    assert_fatal(&outcome, "Invalid_argument(\"Array.init\")");
}

#[test]
fn a_non_numeric_mo_num_is_int_of_string() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", "abc\n");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "ezfio"]);
    assert_fatal(&outcome, "Failure(\"int_of_string\")");
}

// -------------------------------------------------------------------
// `set_file` (the tree open/create path)
// -------------------------------------------------------------------

#[test]
fn a_tree_without_mo_num_fails_the_guard() {
    let work = workdir();
    ezfio(work.path(), "h2");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "ezfio"]);
    assert_fatal(&outcome, "Failure(\"mo_basis/mo_num not found\")");
}

#[test]
fn a_missing_directory_is_created_then_guarded() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "nosuch"]);
    assert_fatal(&outcome, "Failure(\"mo_basis/mo_num not found\")");
    // Side effect of `Ezfio.set_file`: the directory and the `ezfio`
    // group were created before the guard failed.
    assert!(work.path().join("nosuch/.version").exists());
    assert!(work.path().join("nosuch/ezfio").is_dir());
    assert!(work.path().join("nosuch/ezfio/creation").exists());
    assert!(work.path().join("nosuch/ezfio/user").exists());
    assert!(work.path().join("nosuch/ezfio/library").exists());
}

#[test]
fn an_existing_regular_file_is_mkdir_eexist() {
    let work = workdir();
    fs::write(work.path().join("ezfio"), b"").unwrap();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "ezfio"]);
    assert_fatal(
        &outcome,
        "Unix.Unix_error(Unix.EEXIST, \"mkdir\", \"ezfio\")",
    );
}

#[test]
fn a_missing_parent_is_mkdir_enoent() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["ezfio/sub"]);
    assert_fatal(
        &outcome,
        "Unix.Unix_error(Unix.ENOENT, \"mkdir\", \"ezfio/sub\")",
    );
}

// -------------------------------------------------------------------
// A valid partition: the debug lines and the written array
// -------------------------------------------------------------------

#[test]
fn a_valid_partition_writes_the_bare_keywords() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(
        &root,
        work.path(),
        &[
            "-c", "[1]", "-i", "[2,3]", "-a", "[4-6]", "-v", "[7-55]", "ezfio",
        ],
    );
    assert_eq!(outcome.exit, 0);
    assert_eq!(
        outcome.stdout,
        format!(
            "Core     : (1)\n\
             Inactive : (2, 3)\n\
             Active   : (4, 5, 6)\n\
             Virtual  : ({})\n\
             Deleted\n",
            (7..=55)
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )
        .as_bytes()
    );
    assert_eq!(outcome.stderr, "");
    let mut expected = vec![
        "Core".to_owned(),
        "Inactive".to_owned(),
        "Inactive".to_owned(),
    ];
    expected.extend(std::iter::repeat_n("Active".to_owned(), 3));
    expected.extend(std::iter::repeat_n("Virtual".to_owned(), 49));
    assert_eq!(payload(work.path()), expected);
}

#[test]
fn re_specifying_the_current_classes_keeps_the_tree_content() {
    // `hcn` is all `Active`: re-specifying that partition writes the same
    // class per MO (the stored bytes become the bare keyword, as the
    // reference's `write_string_array` does — measured `runs/v-current`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-a", "[1-55]", "ezfio"]);
    assert_eq!(outcome.exit, 0);
    assert_eq!(
        outcome.stdout,
        format!(
            "Core\nInactive\nActive   : ({})\nVirtual\nDeleted\n",
            (1..=55)
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )
        .as_bytes()
    );
    assert_eq!(payload(work.path()), vec!["Active".to_owned(); 55]);
}

// -------------------------------------------------------------------
// The `-q` path
// -------------------------------------------------------------------

#[test]
fn q_on_hcn_is_the_golden_stdout() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_eq!(outcome.exit, 0);
    assert_eq!(outcome.stdout, GOLDEN_Q.as_bytes());
    assert_eq!(outcome.stderr, "");
    // The run re-writes the array (bare keywords, as the reference).
    assert_eq!(payload(work.path()), vec!["Active".to_owned(); 55]);
}

#[test]
fn q_reads_determinants_n_int_when_present() {
    for (value, line) in [("3", "n_int: 3"), ("200", "n_int: 200")] {
        let work = workdir();
        ezfio(work.path(), "hcn");
        set_scalar(work.path(), "determinants/n_int", &format!("{:20}", value));
        let root = fake_root(work.path(), Some(BK64));
        let outcome = run(&root, work.path(), &["-q", "ezfio"]);
        assert_eq!(outcome.exit, 0, "n_int {value}");
        let expected = format!(
            "{}Electrons: 7 7\nMO  : 55\n{line}\n",
            if value == "200" { "N_int > 128" } else { "" }
        );
        assert!(
            outcome.stdout.starts_with(expected.as_bytes()),
            "stdout starts with the summary (n_int {value}): {}",
            String::from_utf8_lossy(&outcome.stdout)
        );
    }
}

#[test]
fn q_without_n_int_uses_the_bit_kind_size_fallback() {
    // The fixture has no `determinants/n_int`: `n_int = (55-1)/bk + 1`.
    for (bitmask, line) in [(Some(BK64), "n_int: 1"), (Some(BK16), "n_int: 4")] {
        let work = workdir();
        ezfio(work.path(), "hcn");
        let root = fake_root(work.path(), bitmask);
        let outcome = run(&root, work.path(), &["-q", "ezfio"]);
        assert_eq!(outcome.exit, 0);
        assert!(
            outcome
                .stdout
                .windows(line.len())
                .any(|w| w == line.as_bytes()),
            "stdout contains {line}: {}",
            String::from_utf8_lossy(&outcome.stdout)
        );
    }
}

#[test]
fn q_with_a_large_mo_num_warns_and_fails_the_maximum() {
    // `mo_num = 9000`, no `n_int`, `bk = 64`: fallback `n_int = 141`
    // (measured `runs3/qmo9000`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", &format!("{:20}", 9000));
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_eq!(outcome.exit, 2);
    assert_eq!(outcome.stdout, b"N_int > 128");
    assert_fatal(&outcome, "Invalid_argument(\"N_int_number: 141\")");
}

#[test]
fn q_with_a_zero_or_negative_n_int_re_fails_the_maximum() {
    // `n_int = 0` (or negative): the direct assert is caught by the
    // `try`, the fallback's `max` re-read is the bad value, and
    // `x = 1 > max` (measured `runs4/nint0`, `runs4/nintneg`).
    for value in ["0", "-5"] {
        let work = workdir();
        ezfio(work.path(), "hcn");
        set_scalar(work.path(), "determinants/n_int", &format!("{:20}", value));
        let root = fake_root(work.path(), Some(BK64));
        let outcome = run(&root, work.path(), &["-q", "ezfio"]);
        assert_eq!(outcome.exit, 2, "n_int {value}");
        assert_eq!(outcome.stdout, b""[..]);
        assert_fatal(&outcome, "Invalid_argument(\"N_int_number: 1\")");
    }
}

#[test]
fn q_with_a_garbage_n_int_reraises_uncaught() {
    // The fallback's `max` default re-reads the file: `int_of_string`
    // there is outside the `try` (measured `runs4/nintabc`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "determinants/n_int", "abc\n");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "Failure(\"int_of_string\")");
}

#[test]
fn q_fallback_errors_from_the_bitmask_file() {
    // Bad value: the `Bit_kind_size` validation (measured `runs7/bkbad2`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(
        work.path(),
        Some("  integer, parameter :: bit_kind_size = 42\n"),
    );
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(
        &outcome,
        "Invalid_argument(\"Bit_kind_size should be (8|16|32|64).\")",
    );
    // Missing file (measured `runs7/bkmissing`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    let expected = format!(
        "Failure(\"File {} not found\")",
        root.join("src/bitmask/bitmasks_module.f90").display()
    );
    assert_fatal(&outcome, &expected);
    // Present but no matching line (measured `runs7/bknoline`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), Some("  integer :: other = 16\n"));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    let expected = format!(
        "Failure(\"bit_kind_size not found in {}\")",
        root.join("src/bitmask/bitmasks_module.f90").display()
    );
    assert_fatal(&outcome, &expected);
}

#[test]
fn q_electrons_edge_cases() {
    // Missing `elec_alpha_num`: the section is `None` (measured
    // `runs6/ealphamissing`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    remove(work.path(), "electrons");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "Failure(\"Unable to read electrons\")");

    // A present-but-empty alpha file: `input_line` → `End_of_file`
    // (measured `runs8/eempty`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "electrons/elec_alpha_num", "");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "End_of_file");

    // A garbage alpha file: `int_of_string` (measured `runs8/eabc`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "electrons/elec_alpha_num", "abc\n");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "Failure(\"int_of_string\")");

    // `elec_alpha_num = 0`: the positivity assert (measured
    // `runs6/ealpha0`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(
        work.path(),
        "electrons/elec_alpha_num",
        &format!("{:20}", 0),
    );
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, ELEC_ALPHA_MIN);

    // A missing beta file: `open_in` on the not `has`-guarded key
    // (measured `runs6/ebetamissing`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    remove(work.path(), "electrons/elec_beta_num");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(
        &outcome,
        "Sys_error(\"ezfio/electrons/elec_beta_num: No such file or directory\")",
    );

    // A negative beta: the `x >= 0` assert (measured `runs6/ebetaneg`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(
        work.path(),
        "electrons/elec_beta_num",
        &format!("{:20}", -1),
    );
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, ELEC_BETA_MIN);

    // Beta zero is valid (measured `runs6/ebetazero`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "electrons/elec_beta_num", &format!("{:20}", 0));
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_eq!(outcome.exit, 0);
    assert!(
        outcome
            .stdout
            .starts_with(b"Electrons: 7 0\nMO  : 55\nn_int: 1\n")
    );
}

#[test]
fn q_mo_edge_cases() {
    // `mo_num = 0` with 7 alpha electrons: the shared handler's guard
    // (measured `runs6/qmo0`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", &format!("{:20}", 0));
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "Failure(\"More alpha electrons than MOs\")");

    // A non-numeric `mo_num`: the driver's own scalar read classifies it
    // (measured `runs3/qbadnum`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", "abc\n");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "Failure(\"int_of_string\")");

    // A rejected `mo_label` spelling is echoed to stdout, then fails
    // (measured `runs6/qlabelbad`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_label", "bogus\n");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_eq!(outcome.stdout, b"bogus\n");
    assert_fatal(
        &outcome,
        "Failure(\"MO_label should be one of:\\nGuess | Orthonormalized | Canonical | Natural | Localized | MCSCF | MOM | None.\")",
    );
}

#[test]
fn q_with_a_missing_mo_class_materializes_the_default() {
    // The reference's read materialises `Active` per MO (and rewrites the
    // file), then `set` re-writes it (measured `runs3/qnoclass`).
    let work = workdir();
    ezfio(work.path(), "hcn");
    remove(work.path(), "mo_basis/mo_class.gz");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_eq!(outcome.exit, 0);
    assert_eq!(outcome.stdout, GOLDEN_Q.as_bytes());
    assert_eq!(payload(work.path()), vec!["Active".to_owned(); 55]);
}

#[test]
fn q_with_mo_num_over_ten_thousand_warns_per_conversion() {
    // Each `MO_number.of_int` conversion inside `Mo_basis.read` warns
    // (measured: the warning is printed *twice*), and the `n_int`
    // fallback (157 with `bk = 64`) warns again before its
    // maximum-check failure.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", &format!("{:20}", 10001));
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_eq!(outcome.exit, 2);
    assert_eq!(
        outcome.stdout,
        b"More than 10_000 MOsMore than 10_000 MOsN_int > 128"
    );
    assert_fatal(&outcome, "Invalid_argument(\"N_int_number: 157\")");
}
