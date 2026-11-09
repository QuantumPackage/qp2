//! End-to-end matrix for `qp_set_mo_class` (wave S2: clap + the shared
//! `ToolError`).
//!
//! Every failure case keeps the reference's scenario (same argv, same
//! fixture mutation) and asserts the exit code (2) plus the one new
//! `qp_set_mo_class: <message>` line on standard error; the success cases
//! keep the functional stdout byte for byte. Write cases additionally
//! check the decompressed `mo_basis/mo_class` payload through
//! `qp_ezfio` (the `.gz` bytes are not comparable across gzip
//! implementations; the decompressed payload is).
//!
//! Fixtures are copied from `rust/tests/golden/fixtures` into per-test
//! temp dirs and never mutated in place. `QP_ROOT` is a per-test fake
//! root (the `src/bitmask/bitmasks_module.f90` file is created where the
//! `-q` fallback path forces `Qpackage.bit_kind_size`).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

#[cfg(unix)]
use std::os::unix::process::CommandExt;

use qp_ezfio::EzFio;

const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/golden/fixtures");

/// The tool name of the one-line standard-error error.
const TOOL: &str = "qp_set_mo_class";

/// The golden `-q` stdout on the `hcn` fixture with the `bit_kind_size = 64`
/// fallback (the fixture has no `determinants/n_int`).
const GOLDEN_Q: &str = "Electrons: 7 7\n\
                         MO  : 55\n\
                         n_int: 1\n\
                         Core\n\
                         Inactive\n\
                         Active   : (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55)\n\
                         Virtual\n\
                         Deleted\n";

/// The `bit_kind_size` line of the reference's bitmask module (64/16).
const BK64: &str = "  integer, parameter             :: bit_kind_size = 64\n";
const BK16: &str = "  integer, parameter             :: bit_kind_size = 16\n";

// -------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------

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
    // The usage line of the help is checked through `argv[0]`; `arg0`
    // keeps it the bare program name.
    #[cfg(unix)]
    command.arg0(TOOL);
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
    command.arg0(TOOL);
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

/// Every failure: exit 2, and the single line
/// `qp_set_mo_class: <message>` on standard error (the message is the
/// stable string named by each test).
fn assert_fatal(outcome: &Outcome, message: &str) {
    assert_eq!(outcome.exit, 2, "stderr: {}", outcome.stderr);
    assert_eq!(
        outcome.stderr,
        format!("{TOOL}: {message}\n"),
        "the one standard-error line"
    );
}

/// Exit 2 with a clap diagnostic containing `needle` (the stable string
/// named by the test).
fn assert_clap(outcome: &Outcome, needle: &str) {
    assert_eq!(outcome.exit, 2, "stderr: {}", outcome.stderr);
    assert!(
        outcome.stderr.contains(needle),
        "stderr:\n{} was expected to contain {needle:?}",
        outcome.stderr
    );
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
        format!("{TOOL}: QP_ROOT is not set (source quantum_package.rc)\n")
    );
}

// -------------------------------------------------------------------
// CLI: clap grammar, the help flag and the anonymous argument
// -------------------------------------------------------------------

#[test]
fn no_argument_fails_with_the_clap_usage() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &[]);
    assert_clap(
        &outcome,
        "the following required arguments were not provided",
    );
    assert_eq!(outcome.stdout, b""[..]);
}

#[test]
fn the_help_flag_prints_help_and_exits_zero() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    for flag in ["-h", "--help"] {
        let outcome = run(&root, work.path(), &[flag]);
        assert_eq!(outcome.exit, 0, "flag {flag}");
        assert!(
            outcome
                .stdout
                .windows(b"Usage:".len())
                .any(|w| w == b"Usage:"),
            "flag {flag}; stdout: {}",
            String::from_utf8_lossy(&outcome.stdout)
        );
        assert_eq!(outcome.stderr, "", "flag {flag}");
    }
}

#[test]
fn two_positionals_are_rejected() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["a", "b"]);
    assert_clap(&outcome, "unexpected argument 'b' found");
}

#[test]
fn unknown_options_are_rejected() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-z", "ezfio"]);
    assert_clap(&outcome, "unexpected argument '-z' found");
    let outcome = run(&root, work.path(), &["--nosuch", "ezfio"]);
    assert_clap(&outcome, "unexpected argument '--nosuch' found");
}

#[test]
fn a_missing_value_names_the_option() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c"]);
    assert_clap(&outcome, "a value is required for '--core");
}

#[test]
fn an_empty_range_is_rejected() {
    // `-c ""`: clap passes the empty value through; the range grammar
    // rejects it before any `Number of MOs:` line.
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "", "ezfio"]);
    assert_eq!(outcome.stdout, b""[..]);
    assert_fatal(&outcome, "invalid range: empty string");
}

#[test]
fn short_option_forms_follow_the_clap_grammar() {
    // clap: `-c=1` and `-c3` both give the value `1`/`3` (the reference
    // had no `=` syntax for short options).
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c=1", "ezfio"]);
    assert_fatal(&outcome, "Orbital 2 is not specified (mo_num = 55)");
    let outcome = run(&root, work.path(), &["-c3", "ezfio"]);
    assert_fatal(&outcome, "Orbital 1 is not specified (mo_num = 55)");
    let outcome = run(&root, work.path(), &["--core=5", "ezfio"]);
    assert_fatal(&outcome, "Orbital 1 is not specified (mo_num = 55)");
    let outcome = run(&root, work.path(), &["-qz", "ezfio"]);
    assert_clap(&outcome, "unexpected argument '-z' found");
}

#[test]
fn the_qc_cluster_is_accepted_by_clap() {
    // The reference's Getopt rejected `-qc` (a value-taking option in a
    // cluster); clap accepts it, takes the next token as the `-c` value
    // and runs the query (the value is then ignored in `-q` mode).
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-qc", "1", "ezfio"]);
    assert_eq!(outcome.exit, 0, "stderr: {}", outcome.stderr);
    assert_eq!(outcome.stdout, GOLDEN_Q.as_bytes());
}

#[test]
fn double_dash_switches_everything_to_anonymous() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    // `-- ezfio`: one positional, the run proceeds with all classes
    // defaulting to `[]`.
    let outcome = run(&root, work.path(), &["--", "ezfio"]);
    assert_eq!(outcome.stdout, b""[..]);
    assert_fatal(&outcome, "Orbital 1 is not specified (mo_num = 55)");
}

// -------------------------------------------------------------------
// Range grammar and `MO_number` validation (set mode, `hcn` = 55 MOs)
// -------------------------------------------------------------------

#[test]
fn range_error_classes() {
    let cases = [
        ("x", "invalid range 'x': expected a digit or '['"),
        ("]", "invalid range ']': expected a digit or '['"),
        ("[1", "invalid range '[1': missing ']'"),
        ("[3-1]", "range 3-1 is reversed"),
        ("0", "orbital number 0 must be at least 1"),
        ("2a", "invalid integer '2a' in range '2a'"),
        // `-1` splits at its `-` into `("", "-1")`: the empty start piece
        // fails first, as in the reference.
        ("[-1]", "invalid integer '' in range '[-1]'"),
        ("[1-2,a]", "invalid integer 'a' in range '[1-2,a]'"),
        ("[a,1-2]", "invalid integer 'a' in range '[a,1-2]'"),
    ];
    for (range_value, message) in cases {
        let work = workdir();
        ezfio(work.path(), "hcn");
        let root = fake_root(work.path(), None);
        let outcome = run(&root, work.path(), &["-c", range_value, "ezfio"]);
        assert_eq!(outcome.stdout, b""[..], "range value {range_value}");
        assert_fatal(&outcome, message);
    }
}

#[test]
fn an_unspecified_orbital_is_reported_in_order() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[]", "ezfio"]);
    assert_fatal(&outcome, "Orbital 1 is not specified (mo_num = 55)");
    let outcome = run(&root, work.path(), &["-c", "[1,2,]", "ezfio"]);
    assert_fatal(&outcome, "Orbital 3 is not specified (mo_num = 55)");
}

#[test]
fn the_maximum_check_prints_the_mo_count_first() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[56]", "ezfio"]);
    assert_eq!(outcome.stdout, b"Number of MOs: 55\n");
    assert_fatal(&outcome, "orbital 56 is above mo_num = 55");
}

#[test]
fn a_duplicate_orbital_names_both_spaces() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "-a", "[1]", "ezfio"]);
    assert_fatal(
        &outcome,
        "Orbital 1 is defined both in the Core and Active spaces",
    );
    let outcome = run(&root, work.path(), &["-c", "[1-55]", "-a", "[55]", "ezfio"]);
    assert_fatal(
        &outcome,
        "Orbital 55 is defined both in the Core and Active spaces",
    );
}

#[test]
fn the_ten_thousand_warning_goes_to_standard_error() {
    // `mo_num = 10001`: `MO_number.of_int(10001)` warns (standard error,
    // wave S2) before the `not specified` failure.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", &format!("{:20}", 10001));
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[10001]", "ezfio"]);
    assert_eq!(outcome.stdout, b""[..]);
    assert_eq!(
        outcome.stderr,
        format!("More than 10_000 MOs\n{TOOL}: Orbital 1 is not specified (mo_num = 10001)\n")
    );
}

#[test]
fn mo_num_one_skips_the_maximum_check() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", &format!("{:20}", 1));
    let root = fake_root(work.path(), None);
    // `max = 1` skips the check; `apply_class` then indexes a
    // one-element array.
    let outcome = run(&root, work.path(), &["-c", "[2]", "ezfio"]);
    assert_eq!(outcome.stdout, b""[..]);
    assert_fatal(&outcome, "orbital index out of bounds: 2");
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
fn a_negative_mo_num_is_rejected() {
    // The reference's `Array.init mo_num` on a negative count.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", "-3");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "ezfio"]);
    assert_fatal(&outcome, "mo_num must be non-negative: -3");
}

#[test]
fn a_non_numeric_mo_num_is_rejected() {
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", "abc\n");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "ezfio"]);
    assert_fatal(&outcome, "invalid integer 'abc' in mo_basis/mo_num");
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
    assert_fatal(&outcome, "mo_basis/mo_num not found");
}

#[test]
fn a_missing_directory_is_created_then_guarded() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "nosuch"]);
    assert_fatal(&outcome, "mo_basis/mo_num not found");
    // Side effect of `Ezfio.set_file`: the directory and the `ezfio`
    // group were created before the guard failed.
    assert!(work.path().join("nosuch/.version").exists());
    assert!(work.path().join("nosuch/ezfio").is_dir());
    assert!(work.path().join("nosuch/ezfio/creation").exists());
    assert!(work.path().join("nosuch/ezfio/user").exists());
    assert!(work.path().join("nosuch/ezfio/library").exists());
}

#[test]
fn an_existing_regular_file_cannot_be_created() {
    let work = workdir();
    fs::write(work.path().join("ezfio"), b"").unwrap();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-c", "[1]", "ezfio"]);
    assert_fatal(&outcome, "cannot create ezfio: File exists (os error 17)");
}

#[test]
fn a_missing_parent_cannot_be_created() {
    let work = workdir();
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["ezfio/sub"]);
    assert_fatal(
        &outcome,
        "cannot create ezfio/sub: No such file or directory (os error 2)",
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
    // class per MO (the stored bytes become the bare keyword).
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
    // The run re-writes the array (bare keywords).
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
    // `mo_num = 9000`, no `n_int`, `bk = 64`: fallback `n_int = 141`.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", &format!("{:20}", 9000));
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_eq!(outcome.exit, 2);
    assert_eq!(outcome.stdout, b"N_int > 128");
    assert_fatal(&outcome, "n_int 141 is above the maximum of 128");
}

#[test]
fn q_with_a_zero_or_negative_n_int_re_fails_the_maximum() {
    // `n_int = 0` (or negative): the direct check takes the fallback,
    // the fallback's `max` re-read is the bad value, and `x = 1 > max`.
    for value in ["0", "-5"] {
        let work = workdir();
        ezfio(work.path(), "hcn");
        set_scalar(work.path(), "determinants/n_int", &format!("{:20}", value));
        let root = fake_root(work.path(), Some(BK64));
        let outcome = run(&root, work.path(), &["-q", "ezfio"]);
        assert_eq!(outcome.stdout, b""[..], "n_int {value}");
        assert_fatal(
            &outcome,
            &format!("n_int 1 is above the maximum of {value}"),
        );
    }
}

#[test]
fn q_with_a_garbage_n_int_reraises() {
    // The fallback's `max` default re-reads the file: the bad value is
    // outside the caught `try`.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "determinants/n_int", "abc\n");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "invalid integer 'abc' in determinants/n_int");
}

#[test]
fn q_fallback_errors_from_the_bitmask_file() {
    // Bad value: the `Bit_kind_size` validation.
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(
        work.path(),
        Some("  integer, parameter :: bit_kind_size = 42\n"),
    );
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "Bit_kind_size should be (8|16|32|64).");
    // Missing file.
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), None);
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    let expected = format!(
        "File {} not found",
        root.join("src/bitmask/bitmasks_module.f90").display()
    );
    assert_fatal(&outcome, &expected);
    // Present but no matching line.
    let work = workdir();
    ezfio(work.path(), "hcn");
    let root = fake_root(work.path(), Some("  integer :: other = 16\n"));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    let expected = format!(
        "bit_kind_size not found in {}",
        root.join("src/bitmask/bitmasks_module.f90").display()
    );
    assert_fatal(&outcome, &expected);
}

#[test]
fn q_electrons_edge_cases() {
    // Missing `elec_alpha_num`: the section is `None`.
    let work = workdir();
    ezfio(work.path(), "hcn");
    remove(work.path(), "electrons");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "Unable to read electrons");

    // A present-but-empty alpha file.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "electrons/elec_alpha_num", "");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(
        &outcome,
        "cannot read ezfio/electrons/elec_alpha_num: the file is empty",
    );

    // A garbage alpha file.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "electrons/elec_alpha_num", "abc\n");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(
        &outcome,
        "invalid integer 'abc' in electrons/elec_alpha_num",
    );

    // `elec_alpha_num = 0`: the positivity check.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(
        work.path(),
        "electrons/elec_alpha_num",
        &format!("{:20}", 0),
    );
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "alpha electron count 0 must be positive");

    // A missing beta file: the read is not `has`-guarded.
    let work = workdir();
    ezfio(work.path(), "hcn");
    remove(work.path(), "electrons/elec_beta_num");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(
        &outcome,
        "cannot read ezfio/electrons/elec_beta_num: No such file or directory (os error 2)",
    );

    // A negative beta: the `x >= 0` check.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(
        work.path(),
        "electrons/elec_beta_num",
        &format!("{:20}", -1),
    );
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "beta electron count -1 must be non-negative");

    // Beta zero is valid.
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
    // `mo_num = 0` with 7 alpha electrons: the shared handler's guard.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", &format!("{:20}", 0));
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "More alpha electrons than MOs");

    // A non-numeric `mo_num`: the driver's own scalar read classifies it.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_num", "abc\n");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_fatal(&outcome, "invalid integer 'abc' in mo_basis/mo_num");

    // A rejected `mo_label` spelling: the echo moved to standard error
    // (wave S2), then the failure.
    let work = workdir();
    ezfio(work.path(), "hcn");
    set_scalar(work.path(), "mo_basis/mo_label", "bogus\n");
    let root = fake_root(work.path(), Some(BK64));
    let outcome = run(&root, work.path(), &["-q", "ezfio"]);
    assert_eq!(outcome.exit, 2);
    assert_eq!(outcome.stdout, b""[..], "the echo no longer goes to stdout");
    assert!(
        outcome.stderr.starts_with("bogus\n"),
        "stderr:\n{}",
        outcome.stderr
    );
    assert!(
        outcome.stderr.contains("MO_label should be one of:"),
        "stderr:\n{}",
        outcome.stderr
    );
}

#[test]
fn q_with_a_missing_mo_class_materializes_the_default() {
    // The reference's read materialises `Active` per MO (and rewrites the
    // file), then `set` re-writes it.
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
    // Each `MO_number.of_int` conversion inside the shared `Mo_basis`
    // reader warns (stdout, kept) *twice*, and the `n_int` fallback
    // (157 with `bk = 64`) warns again (`N_int > 128`, stdout, kept)
    // before its maximum-check failure.
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
    assert_fatal(&outcome, "n_int 157 is above the maximum of 128");
}
