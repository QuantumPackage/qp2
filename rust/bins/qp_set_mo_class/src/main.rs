//! `qp_set_mo_class` — Rust port of `ocaml/qp_set_mo_class.ml` (P4, wave 3).
//!
//! The reference driver (`let ()` at the end of the file) has **no** `with`
//! handler: every exception escapes to the OCaml runtime, which prints
//! `Fatal error: exception <exception>` on stderr and exits 2. The flow is
//! the reference's statement order:
//!
//! 1. the `QP_ROOT` guard, raised by a module initialiser in the
//!    reference (before the option parser runs — even a bare invocation
//!    dies here, measured);
//! 2. the `Getopt`/`Command_line` parse ([`cli`], a statement-by-statement
//!    port including its quirks), the `-h` help, and the mandatory
//!    anonymous `EZFIO_DIR`;
//! 3. `Ezfio.set_file` (create the tree when absent: `mkdir`, the `ezfio`
//!    group, the `creation`/`user`/`library` stamps);
//! 4. `has mo_basis/mo_num`, then [`set_with`] or [`get`] (the `-q` path:
//!    electrons, MOs, the `n_int` resolution, then `set` re-run on the
//!    current classes).
//!
//! `build_mask` (`qp_set_mo_class.ml:9-28`) is dead code in the reference
//! and is not ported.

mod cli;
mod error;
mod ocaml_num;
mod qpackage;
mod range;

use std::fs;
use std::io::{self, Write};
use std::path::Path;
use std::process::ExitCode;

use qp_ezfio::{Array, EzFio};
use qp_input::handlers;
use qp_input::handlers::mo_basis::MoClass;

use crate::error::{AppError, EXIT_CODE, OcamlException, QP_ROOT_MISSING};

/// The five class keywords, in the reference's print order
/// (core, inact, act, virt, del) — the bare `MO_class.to_string` forms.
const KEYWORDS: [&str; 5] = ["Core", "Inactive", "Active", "Virtual", "Deleted"];

/// The padded debug labels, from `MO_class.to_string`
/// (`ocaml/MO_class.ml:15-30`).
const LABELS: [&str; 5] = [
    "Core     ",
    "Inactive ",
    "Active   ",
    "Virtual  ",
    "Deleted  ",
];

fn main() -> ExitCode {
    let mut stdout = io::stdout();
    let result = run(&mut stdout);
    // The reference writes help/warnings to stdout and the fatal line to
    // stderr; flush stdout first so the order holds when both go to one
    // terminal.
    let _ = stdout.flush();
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            let _ = io::stderr().write_all(error.render().as_bytes());
            ExitCode::from(EXIT_CODE)
        }
    }
}

/// The reference's top-level flow: guard, parse, act.
fn run(out: &mut dyn Write) -> Result<(), AppError> {
    // Module-initialiser semantics: the failure precedes the parse and any
    // output (measured: bare invocation, no help, empty stdout).
    if std::env::var_os("QP_ROOT").is_none() {
        return Err(AppError(OcamlException::Failure(
            QP_ROOT_MISSING.to_owned(),
        )));
    }
    let argv = collect_argv()?;
    let parsed = cli::parse(&argv)?;
    parsed.after_specs(&argv[0], out)?;
    let dir = cli::require_one_anonymous(&parsed, &argv[0], out)?;
    run_action(&dir, &parsed)
}

/// Collect `argv` (including `argv[0]`) as owned `String`s.
///
/// The reference sees raw bytes. This port's parser is `String`-based, so
/// a non-UTF-8 argument is rejected instead of panicking; that message is
/// a deviation with no reference counterpart (argument 0 is included, so
/// the index is the `argv` index).
fn collect_argv() -> Result<Vec<String>, AppError> {
    std::env::args_os()
        .enumerate()
        .map(|(index, argument)| {
            argument.into_string().map_err(|_| {
                AppError(OcamlException::SysError(format!(
                    "command line argument {index} is not valid UTF-8"
                )))
            })
        })
        .collect()
}

/// `run` (`qp_set_mo_class.ml:176-185`): `set_file`, the `mo_num` guard,
/// then `get` or `set`.
fn run_action(dir: &str, p: &cli::Parsed) -> Result<(), AppError> {
    let ezfio = set_file(dir)?;
    if !ezfio.has("mo_basis", "mo_num").map_err(app)? {
        return Err(failure("mo_basis/mo_num not found"));
    }
    if p.query {
        get(&ezfio)
    } else {
        set_with(
            &ezfio,
            p.core.as_deref().unwrap_or("[]"),
            p.inact.as_deref().unwrap_or("[]"),
            p.act.as_deref().unwrap_or("[]"),
            p.virt.as_deref().unwrap_or("[]"),
            p.del.as_deref().unwrap_or("[]"),
        )
    }
}

/// A registry that knows no keys.
///
/// The reference `ezfio.ml` client validates nothing against a declared
/// table: its readers check only the per-line syntax (the OCaml
/// `int_of_string`/`float_of_string`/string parse), and `write_array`
/// writes whatever it is given. The shape/type specs of `StaticRegistry`
/// are a Rust-side defence for the array-heavy `qp_create_ezfio`; this
/// binary opens the tree without them so an inconsistent tree (e.g.
/// `mo_num` changed after the arrays were written) reads and writes
/// exactly as the reference does (measured: `mo_num = 9000` over the
/// `hcn` 55-MO arrays proceeds to the `n_int` fallback).
#[derive(Debug)]
struct NoRegistry;

impl qp_ezfio::Registry for NoRegistry {
    fn lookup(&self, _group: &str, _key: &str) -> Option<qp_ezfio::KeySpec> {
        None
    }

    fn specs(&self) -> &'static [qp_ezfio::KeySpec] {
        &[]
    }
}

/// `Ezfio.set_file` (`ezfio.ml:426-439`): open the tree when
/// `<dir>/.version` exists; otherwise `mkdir` it (single level), `mkdir`
/// the `ezfio` group, and stamp it.
fn set_file(dir: &str) -> Result<EzFio, AppError> {
    let no_registry = || Box::new(NoRegistry) as Box<dyn qp_ezfio::Registry>;
    let root = Path::new(dir);
    if root.join(".version").exists() {
        return EzFio::open(root)
            .map_err(app)
            .map(|ezfio| ezfio.with_registry(no_registry()));
    }
    // The reference's `mkdir` is a single-level `Unix.mkdir`: a missing
    // parent is `ENOENT`, an existing path (file, or a directory without
    // the `.version` marker) is `EEXIST`.
    let mkdir_error = |errno: error::UnixErrno| {
        AppError(OcamlException::UnixError {
            errno,
            function: "mkdir",
            argument: dir.to_owned(),
        })
    };
    if let Some(parent) = root.parent()
        && !parent.as_os_str().is_empty()
        && !parent.exists()
    {
        return Err(mkdir_error(error::UnixErrno::Enoent));
    }
    if root.exists() {
        return Err(mkdir_error(error::UnixErrno::Eexist));
    }
    let ezfio = EzFio::create(root)
        .map_err(|_| mkdir_error(error::UnixErrno::Eexist))?
        .with_registry(no_registry());
    // The reference leaves a partial tree on a later failure (no cleanup).
    if let Err(error) = ezfio.create_group("ezfio") {
        return Err(app(error));
    }
    write_stamps(root)?;
    Ok(ezfio)
}

/// The `ezfio/{creation,user,library}` stamps of `Ezfio.set_file`
/// (`ezfio.ml:431-434`), the `qp_create_ezfio` pattern.
fn write_stamps(ezfio_file: &Path) -> Result<(), AppError> {
    let failure = || {
        AppError(OcamlException::Failure(format!(
            "Unable to create new ezfio file:\n{}",
            ezfio_file.display()
        )))
    };
    let dir = ezfio_file.join("ezfio");
    // `LANG= date`: the local time in the C locale.
    let date = std::process::Command::new("date")
        .env_remove("LANG")
        .output()
        .map_err(|_| failure())?;
    if !date.status.success() {
        return Err(failure());
    }
    // `echo $USER`: empty line when the variable is unset.
    let user = std::env::var("USER").unwrap_or_default();
    // The reference stamps its configure-time ezfio library path
    // (`ezfio.ml:28`); a static pure-Rust build has no such path, and
    // nothing in the tool layer reads the file back — the built-in
    // implementation's name stands in (documented deviation, T4.1).
    const LIBRARY: &str = "qp-ezfio";
    if fs::write(dir.join("creation"), &date.stdout).is_err()
        || fs::write(dir.join("user"), format!("{user}\n").as_bytes()).is_err()
        || fs::write(dir.join("library"), format!("{LIBRARY}\n").as_bytes()).is_err()
    {
        return Err(failure());
    }
    Ok(())
}

/// `set` (`qp_set_mo_class.ml:36-112`): read `mo_num` raw, parse the five
/// ranges in order (first failure wins, through the `check` wrapper),
/// apply them (duplicate → `Failure`), require full coverage, print the
/// five debug lines, write `mo_basis/mo_class`.
fn set_with(
    ezfio: &EzFio,
    core: &str,
    inact: &str,
    act: &str,
    virt: &str,
    del: &str,
) -> Result<(), AppError> {
    let mo_num = read_int_scalar(ezfio, "mo_basis", "mo_num")?;
    // `Array.init mo_num`: a negative count is the uncaught
    // `Invalid_argument("Array.init")` before any range is read.
    if mo_num < 0 {
        return Err(AppError(OcamlException::InvalidArgument(
            "Array.init".into(),
        )));
    }
    let mo_num_us = mo_num as usize;

    // Parse in the reference's order; `check` (inside `range::create_class`)
    // catches only `Invalid_argument`.
    let core = range::create_class(core, mo_num)?;
    let inact = range::create_class(inact, mo_num)?;
    let act = range::create_class(act, mo_num)?;
    let virt = range::create_class(virt, mo_num)?;
    let del = range::create_class(del, mo_num)?;

    // `apply_class` in the same order. The stored value is the *bare*
    // class (the reference stores the empty-list form).
    let mut assigned: Vec<Option<u8>> = vec![None; mo_num_us];
    for (kind, list) in [(0u8, &core), (1, &inact), (2, &act), (3, &virt), (4, &del)] {
        for &i in list {
            // `mo_class.(i-1)`: array indexing is `Invalid_argument`
            // out of bounds. Reachable only when `mo_num = 1` (the
            // `max = 1` maximum-check skip) — measured.
            if i < 1 || i as usize > mo_num_us {
                return Err(AppError(OcamlException::InvalidArgument(
                    "index out of bounds".into(),
                )));
            }
            let idx = i as usize - 1;
            match assigned[idx] {
                None => assigned[idx] = Some(kind),
                Some(prev) => {
                    return Err(failure(format!(
                        "Orbital {i} is defined both in the {} and {} spaces",
                        KEYWORDS[prev as usize], KEYWORDS[kind as usize]
                    )));
                }
            }
        }
    }

    for (idx, slot) in assigned.iter().enumerate() {
        if slot.is_none() {
            return Err(failure(format!(
                "Orbital {} is not specified (mo_num = {mo_num})",
                idx + 1,
            )));
        }
    }

    // Debug output, in the reference's order, full `MO_class.to_string`.
    let mut out = io::stdout();
    for (kind, list) in [(0u8, &core), (1, &inact), (2, &act), (3, &virt), (4, &del)] {
        let line = if list.is_empty() {
            KEYWORDS[kind as usize].to_owned()
        } else {
            format!(
                "{label}: ({numbers})",
                label = LABELS[kind as usize],
                numbers = list
                    .iter()
                    .map(i64::to_string)
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        };
        writeln!(out, "{line}").map_err(io_error)?;
    }

    let data: Vec<String> = assigned
        .iter()
        .map(|slot| match slot {
            Some(kind) => KEYWORDS[*kind as usize].to_owned(),
            None => unreachable!("full coverage checked above"),
        })
        .collect();
    let array = Array::new(vec![mo_num_us], data).map_err(app)?;
    ezfio
        .set_array::<String>("mo_basis", "mo_class", &array)
        .map_err(app)
}

/// `get` (`qp_set_mo_class.ml:116-172`, the `-q` path): read the electrons
/// and the MOs, resolve `n_int`, print the three summary lines, then call
/// `set` on the current classes (which re-parses, re-prints and re-writes
/// `mo_basis/mo_class`).
fn get(ezfio: &EzFio) -> Result<(), AppError> {
    // `Input.Electrons.read`: `None` only when `elec_alpha_num` is absent;
    // a read failure of the *present* alpha propagates, as does the beta
    // read (which is not `has`-guarded).
    if !ezfio.has("electrons", "elec_alpha_num").map_err(app)? {
        return Err(failure("Unable to read electrons"));
    }
    let alpha = read_int_scalar(ezfio, "electrons", "elec_alpha_num")?;
    let beta = read_int_scalar(ezfio, "electrons", "elec_beta_num")?;
    // `Elec_alpha_number.of_int` / `Elec_beta_number.of_int` asserts.
    if alpha <= 0 {
        return Err(AppError(OcamlException::ELEC_ALPHA_MIN));
    }
    if beta < 0 {
        return Err(AppError(OcamlException::ELEC_BETA_MIN));
    }

    // `Input.Mo_basis.read`: classify the `mo_num` read ourselves (the
    // reference's `read_int` failures are uncaught here), then let the
    // shared handler do the reference's field reads — the `mo_label` and
    // `mo_class` materialisation, the `More alpha electrons than MOs`
    // guard, the `MO_number` conversion (its `More than 10_000 MOs`
    // warning included), the `ao_md5` check.
    if !ezfio.has("mo_basis", "mo_num").map_err(app)? {
        return Err(failure("Unable to read MOs"));
    }
    read_int_scalar(ezfio, "mo_basis", "mo_num")?;
    let mo = match handlers::mo_basis::read(ezfio) {
        // Unreachable after the `has` guard above; the reference's
        // `failwith "Unable to read MOs"` covers it.
        Ok(None) => return Err(failure("Unable to read MOs")),
        Ok(Some(mo)) => mo,
        Err(error) => return Err(AppError(OcamlException::Failure(error.to_string()))),
    };
    let mo_num = mo.mo_num;

    let n_int = resolve_n_int(ezfio, mo_num)?;

    let mut out = io::stdout();
    writeln!(out, "Electrons: {alpha} {beta}").map_err(io_error)?;
    writeln!(out, "MO  : {mo_num}").map_err(io_error)?;
    writeln!(out, "n_int: {n_int}").map_err(io_error)?;

    // `work`: per-class accumulators starting at `"["`, one `,<i>` per
    // orbital (1-based), closed by `"]"` — the empty class is `"[]"`, a
    // non-empty one `"[,1,2]"` (the leading empty piece parses to
    // nothing).
    let mut acc: [String; 5] = [
        String::from("["),
        String::from("["),
        String::from("["),
        String::from("["),
        String::from("["),
    ];
    for (i, class) in mo.mo_class.iter().enumerate() {
        let kind = match class {
            MoClass::Core(_) => 0,
            MoClass::Inactive(_) => 1,
            MoClass::Active(_) => 2,
            MoClass::Virtual(_) => 3,
            MoClass::Deleted(_) => 4,
        };
        acc[kind].push_str(&format!(",{}", i + 1));
    }
    for s in acc.iter_mut() {
        s.push(']');
    }
    set_with(ezfio, &acc[0], &acc[1], &acc[2], &acc[3], &acc[4])
}

/// The `n_int` resolution of `get` (`qp_set_mo_class.ml:140-143`):
/// `try N_int_number.of_int (get_determinants_n_int()) with _ ->
/// Bitlist.n_int_of_mo_num mo_num`.
///
/// Measured subtleties reproduced:
/// * the `try` catches **everything** the direct path raises (missing
///   file, `int_of_string`, the `x >= 1` assert) and takes the fallback;
/// * the fallback's own `N_int_number.of_int` re-evaluates its `max`
///   default (`get_max`), which *re-reads* the file — a garbage value
///   there raises the **uncaught** `Failure("int_of_string")`;
/// * `n_int = 0` (or negative) takes the fallback and then fails the
///   maximum check against the re-read value (`N_int_number: <x>`);
/// * the `N_int > 128` warning is a raw `print_string` on both paths;
/// * the fallback forces `Qpackage.bit_kind_size` (see [`qpackage`]).
fn resolve_n_int(ezfio: &EzFio, mo_num: i64) -> Result<i64, AppError> {
    // The direct path.
    match read_n_int(ezfio) {
        Direct::Value(v) if v >= 1 => {
            // `N_int_number.of_int v` with `max = v`: the only possible
            // failure was the `x >= 1` assert (caught → fallback), so a
            // surviving value only warns.
            if v > 128 {
                print!("N_int > 128");
            }
            return Ok(v);
        }
        _ => {}
    }
    // The fallback: `Bitlist.n_int_of_mo_num mo_num` — uncaught from here.
    let root = qpackage::root()?;
    let bit_kind_size = qpackage::bit_kind_size(&root)?;
    // OCaml integer division truncates toward zero, like Rust's.
    let x = (mo_num - 1) / bit_kind_size + 1;
    // `N_int_number.of_int x` with the `max` default re-read:
    let max = match read_n_int(ezfio) {
        Direct::Absent => 128,
        Direct::Garbage(error) => return Err(error),
        Direct::Value(v) => v,
    };
    if x > 128 {
        print!("N_int > 128");
    }
    if max != 1 && x > max {
        return Err(AppError(OcamlException::InvalidArgument(format!(
            "N_int_number: {x}"
        ))));
    }
    Ok(x)
}

/// One read of `determinants/n_int`, in the reference's terms.
enum Direct {
    /// The `has` check failed: `get_max` returns 128.
    Absent,
    /// The file exists but does not parse: the direct call raises (caught
    /// by the `try`), the fallback's `max` re-read raises (uncaught).
    Garbage(AppError),
    /// A parsed value.
    Value(i64),
}

fn read_n_int(ezfio: &EzFio) -> Direct {
    if !ezfio.has("determinants", "n_int").unwrap_or(false) {
        return Direct::Absent;
    }
    match read_int_scalar(ezfio, "determinants", "n_int") {
        Ok(v) => Direct::Value(v),
        Err(e) => Direct::Garbage(e),
    }
}

/// `Ezfio.read_int` (`ezfio.ml:138-142,160`): `open_in`, first line,
/// `String.trim`, OCaml `int_of_string`. The failures are uncaught
/// wherever the reference does not wrap them:
/// * missing file → `Sys_error("<path>: No such file or directory")`;
/// * empty file → `End_of_file`;
/// * non-integer first line → `Failure("int_of_string")`.
fn read_int_scalar(ezfio: &EzFio, group: &str, key: &str) -> Result<i64, AppError> {
    let path = ezfio.root().join(group).join(key);
    let text = match fs::read(&path) {
        Ok(bytes) => String::from_utf8_lossy(&bytes).into_owned(),
        Err(error) => {
            let detail = match error.raw_os_error() {
                Some(code) => format!("{}: {}", path.display(), strerror(code)),
                None => error.to_string(),
            };
            return Err(AppError(OcamlException::SysError(detail)));
        }
    };
    match text.lines().next() {
        None => Err(AppError(OcamlException::EndOfFile)),
        Some(line) => ocaml_num::int_of_string(line.trim())
            .map_err(|_| AppError(OcamlException::Failure("int_of_string".into()))),
    }
}

/// The C library's `strerror` text for the errnos a scalar read can hit
/// (the `qp_create_ezfio` table).
fn strerror(code: i32) -> String {
    match code {
        1 => "Operation not permitted".to_owned(),
        2 => "No such file or directory".to_owned(),
        13 => "Permission denied".to_owned(),
        20 => "Not a directory".to_owned(),
        21 => "Is a directory".to_owned(),
        40 => "Too many levels of symbolic links".to_owned(),
        other => format!("os error {other}"),
    }
}

fn failure(message: impl Into<String>) -> AppError {
    AppError(OcamlException::Failure(message.into()))
}

fn io_error(error: io::Error) -> AppError {
    AppError(OcamlException::SysError(error.to_string()))
}

/// An `qp-ezfio` IO failure, rendered as the closest uncaught `Sys_error`.
///
/// The reference's write path wraps `gzip` in `open_process_out` and
/// raises `Failure("Failed writing compressed file …")` on a nonzero exit;
/// nothing in the test matrix reaches those paths, so the message text is
/// not pinned by a gate.
fn app(error: impl std::fmt::Display) -> AppError {
    AppError(OcamlException::SysError(error.to_string()))
}
