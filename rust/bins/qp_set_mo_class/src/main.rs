//! `qp_set_mo_class` — Rust port of `ocaml/qp_set_mo_class.ml` (P4, wave 3).
//!
//! The reference driver (`let ()` at the end of the file) has **no** `with`
//! handler: every exception escaped to the OCaml runtime, which printed
//! `Fatal error: exception <exception>` on stderr and exited 2. Wave S2
//! replaces that with one plain line, `qp_set_mo_class: <message>` on
//! standard error, exit 2. The flow is the reference's statement order:
//!
//! 1. the `QP_ROOT` guard, raised by a module initialiser in the
//!    reference (before the option parser runs — even a bare invocation,
//!    and even `-h`, dies here);
//! 2. the clap parse: the five range options `-c`/`-i`/`-a`/`-v`/`-d`,
//!    `-q`, and the mandatory single positional `EZFIO_DIR`; `-h`/`--help`
//!    print the help once and exit 0;
//! 3. `Ezfio.set_file` (create the tree when absent: `mkdir`, the `ezfio`
//!    group, the `creation`/`user`/`library` stamps);
//! 4. `has mo_basis/mo_num`, then [`set_with`] or [`get`] (the `-q` path:
//!    electrons, MOs, the `n_int` resolution, then `set` re-run on the
//!    current classes).
//!
//! `build_mask` (`qp_set_mo_class.ml:9-28`) is dead code in the reference
//! and is not ported.
//!
//! # stdout (functional, kept byte for byte)
//!
//! The five debug lines, the `Electrons:`/`MO  :`/`n_int:` summary, the
//! `check` wrapper's `Number of MOs: {n}` line, the raw `N_int > 128`
//! warning and the `More than 10_000 MOs` warning of the shared `Mo_basis`
//! reader. Only the `More than 10_000 MOs` warning of
//! [`range::create_class`] moved to standard error (wave S2).
//!
//! # Deviations
//!
//! A non-UTF-8 `argv` byte sequence is rejected by clap with its own
//! plain diagnostic (the reference walked raw bytes and rejected it with
//! its own message).

mod qpackage;
mod range;

use std::fs;
use std::io;
use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use clap::Parser;
use qp_ezfio::{Array, EzFio};
use qp_input::handlers;
use qp_input::handlers::mo_basis::MoClass;
use qp_types::error::ToolError;

/// The tool name in the `qp_set_mo_class: <message>` error line.
pub(crate) const TOOL: &str = "qp_set_mo_class";

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

/// The command line: the five class ranges, `-q` and `EZFIO_DIR`.
#[derive(Parser)]
struct Args {
    /// The core orbital range.
    #[arg(short = 'c', long = "core")]
    core: Option<String>,
    /// The inactive orbital range.
    #[arg(short = 'i', long = "inact")]
    inact: Option<String>,
    /// The active orbital range.
    #[arg(short = 'a', long = "act")]
    act: Option<String>,
    /// The virtual orbital range.
    #[arg(short = 'v', long = "virt")]
    virt: Option<String>,
    /// The deleted orbital range.
    #[arg(short = 'd', long = "del")]
    del: Option<String>,
    /// Query mode: read the current classes and re-write them.
    #[arg(short = 'q', long = "query")]
    query: bool,
    /// The EZFIO directory (created when absent).
    ezfio_dir: String,
}

fn main() -> ExitCode {
    // Module-initialiser semantics: the reference's `Qpackage.root` guard
    // runs before the parser, so even a bare invocation (and even `-h`)
    // dies here. Presence is what is checked: `QP_ROOT=` is accepted.
    if std::env::var_os("QP_ROOT").is_none() {
        let error = ToolError::new(TOOL, "QP_ROOT is not set (source quantum_package.rc)");
        eprintln!("{error}");
        return ExitCode::from(2);
    }

    let args = Args::parse();
    match run_action(&args.ezfio_dir, &args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::from(2)
        }
    }
}

/// `run` (`qp_set_mo_class.ml:176-185`): `set_file`, the `mo_num` guard,
/// then `get` or `set`.
fn run_action(dir: &str, args: &Args) -> Result<(), ToolError> {
    let ezfio = set_file(dir)?;
    if !ezfio
        .has("mo_basis", "mo_num")
        .map_err(|error| tool(&error))?
    {
        return Err(ToolError::new(TOOL, "mo_basis/mo_num not found"));
    }
    if args.query {
        get(&ezfio)
    } else {
        set_with(
            &ezfio,
            args.core.as_deref().unwrap_or("[]"),
            args.inact.as_deref().unwrap_or("[]"),
            args.act.as_deref().unwrap_or("[]"),
            args.virt.as_deref().unwrap_or("[]"),
            args.del.as_deref().unwrap_or("[]"),
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
///
/// The reference's `mkdir` is a single-level `Unix.mkdir`: a missing
/// parent is `ENOENT`, an existing path (file, or a directory without the
/// `.version` marker) is `EEXIST`. Wave S2: plain `io::Error` text.
fn set_file(dir: &str) -> Result<EzFio, ToolError> {
    let no_registry = || Box::new(NoRegistry) as Box<dyn qp_ezfio::Registry>;
    let root = Path::new(dir);
    if root.join(".version").exists() {
        return EzFio::open(root)
            .map_err(|error| ToolError::new(TOOL, format!("cannot open {dir}: {error}")))
            .map(|ezfio| ezfio.with_registry(no_registry()));
    }
    if let Some(parent) = root.parent()
        && !parent.as_os_str().is_empty()
        && !parent.exists()
    {
        // Raw errno 2 (`ENOENT`), the same value on Linux and macOS.
        let error = io::Error::from_raw_os_error(2);
        return Err(ToolError::new(
            TOOL,
            format!("cannot create {dir}: {error}"),
        ));
    }
    if root.exists() {
        // Raw errno 17 (`EEXIST`), the same value on Linux and macOS.
        let error = io::Error::from_raw_os_error(17);
        return Err(ToolError::new(
            TOOL,
            format!("cannot create {dir}: {error}"),
        ));
    }
    let ezfio = EzFio::create(root)
        .map_err(|error| ToolError::new(TOOL, format!("cannot create {dir}: {error}")))?
        .with_registry(no_registry());
    // The reference leaves a partial tree on a later failure (no cleanup).
    ezfio.create_group("ezfio").map_err(|error| {
        ToolError::new(
            TOOL,
            format!("cannot create the ezfio group in {dir}: {error}"),
        )
    })?;
    write_stamps(root)?;
    Ok(ezfio)
}

/// The `ezfio/{creation,user,library}` stamps of `Ezfio.set_file`
/// (`ezfio.ml:431-434`), the `qp_create_ezfio` pattern.
fn write_stamps(ezfio_file: &Path) -> Result<(), ToolError> {
    let dir = ezfio_file.join("ezfio");
    let failure = || {
        ToolError::new(
            TOOL,
            format!("cannot write the ezfio stamps in {}", ezfio_file.display()),
        )
    };
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
/// apply them (duplicate → error), require full coverage, print the five
/// debug lines, write `mo_basis/mo_class`.
fn set_with(
    ezfio: &EzFio,
    core: &str,
    inact: &str,
    act: &str,
    virt: &str,
    del: &str,
) -> Result<(), ToolError> {
    let mo_num = read_int_scalar(ezfio, "mo_basis", "mo_num")?;
    // `Array.init mo_num`: a negative count fails before any range is read.
    if mo_num < 0 {
        return Err(ToolError::new(
            TOOL,
            format!("mo_num must be non-negative: {mo_num}"),
        ));
    }
    let mo_num_us = mo_num as usize;

    // Parse in the reference's order; the first failure wins.
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
            // `mo_class.(i-1)`: array indexing is out of bounds.
            // Reachable only when `mo_num = 1` (the `max = 1`
            // maximum-check skip) — measured.
            if i < 1 || i as usize > mo_num_us {
                return Err(ToolError::new(
                    TOOL,
                    format!("orbital index out of bounds: {i}"),
                ));
            }
            let idx = i as usize - 1;
            match assigned[idx] {
                None => assigned[idx] = Some(kind),
                Some(prev) => {
                    return Err(ToolError::new(
                        TOOL,
                        format!(
                            "Orbital {i} is defined both in the {} and {} spaces",
                            KEYWORDS[prev as usize], KEYWORDS[kind as usize]
                        ),
                    ));
                }
            }
        }
    }

    for (idx, slot) in assigned.iter().enumerate() {
        if slot.is_none() {
            return Err(ToolError::new(
                TOOL,
                format!("Orbital {} is not specified (mo_num = {mo_num})", idx + 1),
            ));
        }
    }

    // Debug output, in the reference's order, full `MO_class.to_string`.
    let mut out = io::stdout().lock();
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
        writeln!(out, "{line}").map_err(|error| tool(&error))?;
    }
    drop(out);

    let data: Vec<String> = assigned
        .iter()
        .map(|slot| match slot {
            Some(kind) => KEYWORDS[*kind as usize].to_owned(),
            None => unreachable!("full coverage checked above"),
        })
        .collect();
    let array = Array::new(vec![mo_num_us], data).map_err(|error| {
        ToolError::new(TOOL, format!("cannot build the mo_class array: {error}"))
    })?;
    ezfio
        .set_array::<String>("mo_basis", "mo_class", &array)
        .map_err(|error| ToolError::new(TOOL, format!("cannot write mo_basis/mo_class: {error}")))
}

/// `get` (`qp_set_mo_class.ml:116-172`, the `-q` path): read the electrons
/// and the MOs, resolve `n_int`, print the three summary lines, then call
/// `set` on the current classes (which re-parses, re-prints and re-writes
/// `mo_basis/mo_class`).
fn get(ezfio: &EzFio) -> Result<(), ToolError> {
    // `Input.Electrons.read`: `None` only when `elec_alpha_num` is absent;
    // a read failure of the *present* alpha propagates, as does the beta
    // read (which is not `has`-guarded).
    if !ezfio
        .has("electrons", "elec_alpha_num")
        .map_err(|error| tool(&error))?
    {
        return Err(ToolError::new(TOOL, "Unable to read electrons"));
    }
    let alpha = read_int_scalar(ezfio, "electrons", "elec_alpha_num")?;
    let beta = read_int_scalar(ezfio, "electrons", "elec_beta_num")?;
    // `Elec_alpha_number.of_int` / `Elec_beta_number.of_int` asserts.
    if alpha <= 0 {
        return Err(ToolError::new(
            TOOL,
            format!("alpha electron count {alpha} must be positive"),
        ));
    }
    if beta < 0 {
        return Err(ToolError::new(
            TOOL,
            format!("beta electron count {beta} must be non-negative"),
        ));
    }

    // `Input.Mo_basis.read`: classify the `mo_num` read ourselves, then let
    // the shared handler do the reference's field reads — the `mo_label`
    // and `mo_class` materialisation, the `More alpha electrons than MOs`
    // guard, the `MO_number` conversion (its `More than 10_000 MOs`
    // warning included, still on stdout), the `ao_md5` check.
    if !ezfio
        .has("mo_basis", "mo_num")
        .map_err(|error| tool(&error))?
    {
        return Err(ToolError::new(TOOL, "Unable to read MOs"));
    }
    read_int_scalar(ezfio, "mo_basis", "mo_num")?;
    let mo = match handlers::mo_basis::read(ezfio) {
        // Unreachable after the `has` guard above; the reference's
        // `failwith "Unable to read MOs"` covers it.
        Ok(None) => return Err(ToolError::new(TOOL, "Unable to read MOs")),
        Ok(Some(mo)) => mo,
        Err(error) => return Err(tool(&error)),
    };
    let mo_num = mo.mo_num;

    let n_int = resolve_n_int(ezfio, mo_num)?;

    let mut out = io::stdout().lock();
    writeln!(out, "Electrons: {alpha} {beta}").map_err(|error| tool(&error))?;
    writeln!(out, "MO  : {mo_num}").map_err(|error| tool(&error))?;
    writeln!(out, "n_int: {n_int}").map_err(|error| tool(&error))?;
    drop(out);

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
///   file, a non-integer value, the `x >= 1` check) and takes the
///   fallback;
/// * the fallback's own `N_int_number` re-evaluates its `max`
///   default (`get_max`), which *re-reads* the file — a garbage value
///   there is the uncaught failure;
/// * `n_int = 0` (or negative) takes the fallback and then fails the
///   maximum check against the re-read value;
/// * the `N_int > 128` warning is a raw `print_string` on both paths
///   (stdout, functional);
/// * the fallback forces `Qpackage.bit_kind_size` (see [`qpackage`]).
fn resolve_n_int(ezfio: &EzFio, mo_num: i64) -> Result<i64, ToolError> {
    // The direct path.
    match read_n_int(ezfio) {
        Direct::Value(v) if v >= 1 => {
            // `N_int_number.of_int v` with `max = v`: the only possible
            // failure was the `x >= 1` check (caught → fallback), so a
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
        return Err(ToolError::new(
            TOOL,
            format!("n_int {x} is above the maximum of {max}"),
        ));
    }
    Ok(x)
}

/// One read of `determinants/n_int`, in the reference's terms.
enum Direct {
    /// The `has` check failed: `get_max` returns 128.
    Absent,
    /// The file exists but does not parse: the direct call raises (caught
    /// by the `try`), the fallback's `max` re-read raises (uncaught).
    Garbage(ToolError),
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
/// `String.trim`, `int_of_string`. Wave S2: plain messages, decimal-only
/// `i64::from_str`.
fn read_int_scalar(ezfio: &EzFio, group: &str, key: &str) -> Result<i64, ToolError> {
    let path = ezfio.root().join(group).join(key);
    let bytes = fs::read(&path).map_err(|error| {
        ToolError::new(TOOL, format!("cannot read {}: {error}", path.display()))
    })?;
    let text = String::from_utf8_lossy(&bytes).into_owned();
    let line = text.lines().next().unwrap_or("").trim();
    if line.is_empty() {
        return Err(ToolError::new(
            TOOL,
            format!("cannot read {}: the file is empty", path.display()),
        ));
    }
    line.parse::<i64>()
        .map_err(|_| ToolError::new(TOOL, format!("invalid integer '{line}' in {group}/{key}")))
}

/// A failure that carries a displayable reason (0MQ/EZFIO/anyhow), kept
/// in its own text.
fn tool(error: &impl std::fmt::Display) -> ToolError {
    ToolError::new(TOOL, error.to_string())
}
