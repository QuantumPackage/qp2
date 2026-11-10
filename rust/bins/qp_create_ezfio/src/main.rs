//! `qp_create_ezfio` — Rust port of `ocaml/qp_create_ezfio.ml` (P4, wave 3;
//! wave S2 text release).
//!
//! This crate reproduces the reference end to end: the `QP_ROOT` guard, the
//! clap command line (wave S2 replaces the hand-rolled `Getopt` port), the
//! geometry reader (see [`geometry`]), the basis/pseudopotential input
//! (see [`basis_read`]), and the EZFIO write (see [`write`], whose value
//! computations live in [`specs`]).
//!
//! The reference checked `QP_ROOT` from a module initialiser, i.e.
//! **before** the option parser ran: even `-h` died there. An empty value is
//! accepted; the check is presence only.
//!
//! # Flow
//!
//! 1. the `QP_ROOT` guard, before any parsing;
//! 2. the clap parse: `-o`/`-b`/`-c`/`-g`/`-m`/`-p`/`-q`/`-x`/`-a` plus the
//!    required single positional `FILE`; `-h`/`--help` print the help once
//!    and exit 0 (the reference printed a captured block and continued);
//! 3. the manual checks, in the reference's order (first failure wins):
//!    `--basis` mandatory, then the `charge` → `ghost` → `multiplicity`
//!    conversions (decimal-only `FromStr`, wave S2.3; the `x >= 1`
//!    multiplicity value check stays in [`geometry::read`]);
//! 4. the action: geometry (including ghosts), the basis table, the
//!    pseudopotential table, the `-o` "already exists" check, the
//!    per-nucleus basis reads, and only then the write.
//!
//! Every failure prints exactly one line, `qp_create_ezfio: <message>`, on
//! standard error and exits 2 (both classes of the reference's driver — the
//! caught `Failure` and the escaped exceptions — exited 2). Success prints
//! the directory name on stdout and exits 0.
//!
//! # Units
//!
//! [`Options`] carries the raw strings and numbers; [`geometry`] attaches
//! bohr/angstrom and interprets `ghost` as a multiplier of the covalent
//! radii; [`basis_read`] converts nothing.
//!
//! # Deviations
//!
//! * A non-UTF-8 `argv` byte sequence is rejected by clap with its own
//!   plain diagnostic (the reference walked raw bytes and rejected it
//!   with its own message).
//! * A missing `FILE` now fails at the clap parse, before the numeric
//!   conversions; the reference checked it after them.
//! * Long options accept `--opt value` as well as `--opt=value`, and short
//!   options follow the clap grammar (clusters, attached values); the
//!   reference accepted only the `=` form for long options.

mod basis_read;
mod geometry;
mod specs;
mod write;

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use qp_types::error::ToolError;

/// The tool name in the `qp_create_ezfio: <message>` error line.
pub(crate) const TOOL: &str = "qp_create_ezfio";

/// The command line: the reference's option table plus the mandatory
/// `FILE` positional.
#[derive(Parser, Debug)]
struct Args {
    /// The EZFIO directory to create (default: the input stem + `.ezfio`).
    #[arg(short = 'o', long = "output")]
    output: Option<String>,
    /// The basis-set specification (mandatory).
    #[arg(short = 'b', long = "basis")]
    basis: Option<String>,
    /// The input geometry is already in bohr.
    #[arg(short = 'a', long = "au")]
    au: bool,
    /// The total molecular charge; a leading `m` negates the rest (`m1` =
    /// -1).
    #[arg(short = 'c', long = "charge")]
    charge: Option<String>,
    /// Ghost-atom threshold, in units of the covalent radii (`0.0`
    /// disables).
    #[arg(short = 'g', long = "ghost")]
    ghost: Option<String>,
    /// The spin multiplicity `2S+1` (default 1).
    #[arg(short = 'm', long = "multiplicity")]
    multiplicity: Option<String>,
    /// The pseudopotential specification.
    #[arg(short = 'p', long = "pseudo")]
    pseudo: Option<String>,
    /// Use the 20210101 EZFIO convention.
    #[arg(short = 'q', long = "qmc")]
    qmc: bool,
    /// Request Cartesian AOs.
    #[arg(short = 'x', long = "cartesian")]
    cartesian: bool,
    /// The input XYZ or Z-matrix file.
    file: String,
}

/// The parsed and defaulted invocation.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    /// EZFIO directory to create. `-o` value, or the input-derived default.
    pub output: PathBuf,
    /// Basis-set specification (mandatory unless `-h` was given).
    pub basis: String,
    /// Input geometry is already in bohr (`-a`).
    pub au: bool,
    /// Total charge. `m1` on the command line means `-1`.
    pub charge: i64,
    /// Ghost-atom threshold, in units of the covalent radii; `0.0`
    /// disables.
    pub ghost: f64,
    /// Spin multiplicity `2S+1` (`-m`), default 1.
    pub multiplicity: i64,
    /// Pseudopotential specification (`-p`), if any.
    pub pseudo: Option<String>,
    /// Use the 20210101 EZFIO convention (`-q`).
    pub qmc: bool,
    /// Request Cartesian AOs (`-x`).
    pub cartesian: bool,
    /// Input XYZ or Z-matrix file (the `FILE` positional).
    pub input: PathBuf,
}

fn main() -> ExitCode {
    // Module-initialiser semantics: the reference's `Qpackage.root` guard
    // runs before the parser, so even `-h` dies here. Presence is what is
    // checked: `QP_ROOT=` is accepted.
    let Some(root) = std::env::var_os("QP_ROOT") else {
        let error = ToolError::new(TOOL, "QP_ROOT is not set (source quantum_package.rc)");
        eprintln!("{error}");
        return ExitCode::from(2);
    };

    let args = Args::parse();
    let options = match validate(&args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            return ExitCode::from(2);
        }
    };
    match create_ezfio(&options, Path::new(&root)) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::from(2)
        }
    }
}

/// The manual checks, in the reference's order (first failure wins):
/// `--basis` mandatory, then the `charge` → `ghost` → `multiplicity`
/// conversions. The `FILE` positional is required by clap itself, so it
/// fails at the parse.
fn validate(args: &Args) -> Result<Options, ToolError> {
    let basis = args
        .basis
        .as_deref()
        .ok_or_else(|| ToolError::new(TOOL, "--basis option is missing"))?;
    let charge = match args.charge.as_deref() {
        Some(text) => convert_charge(text)?,
        None => 0,
    };
    let ghost = match args.ghost.as_deref() {
        Some(text) => text.parse::<f64>().map_err(|_| {
            ToolError::new(TOOL, format!("invalid float '{text}' in ghost threshold"))
        })?,
        None => 0.0,
    };
    let multiplicity = match args.multiplicity.as_deref() {
        Some(text) => text.parse::<i64>().map_err(|_| {
            ToolError::new(TOOL, format!("invalid integer '{text}' in multiplicity"))
        })?,
        None => 1,
    };
    let output = args
        .output
        .as_ref()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(default_output(&args.file)));
    Ok(Options {
        output,
        basis: basis.to_owned(),
        au: args.au,
        charge,
        ghost,
        multiplicity,
        pseudo: args.pseudo.clone(),
        qmc: args.qmc,
        cartesian: args.cartesian,
        input: PathBuf::from(&args.file),
    })
}

/// The reference's charge parser: a leading `m` negates the rest, allowing
/// a negative charge without a token the parser would mistake for an
/// option. Wave S2.3: the inner integer is decimal-only
/// (`i64::from_str`), so `m+1` is -1, `m-1` is +1, and `0x10` is rejected.
fn convert_charge(text: &str) -> Result<i64, ToolError> {
    let value = match text.strip_prefix('m') {
        Some(rest) => parse_int(rest, "charge")?,
        None => parse_int(text, "charge")?,
    };
    Ok(if text.starts_with('m') { -value } else { value })
}

/// Decimal-only integer parse with the plain `invalid integer` message.
fn parse_int(text: &str, what: &str) -> Result<i64, ToolError> {
    text.parse::<i64>()
        .map_err(|_| ToolError::new(TOOL, format!("invalid integer '{text}' in {what}")))
}

/// The default EZFIO directory name: `stem.ezfio` when the input ends in
/// `.xyz` or `.zmt` (split at the last point), else the whole input name
/// plus `.ezfio`.
fn default_output(input: &str) -> String {
    match input.rsplit_once('.') {
        Some((stem, "xyz" | "zmt")) => format!("{stem}.ezfio"),
        _ => format!("{input}.ezfio"),
    }
}

/// Perform the reference's action: build the EZFIO directory from `options`.
///
/// The input phases run in the reference's order (geometry including
/// ghosts, then the basis table, then the pseudopotential table), then the
/// `-o` "already exists" check, then the per-nucleus basis reads, and only
/// then the write. The check sits between the table builds and the reads
/// (`qp_create_ezfio.ml:273` precedes the per-nucleus reads inside
/// `write_file`), so an existing output directory wins over a basis-read
/// failure, exactly as measured.
///
/// # Errors
///
/// Whatever the first failing phase returns.
fn create_ezfio(options: &Options, root: &Path) -> Result<(), ToolError> {
    let geometry = geometry::read(
        &options.input,
        options.au,
        options.charge,
        options.multiplicity,
        options.ghost,
    )?;
    let mut basis_table = basis_read::build_basis_table(&options.basis, &geometry.nuclei, root)?;
    let pseudo_table = basis_read::read_pseudo(options.pseudo.as_deref(), &geometry.nuclei, root)?;
    write::check_output_available(&options.output)?;
    let bases = basis_read::read_basis(&mut basis_table, &geometry.nuclei)?;
    write::write(options, &geometry, &bases, &pseudo_table)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_output_names() {
        for (input, output) in [
            ("h2.xyz", "h2.ezfio"),
            ("h2.zmt", "h2.ezfio"),
            ("h2.dat", "h2.dat.ezfio"),
            ("h2", "h2.ezfio"),
            ("noext.tar.xyz", "noext.tar.ezfio"),
            (".xyz", ".ezfio"),
            ("a.b.xyz", "a.b.ezfio"),
            ("h2.", "h2..ezfio"),
        ] {
            assert_eq!(default_output(input), output, "input {input}");
        }
    }

    #[test]
    fn the_charge_m_prefix_negates() {
        assert_eq!(convert_charge("5").expect("plain"), 5);
        assert_eq!(convert_charge("+5").expect("signed"), 5);
        assert_eq!(convert_charge("-5").expect("negative"), -5);
        assert_eq!(convert_charge("m1").expect("m1"), -1);
        assert_eq!(convert_charge("m0").expect("m0"), 0);
        // The `m` prefix negates the parsed integer, including a signed
        // one: `m+1` is -1 and `m-1` is +1.
        assert_eq!(convert_charge("m+1").expect("m+1"), -1);
        assert_eq!(convert_charge("m-1").expect("m-1"), 1);
    }

    #[test]
    fn the_charge_is_decimal_only() {
        // OCaml `int_of_string` accepted hex and underscores; the port is
        // decimal-only (wave S2.3, recorded deviation).
        let error = convert_charge("0x10").expect_err("hex");
        assert_eq!(
            error.to_string(),
            "qp_create_ezfio: invalid integer '0x10' in charge"
        );
        let error = convert_charge("1_000").expect_err("underscore");
        assert_eq!(
            error.to_string(),
            "qp_create_ezfio: invalid integer '1_000' in charge"
        );
        let error = convert_charge("m").expect_err("bare m");
        assert_eq!(
            error.to_string(),
            "qp_create_ezfio: invalid integer '' in charge"
        );
    }

    #[test]
    fn validate_checks_basis_charge_ghost_multiplicity_in_order() {
        let base = || Args {
            output: None,
            basis: Some("x".to_owned()),
            au: false,
            charge: None,
            ghost: None,
            multiplicity: None,
            pseudo: None,
            qmc: false,
            cartesian: false,
            file: "f.xyz".to_owned(),
        };
        // No basis: the mandatory check wins over any conversion.
        let mut args = base();
        args.basis = None;
        args.charge = Some("abc".to_owned());
        let error = validate(&args).expect_err("basis first");
        assert_eq!(
            error.to_string(),
            "qp_create_ezfio: --basis option is missing"
        );
        // Then the conversions, charge before ghost before multiplicity.
        let mut args = base();
        args.charge = Some("abc".to_owned());
        args.ghost = Some("xx".to_owned());
        let error = validate(&args).expect_err("charge before ghost");
        assert_eq!(
            error.to_string(),
            "qp_create_ezfio: invalid integer 'abc' in charge"
        );
        let mut args = base();
        args.ghost = Some("xx".to_owned());
        args.multiplicity = Some("yy".to_owned());
        let error = validate(&args).expect_err("ghost before multiplicity");
        assert_eq!(
            error.to_string(),
            "qp_create_ezfio: invalid float 'xx' in ghost threshold"
        );
        let mut args = base();
        args.multiplicity = Some("yy".to_owned());
        let error = validate(&args).expect_err("multiplicity");
        assert_eq!(
            error.to_string(),
            "qp_create_ezfio: invalid integer 'yy' in multiplicity"
        );
    }
}
