//! The reference command line, parsed by a literal port of the OCaml
//! `Getopt` library plus the driver's own post-processing.
//!
//! `clap` cannot express this surface: option arguments may be attached
//! (`-bcc-pvdz`), long options only take arguments with `=` (`--basis=a`,
//! never `--basis a`), short options can be concatenated (`-ax`), a lone
//! `-` and an empty token are anonymous, `--` ends parsing, and `-h` prints
//! the help and **continues** instead of exiting. The error messages are
//! part of the observable contract (exit 2), so the scan below follows
//! `getopt.ml` statement by statement:
//!
//! * a long option without `=` never consumes the next token;
//! * a short option that needs an argument takes the rest of its own token
//!   verbatim, else the next token if the next token does not start with
//!   `-`; an empty next token raises `Invalid_argument("index out of
//!   bounds")` (the reference indexes `args.(no+1).[0]`);
//! * every handler runs while scanning, so a later syntax error suppresses
//!   an earlier `-h`;
//! * after a successful scan the driver prints the help once if `-h` was
//!   seen, and only then checks the mandatory `--basis`; with `-h` the
//!   mandatory check is skipped entirely;
//! * numeric conversions run charge → ghost → multiplicity, before the
//!   anonymous `FILE` check.
//!
//! # Deviations
//!
//! The help block is the reference capture for `argv[0] = "qp_create_ezfio"`;
//! the reference reflows it to width 79 for a longer name, while this port
//! only substitutes the name into the captured lines. A non-ASCII short
//! option is reported as `Unknown option -<char>` (valid UTF-8); the
//! reference works on bytes and would print the first raw byte. No option in
//! the table is non-ASCII, so the message is only reachable for typos.
//!
//! # Units
//!
//! `ghost` is the multiplier of the covalent radii (dimensionless). The
//! other fields are strings and flags; geometry handling attaches units.

use std::io::Write;
use std::path::PathBuf;

use crate::error::{AppError, OcamlException};
use crate::ocaml_num::{float_of_string, int_of_string};

/// Name used in diagnostics and as the fallback `argv[0]`.
pub const PROGRAM_NAME: &str = "qp_create_ezfio";

/// The bare-name help capture from the reference build.
///
/// 1817 bytes, 45 newlines. `qp_create_ezfio` occurs exactly twice — the
/// header and the usage line — so substituting with [`String::replace`]
/// reproduces the reference's two substitution points. The reference
/// reflows the block to width 79 for a longer `argv[0]`; this port
/// reproduces the bare-name form only (documented deviation).
const HELP_TEMPLATE: &str = include_str!("../tests/data/help_bare.txt");

/// Render the help block with `program` substituted for the invoked name.
#[must_use]
pub fn render_help(program: &str) -> String {
    HELP_TEMPLATE.replace(PROGRAM_NAME, program)
}

/// One option in the reference's `set_specs` table (after the `-h` entry
/// that `set_specs` prepends).
struct Spec {
    /// Single-character short name.
    short: char,
    /// Long name, without the leading `--`.
    long: &'static str,
    /// Whether the option takes an argument.
    takes_arg: bool,
}

/// The option table, in `set_specs` order.
const SPECS: [Spec; 10] = [
    Spec {
        short: 'h',
        long: "help",
        takes_arg: false,
    },
    Spec {
        short: 'o',
        long: "output",
        takes_arg: true,
    },
    Spec {
        short: 'b',
        long: "basis",
        takes_arg: true,
    },
    Spec {
        short: 'a',
        long: "au",
        takes_arg: false,
    },
    Spec {
        short: 'c',
        long: "charge",
        takes_arg: true,
    },
    Spec {
        short: 'g',
        long: "ghost",
        takes_arg: true,
    },
    Spec {
        short: 'm',
        long: "multiplicity",
        takes_arg: true,
    },
    Spec {
        short: 'p',
        long: "pseudo",
        takes_arg: true,
    },
    Spec {
        short: 'q',
        long: "qmc",
        takes_arg: false,
    },
    Spec {
        short: 'x',
        long: "cartesian",
        takes_arg: false,
    },
];

/// Find a spec by its long name.
fn find_long(name: &str) -> Result<&'static Spec, AppError> {
    SPECS
        .iter()
        .find(|spec| spec.long == name)
        .ok_or_else(|| getopt_error(format!("Unknown option --{name}")))
}

/// Find a spec by its short name.
fn find_short(name: char) -> Result<&'static Spec, AppError> {
    SPECS
        .iter()
        .find(|spec| spec.short == name)
        .ok_or_else(|| getopt_error(format!("Unknown option -{name}")))
}

/// `Getopt.Error`, which the driver catches and prefixes.
fn getopt_error(message: String) -> AppError {
    AppError::uncaught(OcamlException::Getopt(message))
}

/// `Invalid_argument("index out of bounds")`, raised by the reference when
/// an option argument is looked up in an empty string.
fn index_out_of_bounds() -> AppError {
    AppError::uncaught(OcamlException::InvalidArgument(
        "index out of bounds".into(),
    ))
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
    /// Ghost-atom threshold, in units of the covalent radii; `0.0` disables.
    pub ghost: f64,
    /// Spin multiplicity `2S+1` (`-m`), default 1.
    pub multiplicity: i64,
    /// Pseudopotential specification (`-p`), if any.
    pub pseudo: Option<String>,
    /// Use the 20210101 EZFIO convention (`-q`).
    pub qmc: bool,
    /// Request Cartesian AOs (`-x`).
    pub cartesian: bool,
    /// Input XYZ or Z-matrix file (the first anonymous argument).
    pub input: PathBuf,
}

/// The option values as the `Hashtbl` of the reference holds them, before
/// defaulting and numeric conversion.
#[derive(Default)]
struct Raw {
    output: Option<String>,
    basis: Option<String>,
    au: bool,
    charge: Option<String>,
    ghost: Option<String>,
    multiplicity: Option<String>,
    pseudo: Option<String>,
    qmc: bool,
    cartesian: bool,
    help: bool,
}

impl Raw {
    /// Run the handler for one option, mirroring `Hashtbl.replace`.
    fn apply(&mut self, spec: &Spec, argument: Option<&str>) {
        match spec.short {
            'o' => self.output = argument.map(str::to_owned),
            'b' => self.basis = argument.map(str::to_owned),
            'a' => self.au = true,
            'c' => self.charge = argument.map(str::to_owned),
            'g' => self.ghost = argument.map(str::to_owned),
            'm' => self.multiplicity = argument.map(str::to_owned),
            'p' => self.pseudo = argument.map(str::to_owned),
            'q' => self.qmc = true,
            'x' => self.cartesian = true,
            'h' => self.help = true,
            other => unreachable!("internal invariant: SPECS covers option {other}"),
        }
    }
}

/// Parse `argv` (including `argv[0]`) into [`Options`].
///
/// Help text is written to `help_out` at the point the reference prints it:
/// after a complete scan, before the mandatory-basis check and the numeric
/// conversions. A scan error therefore produces no help, while a conversion
/// or missing-file error can follow help on stdout.
///
/// # Errors
///
/// Any [`AppError`] the reference would print and exit 2 for.
pub fn parse(argv: &[String], help_out: &mut dyn Write) -> Result<Options, AppError> {
    let mut raw = Raw::default();
    let mut anonymous: Vec<String> = Vec::new();

    let mut index = 1;
    while index < argv.len() {
        let token = &argv[index];
        let bytes = token.as_bytes();
        if bytes.is_empty() {
            anonymous.push(token.clone());
            index += 1;
            continue;
        }
        if bytes[0] != b'-' {
            anonymous.push(token.clone());
            index += 1;
            continue;
        }
        if bytes.len() >= 2 && bytes[1] == b'-' {
            if bytes.len() == 2 {
                // `--`: everything that follows is anonymous.
                anonymous.extend(argv[index + 1..].iter().cloned());
                break;
            }
            match token.find('=') {
                Some(equals) => {
                    let name = &token[2..equals];
                    let argument = &token[equals + 1..];
                    let spec = find_long(name)?;
                    if !spec.takes_arg {
                        return Err(getopt_error(format!(
                            "Option --{name} does not accept argument"
                        )));
                    }
                    raw.apply(spec, Some(argument));
                }
                None => {
                    let name = &token[2..];
                    let spec = find_long(name)?;
                    if spec.takes_arg {
                        // The reference raises before looking at the next
                        // token: long options take arguments only with `=`.
                        return Err(getopt_error(format!(
                            "Option {token} must have an argument"
                        )));
                    }
                    raw.apply(spec, None);
                }
            }
            index += 1;
            continue;
        }
        if bytes.len() == 1 {
            // A lone `-` is anonymous.
            anonymous.push(token.clone());
            index += 1;
            continue;
        }

        // Short option(s). The first byte after `-` is ASCII for every
        // entry of the table; an unknown non-ASCII char is reported and
        // returns before any slicing.
        let first = token[1..]
            .chars()
            .next()
            .expect("internal invariant: token has at least one byte after `-`");
        let spec = find_short(first)?;
        if !spec.takes_arg {
            raw.apply(spec, None);
            for extra in token[2..].chars() {
                let spec = find_short(extra)?;
                if spec.takes_arg {
                    return Err(getopt_error(format!(
                        "Only non-argument short-options can be concatenated (error with option {extra} in {token})"
                    )));
                }
                raw.apply(spec, None);
            }
            index += 1;
            continue;
        }

        // An argument-taking short option.
        if token.len() > 2 {
            raw.apply(spec, Some(&token[2..]));
            index += 1;
            continue;
        }
        match argv.get(index + 1) {
            Some(next) => {
                if next.is_empty() {
                    // `args.(no+1).[0]` on an empty string.
                    return Err(index_out_of_bounds());
                }
                if next.starts_with('-') {
                    return Err(getopt_error(format!(
                        "Option -{first} must have an argument"
                    )));
                }
                raw.apply(spec, Some(next));
                index += 2;
            }
            None => {
                return Err(getopt_error(format!(
                    "Option -{first} must have an argument"
                )));
            }
        }
    }

    if raw.help {
        let program = argv.first().map_or(PROGRAM_NAME, String::as_str);
        help_out
            .write_all(render_help(program).as_bytes())
            .map_err(|error| AppError::uncaught(OcamlException::SysError(error.to_string())))?;
    } else if raw.basis.is_none() {
        return Err(AppError::fail("--basis option is missing."));
    }

    let charge = convert_charge(raw.charge.as_deref())?;
    let ghost = match raw.ghost.as_deref() {
        None => 0.0,
        Some(text) => float_of_string(text).map_err(number_error)?,
    };
    let multiplicity = match raw.multiplicity.as_deref() {
        None => 1,
        Some(text) => int_of_string(text).map_err(number_error)?,
    };

    let input = anonymous
        .into_iter()
        .next()
        .ok_or_else(|| AppError::fail("input file is missing"))?;
    let output = match raw.output {
        Some(directory) => PathBuf::from(directory),
        None => PathBuf::from(default_output(&input)),
    };

    Ok(Options {
        output,
        basis: raw.basis.unwrap_or_default(),
        au: raw.au,
        charge,
        ghost,
        multiplicity,
        pseudo: raw.pseudo,
        qmc: raw.qmc,
        cartesian: raw.cartesian,
        input: PathBuf::from(input),
    })
}

/// `Failure("int_of_string")` / `Failure("float_of_string")`, caught and
/// prefixed by the driver.
fn number_error(error: crate::ocaml_num::NumParseError) -> AppError {
    AppError::fail(error.to_string())
}

/// The reference's charge parser: a leading `m` negates the rest, allowing
/// negative charges without a token that `getopt` would mistake for an
/// option.
fn convert_charge(text: Option<&str>) -> Result<i64, AppError> {
    let Some(text) = text else {
        return Ok(0);
    };
    match text.as_bytes().first() {
        // `x.[0]` on an empty string in the reference.
        None => Err(index_out_of_bounds()),
        Some(b'm') => int_of_string(&text[1..])
            .map(|value| -value)
            .map_err(number_error),
        Some(_) => int_of_string(text).map_err(number_error),
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_ok(args: &[&str]) -> Options {
        let argv: Vec<String> = args.iter().map(|argument| (*argument).to_owned()).collect();
        let mut help = Vec::new();
        parse(&argv, &mut help).expect("parse should succeed")
    }

    fn parse_err(args: &[&str]) -> (AppError, Vec<u8>) {
        let argv: Vec<String> = args.iter().map(|argument| (*argument).to_owned()).collect();
        let mut help = Vec::new();
        let error = parse(&argv, &mut help).expect_err("parse should fail");
        (error, help)
    }

    #[test]
    fn help_render_is_byte_identical_for_the_bare_name() {
        assert_eq!(render_help(PROGRAM_NAME), HELP_TEMPLATE);
        assert_eq!(render_help(PROGRAM_NAME).len(), 1817);
        assert_eq!(render_help(PROGRAM_NAME).matches('\n').count(), 45);
        let substituted = render_help("/tmp/bin/qp_create_ezfio");
        assert!(substituted.starts_with("/tmp/bin/qp_create_ezfio - Quantum Package command "));
        assert!(substituted.contains("\n  /tmp/bin/qp_create_ezfio [-a] -b <string>"));
    }

    #[test]
    fn defaults_and_output_naming() {
        let options = parse_ok(&["qp_create_ezfio", "-b", "cc-pvdz", "h2.xyz"]);
        assert_eq!(options.output, PathBuf::from("h2.ezfio"));
        assert_eq!(options.basis, "cc-pvdz");
        assert_eq!(options.charge, 0);
        assert_eq!(options.ghost, 0.0);
        assert_eq!(options.multiplicity, 1);
        assert_eq!(options.pseudo, None);
        assert!(!options.au && !options.qmc && !options.cartesian);
        assert_eq!(options.input, PathBuf::from("h2.xyz"));

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
            assert_eq!(
                parse_ok(&["qp_create_ezfio", "-b", "x", input]).output,
                PathBuf::from(output),
                "input {input}"
            );
        }
    }

    #[test]
    fn short_and_long_argument_forms() {
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-bcc-pvdz", "h2.xyz"]).basis,
            "cc-pvdz"
        );
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "--basis=cc-pvdz", "h2.xyz"]).basis,
            "cc-pvdz"
        );
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "--basis=", "h2.xyz"]).basis,
            ""
        );
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", "-o", "dir", "h2.xyz"]).output,
            PathBuf::from("dir")
        );
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", "-odir", "h2.xyz"]).output,
            PathBuf::from("dir")
        );
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", "--output=dir", "h2.xyz"]).output,
            PathBuf::from("dir")
        );
        // `Hashtbl.replace`: the last occurrence wins.
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", "-o", "a", "-o", "b", "h2.xyz"]).output,
            PathBuf::from("b")
        );
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", "-b", "y", "h2.xyz"]).basis,
            "y"
        );
    }

    #[test]
    fn flags_and_concatenation() {
        let options = parse_ok(&["qp_create_ezfio", "-aqx", "-b", "x", "h2.xyz"]);
        assert!(options.au && options.qmc && options.cartesian);

        let options = parse_ok(&["qp_create_ezfio", "-axq", "-m2", "-b", "x", "h2.xyz"]);
        assert!(options.au && options.cartesian && options.qmc);
        assert_eq!(options.multiplicity, 2);

        // An argument-taking option cannot appear inside a concatenation,
        // even when it has its own attached argument.
        let (error, help) = parse_err(&["qp_create_ezfio", "-axm2"]);
        assert_eq!(
            error.render(),
            "Fatal error: exception Getopt.Error(\"Only non-argument short-options can be concatenated (error with option m in -axm2)\")\n"
        );
        assert!(help.is_empty());

        let (error, help) = parse_err(&["qp_create_ezfio", "-ab"]);
        assert_eq!(
            error.render(),
            "Fatal error: exception Getopt.Error(\"Only non-argument short-options can be concatenated (error with option b in -ab)\")\n"
        );
        assert!(help.is_empty());

        let (error, _) = parse_err(&["qp_create_ezfio", "-aZ"]);
        assert_eq!(
            error.render(),
            "Fatal error: exception Getopt.Error(\"Unknown option -Z\")\n"
        );
        let (error, _) = parse_err(&["qp_create_ezfio", "-a=x"]);
        assert_eq!(
            error.render(),
            "Fatal error: exception Getopt.Error(\"Unknown option -=\")\n"
        );
    }

    #[test]
    fn getopt_error_messages() {
        for (args, expected) in [
            (
                vec!["qp_create_ezfio", "-Z"],
                "Fatal error: exception Getopt.Error(\"Unknown option -Z\")\n",
            ),
            (
                vec!["qp_create_ezfio", "--zzz"],
                "Fatal error: exception Getopt.Error(\"Unknown option --zzz\")\n",
            ),
            (
                vec!["qp_create_ezfio", "--zzz=x"],
                "Fatal error: exception Getopt.Error(\"Unknown option --zzz\")\n",
            ),
            (
                vec!["qp_create_ezfio", "--basis"],
                "Fatal error: exception Getopt.Error(\"Option --basis must have an argument\")\n",
            ),
            (
                vec!["qp_create_ezfio", "--basis", "x"],
                "Fatal error: exception Getopt.Error(\"Option --basis must have an argument\")\n",
            ),
            (
                vec!["qp_create_ezfio", "--au=x"],
                "Fatal error: exception Getopt.Error(\"Option --au does not accept argument\")\n",
            ),
            (
                vec!["qp_create_ezfio", "--help=1"],
                "Fatal error: exception Getopt.Error(\"Option --help does not accept argument\")\n",
            ),
            (
                vec!["qp_create_ezfio", "-b"],
                "Fatal error: exception Getopt.Error(\"Option -b must have an argument\")\n",
            ),
            (
                vec!["qp_create_ezfio", "-b", "-x"],
                "Fatal error: exception Getopt.Error(\"Option -b must have an argument\")\n",
            ),
        ] {
            let (error, help) = parse_err(&args);
            assert_eq!(error.render(), expected, "args {args:?}");
            assert!(help.is_empty(), "help must not print before a scan error");
        }
    }

    #[test]
    fn empty_next_token_is_invalid_argument() {
        let (error, _) = parse_err(&["qp_create_ezfio", "-b", ""]);
        assert_eq!(
            error.render(),
            "Fatal error: exception Invalid_argument(\"index out of bounds\")\n"
        );
        let (error, _) = parse_err(&["qp_create_ezfio", "-o", "", "-b", "x", "f.xyz"]);
        assert_eq!(
            error.render(),
            "Fatal error: exception Invalid_argument(\"index out of bounds\")\n"
        );
    }

    #[test]
    fn mandatory_basis_is_checked_before_conversions_and_input() {
        let (error, _) = parse_err(&["qp_create_ezfio", "-a"]);
        assert_eq!(
            error.render(),
            "qp_create_ezfio: --basis option is missing.\n"
        );
        // The mandatory check precedes the charge parser.
        let (error, _) = parse_err(&["qp_create_ezfio", "-c", "abc"]);
        assert_eq!(
            error.render(),
            "qp_create_ezfio: --basis option is missing.\n"
        );
        let (error, _) = parse_err(&["qp_create_ezfio", "-b", "x"]);
        assert_eq!(error.render(), "qp_create_ezfio: input file is missing\n");
    }

    #[test]
    fn conversions_run_charge_ghost_multiplicity() {
        let (error, _) = parse_err(&["qp_create_ezfio", "-b", "x", "-c", "abc", "-g", "x"]);
        assert_eq!(error.render(), "qp_create_ezfio: int_of_string\n");
        let (error, _) = parse_err(&[
            "qp_create_ezfio",
            "-b",
            "x",
            "-g",
            "x",
            "-m",
            "y",
            "-c",
            "abc",
        ]);
        assert_eq!(error.render(), "qp_create_ezfio: int_of_string\n");
        let (error, _) = parse_err(&["qp_create_ezfio", "-b", "x", "-m", "y", "-g", "x"]);
        assert_eq!(error.render(), "qp_create_ezfio: float_of_string\n");
        let (error, _) = parse_err(&["qp_create_ezfio", "-b", "x", "-m", "y"]);
        assert_eq!(error.render(), "qp_create_ezfio: int_of_string\n");
        // `--charge=` reaches the driver as an empty string, whose `.[0]`
        // is the Invalid_argument (unlike short `-c ""`, rejected by getopt).
        let (error, _) = parse_err(&["qp_create_ezfio", "-b", "x", "--charge=", "f.xyz"]);
        assert_eq!(
            error.render(),
            "Fatal error: exception Invalid_argument(\"index out of bounds\")\n"
        );
        let (error, _) = parse_err(&["qp_create_ezfio", "-b", "x", "--ghost=", "f.xyz"]);
        assert_eq!(error.render(), "qp_create_ezfio: float_of_string\n");
        let (error, _) = parse_err(&["qp_create_ezfio", "-b", "x", "--multiplicity=", "f.xyz"]);
        assert_eq!(error.render(), "qp_create_ezfio: int_of_string\n");
    }

    #[test]
    fn charge_m_prefix_negates() {
        let base = ["qp_create_ezfio", "-b", "x", "f.xyz"];
        let charge = |value: &'static str| {
            let mut args = base.to_vec();
            args.extend(["-c", value]);
            parse_ok(&args).charge
        };
        assert_eq!(charge("5"), 5);
        assert_eq!(charge("m1"), -1);
        assert_eq!(charge("m0x10"), -16);
        // The `m` prefix negates the parsed integer, including a signed one:
        // `m+1` is -1 and `m-1` is +1.
        assert_eq!(charge("m+1"), -1);
        assert_eq!(charge("m-1"), 1);
        assert_eq!(charge("+5"), 5);
        assert_eq!(charge("0x10"), 16);
    }

    #[test]
    fn help_prints_and_parsing_continues() {
        let argv: Vec<String> = ["qp_create_ezfio", "-h", "-b", "cc-pvdz", "h2.xyz"]
            .iter()
            .map(|argument| (*argument).to_owned())
            .collect();
        let mut help = Vec::new();
        let options = parse(&argv, &mut help).expect("help does not stop parsing");
        assert_eq!(help, HELP_TEMPLATE.as_bytes());
        assert_eq!(options.basis, "cc-pvdz");

        // `-h` alone: help then the missing-file failure.
        let (error, help) = parse_err(&["qp_create_ezfio", "-h"]);
        assert_eq!(help, HELP_TEMPLATE.as_bytes());
        assert_eq!(error.render(), "qp_create_ezfio: input file is missing\n");

        // Repeated help prints once (`Hashtbl.replace`).
        let (error, help) = parse_err(&["qp_create_ezfio", "-hh"]);
        assert_eq!(help, HELP_TEMPLATE.as_bytes());
        assert_eq!(error.render(), "qp_create_ezfio: input file is missing\n");

        // Help skips the mandatory basis check but not the conversions.
        let (error, help) = parse_err(&["qp_create_ezfio", "-h", "-c", "abc"]);
        assert_eq!(help, HELP_TEMPLATE.as_bytes());
        assert_eq!(error.render(), "qp_create_ezfio: int_of_string\n");

        // A scan error after `-h` suppresses the help: handlers are not run
        // when the subsequent token fails to parse.
        let (error, help) = parse_err(&["qp_create_ezfio", "-h", "-Z"]);
        assert!(help.is_empty());
        assert_eq!(
            error.render(),
            "Fatal error: exception Getopt.Error(\"Unknown option -Z\")\n"
        );
    }

    #[test]
    fn help_with_help_skips_the_mandatory_basis() {
        let options = parse_ok(&["qp_create_ezfio", "-h", "-q", "f.xyz"]);
        assert_eq!(options.basis, "");
        assert!(options.qmc);
    }

    #[test]
    fn anonymous_arguments_and_terminator() {
        // Extra anonymous arguments are ignored; the first one is FILE.
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", "a.xyz", "b.xyz"]).input,
            PathBuf::from("a.xyz")
        );
        // `--` terminates option parsing.
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", "--", "-h"]).input,
            PathBuf::from("-h")
        );
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", "--", "-b", "f.xyz"]).input,
            PathBuf::from("-b")
        );
        // A lone `-` and an empty token are anonymous.
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", "-"]).input,
            PathBuf::from("-")
        );
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", ""]).input,
            PathBuf::from("")
        );
    }

    #[test]
    fn pseudo_and_flags() {
        let options = parse_ok(&["qp_create_ezfio", "-b", "x", "-p", "bfd", "f.xyz"]);
        assert_eq!(options.pseudo.as_deref(), Some("bfd"));
        assert_eq!(
            parse_ok(&["qp_create_ezfio", "-b", "x", "--pseudo=", "f.xyz"]).pseudo,
            Some(String::new())
        );
    }
}
