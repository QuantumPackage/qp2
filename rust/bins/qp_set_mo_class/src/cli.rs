//! The reference CLI of `qp_set_mo_class`: a statement-by-statement port of
//! `ocaml/Command_line.ml` + `Getopt.parse_cmdline` (opam `getopt`
//! 20120615) for this tool's option set.
//!
//! The reference parser is **not** clap-compatible; the quirks below are
//! measured against `qp_set_mo_class.native` (2026-11-05) and reproduced:
//!
//! * `-c ""` (two tokens): the lookahead `args.(no+1).[0] <> '-'` is
//!   evaluated on the *empty* next argument and raises the uncaught
//!   `Invalid_argument("index out of bounds")` — before any range parsing
//!   (and therefore before the driver's `check` wrapper can print the
//!   `Number of MOs:` line).
//! * `-c=1`: the value is `"=1"` (short options have no `=` syntax);
//! * `-c` with no next argument: `Getopt.Error("Option -c must have an
//!   argument")` (`"--core"` for the long form);
//! * `--help=1`: `Getopt.Error("Option --help does not accept argument")`;
//! * `-qc 1`: `Getopt.Error("Only non-argument short-options can be
//!   concatenated (error with option c in -qc)")`;
//! * unknown option: `Getopt.Error("Unknown option -z")` /
//!   `Getopt.Error("Unknown option --nosuch")`;
//! * `-h` prints the help **and parsing continues** (`set_specs` prints
//!   once; the driver prints again when the anonymous argument is missing,
//!   so `-h` alone shows the help twice);
//! * a single `-` and empty anonymous arguments are anonymous arguments.
//!
//! The help block is the reference capture for `argv[0] =
//! "qp_set_mo_class"` with the program name substituted at the use site.

use std::io::Write;

use crate::error::{AppError, OcamlException};

/// The reference capture of the 639-byte help block (`argv[0]` bare).
const HELP_TEMPLATE: &str = include_str!("../tests/data/help_bare.txt");

/// The captured (bare-name) usage lines, as they appear in the template:
/// the only part of the help that depends on `argv[0]`.
const STATIC_USAGE: &str = "  qp_set_mo_class [-a <range>] [-c <range>] \
                             [-d <range>] [-h] [-i <range>]\n      \
                             [-q] [-v <range>] [--] EZFIO_DIR\n";

/// The usage words after the program name, in the reference's
/// `set_specs` order (sorted by short option).
const USAGE_WORDS: [&str; 9] = [
    "[-a <range>]",
    "[-c <range>]",
    "[-d <range>]",
    "[-h]",
    "[-i <range>]",
    "[-q]",
    "[-v <range>]",
    "[--]",
    "EZFIO_DIR",
];

/// The reference's usage box is a `Format` `@[<hov 4>` (`Command_line.ml`,
/// `help()`), so the line containing `argv[0]` re-wraps with the name.
/// The wrap rule was derived from argv0-length sweeps against renamed
/// copies of the reference binary (2026-11-05, lengths 15–89 for this
/// tool and 9–80 for `qp_tunnel`, both reproducing the same rule): the
/// first usage line starts at column 2, continuation lines are indented
/// 6 columns, and a word following a line of length `L` (separator space
/// at column `s = L + 1`, ending at `e = s + w`) is kept on the line
/// only when **both** `s ≤ 68` and `e ≤ 77`; when it wraps, the
/// separator space stays as a trailing space on the old line only when
/// `e ≤ 77` (the reference commits the space optimistically when the
/// word would fit the 77-column margin, then still wraps a short word
/// that starts past column 69). The two `"  "` guard lines around the
/// box render constant.
fn wrap_usage(program: &str) -> String {
    const MARGIN: usize = 77;
    const START: usize = 68;
    const CONTINUATION: &str = "      ";
    let mut lines = Vec::new();
    let mut line = format!("  {program}");
    for word in USAGE_WORDS {
        let space_col = line.len() + 1;
        let end_col = space_col + word.len();
        if space_col <= START && end_col <= MARGIN {
            line.push(' ');
            line.push_str(word);
        } else {
            if end_col <= MARGIN {
                line.push(' ');
            }
            lines.push(line);
            line = format!("{CONTINUATION}{word}");
        }
    }
    lines.push(line);
    lines.join("\n") + "\n"
}

/// Render the help block with `program` as `argv[0]` (the reference
/// prints `Sys.argv.(0)` on the `Usage:` line, re-wrapped).
#[must_use]
pub fn render_help(program: &str) -> String {
    let (before, after) = HELP_TEMPLATE
        .split_once(STATIC_USAGE)
        .expect("the captured help contains the bare usage block");
    format!("{before}{}{after}", wrap_usage(program))
}

/// One option in the reference's `set_specs` table (after the `-h` entry
/// `set_specs` prepends).
struct Spec {
    short: char,
    long: &'static str,
    takes_arg: bool,
}

/// The reference table, in source order (`qp_set_mo_class.ml:189-206`);
/// lookups are by name, so order only documents the source.
const SPECS: &[Spec] = &[
    Spec {
        short: 'c',
        long: "core",
        takes_arg: true,
    },
    Spec {
        short: 'i',
        long: "inact",
        takes_arg: true,
    },
    Spec {
        short: 'a',
        long: "act",
        takes_arg: true,
    },
    Spec {
        short: 'v',
        long: "virt",
        takes_arg: true,
    },
    Spec {
        short: 'd',
        long: "del",
        takes_arg: true,
    },
    Spec {
        short: 'q',
        long: "query",
        takes_arg: false,
    },
    Spec {
        short: 'h',
        long: "help",
        takes_arg: false,
    },
];

fn find_spec(short: char, long: &str) -> Option<&'static Spec> {
    SPECS.iter().find(|s| s.short == short || s.long == long)
}

/// The parsed command line, in the reference's terms.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Parsed {
    /// `-c` (the option is not `With_opt_arg`, so the value is present
    /// exactly when the option was seen).
    pub core: Option<String>,
    /// `-i`.
    pub inact: Option<String>,
    /// `-a`.
    pub act: Option<String>,
    /// `-v`.
    pub virt: Option<String>,
    /// `-d`.
    pub del: Option<String>,
    /// `-q`.
    pub query: bool,
    /// `-h` was seen (the help flag is printed by the driver, which then
    /// **continues**).
    pub help: bool,
    /// The anonymous arguments, in order.
    pub anonymous: Vec<String>,
}

impl Parsed {
    /// `set_specs`: the mandatory-option check. There are no mandatory
    /// **options** (the anonymous `EZFIO_DIR` is checked by the driver),
    /// so this is the help decision: print the help when `-h` was seen.
    ///
    /// # Errors
    ///
    /// An [`AppError`] when the help output itself fails.
    pub fn after_specs(&self, program: &str, out: &mut dyn Write) -> Result<(), AppError> {
        if self.help {
            out.write_all(render_help(program).as_bytes())
                .map_err(io_error)?;
        }
        Ok(())
    }
}

fn io_error(error: std::io::Error) -> AppError {
    AppError(OcamlException::SysError(error.to_string()))
}

/// `Getopt.parse`: walk `args` (the reference's `Sys.argv`, `argv[0]`
/// included; parsing starts at 1) exactly as `aux` does.
///
/// # Errors
///
/// The uncaught reference exceptions: `Getopt.Error(…)`, and the
/// `Invalid_argument("index out of bounds")` of the empty-next-argument
/// lookahead.
pub fn parse(args: &[String]) -> Result<Parsed, AppError> {
    let mut p = Parsed::default();
    let mut no = 1usize;
    let last = args.len().saturating_sub(1);
    let set = |p: &mut Parsed, spec: &Spec, value: &str| match spec.long {
        "core" => p.core = Some(value.to_owned()),
        "inact" => p.inact = Some(value.to_owned()),
        "act" => p.act = Some(value.to_owned()),
        "virt" => p.virt = Some(value.to_owned()),
        "del" => p.del = Some(value.to_owned()),
        "query" => p.query = true,
        "help" => p.help = true,
        _ => unreachable!("every spec is one of the five ranges, query or help"),
    };
    let unknown_short =
        |c: char| AppError(OcamlException::GetoptError(format!("Unknown option -{c}")));
    let unknown_long =
        |s: &str| AppError(OcamlException::GetoptError(format!("Unknown option --{s}")));

    while no <= last {
        let s = &args[no];
        let l = s.chars().count();
        if l == 0 {
            p.anonymous.push(String::new());
            no += 1;
            continue;
        }
        let bytes = s.as_bytes();
        if bytes[0] == b'-' {
            if l >= 2 && bytes[1] == b'-' {
                if l == 2 {
                    // `--`: consumed (not an anonymous argument itself);
                    // every following token is anonymous.
                    no += 1;
                    while no <= last {
                        p.anonymous.push(args[no].clone());
                        no += 1;
                    }
                    break;
                }
                match s[2..].find('=') {
                    Some(i) => {
                        let opt = &s[2..2 + i];
                        let value = &s[2 + i + 1..];
                        let spec = find_spec('\0', opt).ok_or_else(|| unknown_long(opt))?;
                        if !spec.takes_arg {
                            return Err(AppError(OcamlException::GetoptError(format!(
                                "Option --{} does not accept argument",
                                opt
                            ))));
                        }
                        set(&mut p, spec, value);
                        no += 1;
                    }
                    None => {
                        let opt = &s[2..];
                        let spec = find_spec('\0', opt).ok_or_else(|| unknown_long(opt))?;
                        if spec.takes_arg {
                            return Err(AppError(OcamlException::GetoptError(format!(
                                "Option --{opt} must have an argument"
                            ))));
                        }
                        set(&mut p, spec, "");
                        no += 1;
                    }
                }
            } else if l == 1 {
                // A single `-` is an anonymous argument.
                p.anonymous.push("-".to_owned());
                no += 1;
            } else {
                let opt = bytes[1] as char;
                let spec = find_spec(opt, "").ok_or_else(|| unknown_short(opt))?;
                if !spec.takes_arg {
                    set(&mut p, spec, "");
                    // Combined short flags: the reference looks up every
                    // following character in `find_short` — an unknown
                    // character raises `Unknown option -<c>` (the lookup
                    // itself, not the match arm), a known character that
                    // takes an argument hits the match's catch-all arm.
                    for c in s[2..].chars() {
                        match find_spec(c, "") {
                            None => {
                                return Err(unknown_short(c));
                            }
                            Some(next) if next.takes_arg => {
                                return Err(AppError(OcamlException::GetoptError(format!(
                                    "Only non-argument short-options can be concatenated (error with option {c} in {s})"
                                ))));
                            }
                            Some(next) => set(&mut p, next, ""),
                        }
                    }
                    no += 1;
                } else if l > 2 {
                    // Immediate argument: the whole remainder, `=` included.
                    set(&mut p, spec, &s[2..]);
                    no += 1;
                } else if no < last {
                    // Non-immediate argument. The reference evaluates
                    // `args.(no+1).[0] <> '-'`; on an *empty* next argument
                    // that index itself raises the uncaught
                    // `Invalid_argument("index out of bounds")`.
                    let next = &args[no + 1];
                    if next.is_empty() {
                        return Err(AppError(OcamlException::InvalidArgument(
                            "index out of bounds".into(),
                        )));
                    }
                    if next.starts_with('-') {
                        return Err(AppError(OcamlException::GetoptError(format!(
                            "Option -{opt} must have an argument"
                        ))));
                    }
                    set(&mut p, spec, next);
                    no += 2;
                } else {
                    return Err(AppError(OcamlException::GetoptError(format!(
                        "Option -{opt} must have an argument"
                    ))));
                }
            }
        } else {
            p.anonymous.push(s.clone());
            no += 1;
        }
    }
    Ok(p)
}

/// The driver's anonymous-argument decision (`qp_set_mo_class.ml:229-232`):
/// exactly one anonymous argument, otherwise the help is printed again and
/// the uncaught `Failure("EZFIO_FILE is missing")` is raised.
///
/// # Errors
///
/// [`AppError`] with the `Failure` when the argument is missing (or there
/// is more than one).
pub fn require_one_anonymous(
    parsed: &Parsed,
    program: &str,
    out: &mut dyn Write,
) -> Result<String, AppError> {
    match parsed.anonymous.as_slice() {
        [x] => Ok(x.clone()),
        _ => {
            out.write_all(render_help(program).as_bytes())
                .map_err(io_error)?;
            Err(AppError(OcamlException::Failure(
                "EZFIO_FILE is missing".into(),
            )))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_ok(args: &[&str]) -> Parsed {
        parse(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>()).expect("parse")
    }

    fn parse_err(args: &[&str]) -> OcamlException {
        parse(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>())
            .expect_err("must fail")
            .0
    }

    #[test]
    fn the_bare_usage_help_is_the_reference_capture() {
        let mut out: Vec<u8> = Vec::new();
        let parsed = Parsed {
            help: true,
            ..Default::default()
        };
        parsed.after_specs("qp_set_mo_class", &mut out).unwrap();
        assert_eq!(out, include_str!("../tests/data/help_bare.txt").as_bytes());
        // The invoked name is substituted on the Usage line.
        let mut out2: Vec<u8> = Vec::new();
        parsed
            .after_specs("/x/bin/qp_set_mo_class", &mut out2)
            .unwrap();
        assert!(
            out2.windows(b"/x/bin/qp_set_mo_class [-a <range>]".len())
                .any(|w| w == b"/x/bin/qp_set_mo_class [-a <range>]")
        );
    }

    #[test]
    fn the_usage_line_rewraps_with_long_argv0() {
        // Captures taken from renamed copies of the reference binary
        // (its help embeds `argv[0]` in a Format `hov` box).
        assert_eq!(
            render_help("./qp_set_mo_class_medium_name_xyz"),
            include_str!("../tests/data/help_medium.txt")
        );
        assert_eq!(
            render_help("./qp_set_mo_class_really_really_really_long_name_for_wrapping_tests"),
            include_str!("../tests/data/help_long.txt")
        );
        // And the bare name the capture was taken with reproduces it.
        assert_eq!(
            render_help("qp_set_mo_class"),
            include_str!("../tests/data/help_bare.txt")
        );
    }

    #[test]
    fn options_and_anonymous_args_are_collected() {
        let p = parse_ok(&["qp_set_mo_class", "-c", "1", "-i", "2-3", "ezfio"]);
        assert_eq!(p.core.as_deref(), Some("1"));
        assert_eq!(p.inact.as_deref(), Some("2-3"));
        assert_eq!(p.anonymous, vec!["ezfio".to_string()]);
        assert!(!p.query && !p.help);
    }

    #[test]
    fn a_single_dash_is_anonymous_argument() {
        let p = parse_ok(&["qp_set_mo_class", "-"]);
        assert_eq!(p.anonymous, vec!["-".to_string()]);
    }

    #[test]
    fn dash_dash_switches_everything_to_anonymous() {
        let p = parse_ok(&["qp_set_mo_class", "--", "-c", "1"]);
        assert_eq!(p.anonymous, vec!["-c".to_string(), "1".to_string()]);
        assert_eq!(p.core, None);
    }

    #[test]
    fn long_options_take_the_eq_value() {
        let p = parse_ok(&["qp_set_mo_class", "--core=5", "ezfio"]);
        assert_eq!(p.core.as_deref(), Some("5"));
    }

    #[test]
    fn short_options_take_the_immediate_remainder_including_eq() {
        // The reference has no `=` syntax for short options: `-c=1` means
        // the value is `"=1"`.
        let p = parse_ok(&["qp_set_mo_class", "-c=1", "ezfio"]);
        assert_eq!(p.core.as_deref(), Some("=1"));
        let p = parse_ok(&["qp_set_mo_class", "-c3", "ezfio"]);
        assert_eq!(p.core.as_deref(), Some("3"));
    }

    #[test]
    fn combined_short_flags_must_all_be_no_argument() {
        let p = parse_ok(&["qp_set_mo_class", "-qh", "ezfio"]);
        assert!(p.query && p.help);
        assert_eq!(
            parse_err(&["qp_set_mo_class", "-qc", "1", "ezfio"]),
            OcamlException::GetoptError(
                "Only non-argument short-options can be concatenated (error with option c in -qc)"
                    .into()
            )
        );
        // An *unknown* character in the tail is reported by the lookup
        // itself, not the catch-all arm (measured: `-qz` → Unknown option).
        assert_eq!(
            parse_err(&["qp_set_mo_class", "-qz", "ezfio"]),
            OcamlException::GetoptError("Unknown option -z".into())
        );
    }

    #[test]
    fn the_missing_value_messages_name_the_option_form_used() {
        assert_eq!(
            parse_err(&["qp_set_mo_class", "-c"]),
            OcamlException::GetoptError("Option -c must have an argument".into())
        );
        assert_eq!(
            parse_err(&["qp_set_mo_class", "--core"]),
            OcamlException::GetoptError("Option --core must have an argument".into())
        );
        // A value that looks like an option is not taken.
        assert_eq!(
            parse_err(&["qp_set_mo_class", "-c", "-q", "ezfio"]),
            OcamlException::GetoptError("Option -c must have an argument".into())
        );
    }

    #[test]
    fn no_argument_options_reject_eq_values() {
        assert_eq!(
            parse_err(&["qp_set_mo_class", "--help=1"]),
            OcamlException::GetoptError("Option --help does not accept argument".into())
        );
    }

    #[test]
    fn unknown_options_are_getopt_errors() {
        assert_eq!(
            parse_err(&["qp_set_mo_class", "-z", "ezfio"]),
            OcamlException::GetoptError("Unknown option -z".into())
        );
        assert_eq!(
            parse_err(&["qp_set_mo_class", "--nosuch", "ezfio"]),
            OcamlException::GetoptError("Unknown option --nosuch".into())
        );
    }

    #[test]
    fn the_empty_next_argument_raises_index_out_of_bounds() {
        // `-c ""`: the reference's `args.(no+1).[0]` lookahead indexes an
        // empty string.
        assert_eq!(
            parse_err(&["qp_set_mo_class", "-c", "", "ezfio"]),
            OcamlException::InvalidArgument("index out of bounds".into())
        );
    }

    #[test]
    fn the_next_argument_is_taken_verbatim() {
        // A non-dash, non-empty next token is the value, whatever it is
        // (e.g. it can even look like a directory name).
        let p = parse_ok(&["qp_set_mo_class", "-c", "ezfio"]);
        assert_eq!(p.core.as_deref(), Some("ezfio"));
    }

    #[test]
    fn one_anonymous_argument_passes_more_fail_with_help_and_failure() {
        let mut out: Vec<u8> = Vec::new();
        let parsed = parse_ok(&["qp_set_mo_class", "ezfio"]);
        assert_eq!(
            require_one_anonymous(&parsed, "qp_set_mo_class", &mut out).unwrap(),
            "ezfio"
        );
        assert!(out.is_empty());

        let mut out: Vec<u8> = Vec::new();
        let parsed = parse_ok(&["qp_set_mo_class"]);
        let err = require_one_anonymous(&parsed, "qp_set_mo_class", &mut out).unwrap_err();
        assert_eq!(
            err.0,
            OcamlException::Failure("EZFIO_FILE is missing".into())
        );
        // The help is printed a second time before the failure.
        assert_eq!(out, include_str!("../tests/data/help_bare.txt").as_bytes());

        let mut out: Vec<u8> = Vec::new();
        let parsed = parse_ok(&["qp_set_mo_class", "a", "b"]);
        let err = require_one_anonymous(&parsed, "qp_set_mo_class", &mut out).unwrap_err();
        assert_eq!(
            err.0,
            OcamlException::Failure("EZFIO_FILE is missing".into())
        );
    }
}
