//! The reference CLI of `qp_tunnel`: a statement-by-statement port of
//! `ocaml/Command_line.ml` + `Getopt.parse_cmdline` (opam `getopt`
//! 20120615) for this tool's option set (`qp_tunnel.ml:14-45`: `-g`
//! and the anonymous `(EZFIO_DIR|ADDRESS)`, with `-h` prepended by
//! `set_specs`).
//!
//! The reference parser is **not** clap-compatible; the quirks below are
//! measured against `qp_tunnel.native` (2026-11-05) and reproduced:
//!
//! * `--get-input=1`: `Getopt.Error("Option --get-input does not accept
//!   argument")` (a value on a `Without_arg` long option);
//! * `--nosuch=1`: the `=1` is ignored in the lookup —
//!   `Getopt.Error("Unknown option --nosuch")`;
//! * `--h`: `Getopt.Error("Unknown option --h")` (one-character long
//!   names are unknown);
//! * `-h=x`: the walk continues after the flag, so `=` is looked up as a
//!   short option — `Getopt.Error("Unknown option -=")`;
//! * `-gh`: both flags are set (concatenated no-argument shorts);
//! * `-g1`: `Getopt.Error("Unknown option -1")`;
//! * `-hx`: `Getopt.Error("Unknown option -x")`;
//! * `--` is **consumed** (not an anonymous argument itself); every
//!   following token is anonymous;
//! * a single `-` and empty anonymous arguments are anonymous arguments;
//! * `-h` prints the help **and parsing continues** (`set_specs` prints
//!   once; the driver prints again when the anonymous argument is missing,
//!   so `-h` alone shows the help twice);
//! * there are no `With_arg` options, so the missing-argument and
//!   concatenation-with-argument errors of the full `aux` walk are
//!   unreachable.

use qp_tunnel::OcamlException;

/// The reference capture of the 481-byte help block (`argv[0]` bare).
const HELP_TEMPLATE: &str = include_str!("../tests/data/help_qp_tunnel.txt");

/// The captured (bare-name) usage line, as it appears in the template:
/// the only part of the help that depends on `argv[0]`.
const STATIC_USAGE: &str = "  qp_tunnel [-g] [-h] [--] (EZFIO_DIR|ADDRESS)\n";

/// The usage words after the program name, in the reference's
/// `set_specs` order (sorted by short option).
const USAGE_WORDS: [&str; 4] = ["[-g]", "[-h]", "[--]", "(EZFIO_DIR|ADDRESS)"];

/// The reference's usage box is a `Format` `@[<hov 4>` (`Command_line.ml`,
/// `help()`), so the line containing `argv[0]` re-wraps with the name.
/// The wrap rule was derived from argv0-length sweeps against renamed
/// copies of the reference binary (2026-11-05, lengths 9–80 for tunnel
/// and 15–89 for `qp_set_mo_class`, both reproducing the same rule):
/// the first usage line starts at column 2, continuation lines are
/// indented 6 columns, and a word following a line of length `L`
/// (separator space at column `s = L + 1`, ending at `e = s + w`) is
/// kept on the line only when **both** `s ≤ 68` and `e ≤ 77`; when it
/// wraps, the separator space stays as a trailing space on the old line
/// only when `e ≤ 77` (the reference commits the space optimistically
/// when the word would fit the 77-column margin, then still wraps a
/// short word that starts past column 69). The two `"  "` guard lines
/// around the box render constant.
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
        .expect("the captured help contains the bare usage line");
    format!("{before}{}{after}", wrap_usage(program))
}

/// One option in the reference's `set_specs` table (after the `-h` entry
/// `set_specs` prepends).
struct Spec {
    short: char,
    long: &'static str,
    takes_arg: bool,
}

/// The reference table, in source order (`qp_tunnel.ml:19-29`); lookups
/// are by name, so order only documents the source.
const SPECS: &[Spec] = &[
    Spec {
        short: 'g',
        long: "get-input",
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
    /// `-g` was seen.
    pub get_input: bool,
    /// `-h` was seen (the help flag is printed by `set_specs`, which then
    /// **continues**).
    pub help: bool,
    /// The anonymous arguments, in order.
    pub anonymous: Vec<String>,
}

/// Parse `args` (the reference's `Sys.argv`, `argv[0]` included; parsing
/// starts at 1) exactly as `Getopt.parse` does.
///
/// # Errors
///
/// The uncaught reference exceptions: `Getopt.Error(…)`, and the
/// `Invalid_argument("index out of bounds")` of the empty-next-argument
/// lookahead (unreachable here: no `With_arg` option can be last).
pub fn parse(args: &[String]) -> Result<Parsed, OcamlException> {
    let mut p = Parsed::default();
    let mut no = 1usize;
    let last = args.len().saturating_sub(1);
    let set = |p: &mut Parsed, spec: &Spec| match spec.long {
        "get-input" => p.get_input = true,
        "help" => p.help = true,
        _ => unreachable!("every spec is get-input or help"),
    };
    let unknown_short = |c: char| OcamlException::GetoptError(format!("Unknown option -{c}"));
    let unknown_long = |s: &str| OcamlException::GetoptError(format!("Unknown option --{s}"));

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
                        let spec = find_spec('\0', opt).ok_or_else(|| unknown_long(opt))?;
                        if !spec.takes_arg {
                            return Err(OcamlException::GetoptError(format!(
                                "Option --{opt} does not accept argument"
                            )));
                        }
                        set(&mut p, spec);
                        no += 1;
                    }
                    None => {
                        let opt = &s[2..];
                        let spec = find_spec('\0', opt).ok_or_else(|| unknown_long(opt))?;
                        if spec.takes_arg {
                            return Err(OcamlException::GetoptError(format!(
                                "Option --{opt} must have an argument"
                            )));
                        }
                        set(&mut p, spec);
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
                    set(&mut p, spec);
                    // Combined short flags: the reference looks up every
                    // following character in `find_short` — an unknown
                    // character raises `Unknown option -<c>` (the lookup
                    // itself, not the match arm), a known character that
                    // takes an argument hits the match's catch-all arm
                    // (unreachable: this table has none).
                    for c in s[2..].chars() {
                        match find_spec(c, "") {
                            None => {
                                return Err(unknown_short(c));
                            }
                            Some(next) if next.takes_arg => {
                                return Err(OcamlException::GetoptError(format!(
                                    "Only non-argument short-options can be concatenated (error with option {c} in {s})"
                                )));
                            }
                            Some(next) => set(&mut p, next),
                        }
                    }
                    no += 1;
                } else if no < last {
                    // Unreachable in this table (no `With_arg` option);
                    // kept so the walk mirrors `aux` statement for
                    // statement.
                    let next = &args[no + 1];
                    if next.is_empty() {
                        return Err(OcamlException::InvalidArgument(
                            "index out of bounds".into(),
                        ));
                    }
                    if next.starts_with('-') {
                        return Err(OcamlException::GetoptError(format!(
                            "Option -{opt} must have an argument"
                        )));
                    }
                    set(&mut p, spec);
                    no += 2;
                } else {
                    return Err(OcamlException::GetoptError(format!(
                        "Option -{opt} must have an argument"
                    )));
                }
            }
        } else {
            p.anonymous.push(s.clone());
            no += 1;
        }
    }
    Ok(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(argv: &[&str]) -> Vec<String> {
        argv.iter().map(|s| s.to_string()).collect()
    }

    fn parse_ok(argv: &[&str]) -> Parsed {
        parse(&args(argv)).expect("the case parses in the reference")
    }

    fn parse_fatal(argv: &[&str]) -> OcamlException {
        parse(&args(argv)).expect_err("the case fails in the reference")
    }

    #[test]
    fn the_bare_help_is_the_reference_capture() {
        assert_eq!(
            render_help("qp_tunnel"),
            include_str!("../tests/data/help_qp_tunnel.txt")
        );
    }

    #[test]
    fn the_usage_line_rewraps_with_long_argv0() {
        // Captures taken from renamed copies of the reference binary
        // (its help embeds `argv[0]` in a Format `hov` box).
        assert_eq!(
            render_help("./qp_tunnel_medium_name_xyz"),
            include_str!("../tests/data/help_medium.txt")
        );
        assert_eq!(
            render_help("./qp_tunnel_really_really_really_long_name_for_wrapping_tests"),
            include_str!("../tests/data/help_long.txt")
        );
    }

    #[test]
    fn flags_and_anonymous_args() {
        let p = parse_ok(&["qp_tunnel", "-g", "-h", "ezfio"]);
        assert!(p.get_input);
        assert!(p.help);
        assert_eq!(p.anonymous, vec!["ezfio".to_owned()]);

        let p = parse_ok(&["qp_tunnel", "host:1"]);
        assert!(!p.get_input);
        assert!(!p.help);
        assert_eq!(p.anonymous, vec!["host:1".to_owned()]);
    }

    #[test]
    fn concatenated_no_argument_flags() {
        let p = parse_ok(&["qp_tunnel", "-gh", "ezfio"]);
        assert!(p.get_input);
        assert!(p.help);
        assert_eq!(p.anonymous, vec!["ezfio".to_owned()]);

        let p = parse_ok(&["qp_tunnel", "-hg"]);
        assert!(p.get_input);
        assert!(p.help);
    }

    #[test]
    fn unknown_options() {
        assert_eq!(
            parse_fatal(&["qp_tunnel", "-z", "ezfio"]),
            OcamlException::GetoptError("Unknown option -z".into())
        );
        assert_eq!(
            parse_fatal(&["qp_tunnel", "--nosuch", "ezfio"]),
            OcamlException::GetoptError("Unknown option --nosuch".into())
        );
        // The `=value` part is not part of the lookup.
        assert_eq!(
            parse_fatal(&["qp_tunnel", "--nosuch=1", "ezfio"]),
            OcamlException::GetoptError("Unknown option --nosuch".into())
        );
        // One-character long names are unknown (there is no `h` long).
        assert_eq!(
            parse_fatal(&["qp_tunnel", "--h", "ezfio"]),
            OcamlException::GetoptError("Unknown option --h".into())
        );
    }

    #[test]
    fn a_value_on_a_no_argument_option() {
        assert_eq!(
            parse_fatal(&["qp_tunnel", "--get-input=1", "ezfio"]),
            OcamlException::GetoptError("Option --get-input does not accept argument".into())
        );
        // Same error with the flag spelled as a short cluster.
        assert_eq!(
            parse_fatal(&["qp_tunnel", "-g=x", "ezfio"]),
            OcamlException::GetoptError("Unknown option -=".into())
        );
    }

    #[test]
    fn the_cluster_walk_continues_after_a_flag() {
        // `-h=x`: `=` is looked up as a short option.
        assert_eq!(
            parse_fatal(&["qp_tunnel", "-h=x", "ezfio"]),
            OcamlException::GetoptError("Unknown option -=".into())
        );
        assert_eq!(
            parse_fatal(&["qp_tunnel", "-hx", "ezfio"]),
            OcamlException::GetoptError("Unknown option -x".into())
        );
        assert_eq!(
            parse_fatal(&["qp_tunnel", "-g1", "ezfio"]),
            OcamlException::GetoptError("Unknown option -1".into())
        );
    }

    #[test]
    fn double_dash_is_consumed_and_switches_to_anonymous() {
        let p = parse_ok(&["qp_tunnel", "--", "ezfio"]);
        assert_eq!(p.anonymous, vec!["ezfio".to_owned()]);

        let p = parse_ok(&["qp_tunnel", "--", "ezfio", "extra"]);
        assert_eq!(p.anonymous, vec!["ezfio".to_owned(), "extra".to_owned()]);

        // Options do not work after `--`.
        let p = parse_ok(&["qp_tunnel", "--", "-g"]);
        assert!(!p.get_input);
        assert_eq!(p.anonymous, vec!["-g".to_owned()]);
    }

    #[test]
    fn a_single_dash_and_empty_strings_are_anonymous() {
        let p = parse_ok(&["qp_tunnel", "-", ""]);
        assert_eq!(p.anonymous, vec!["-".to_owned(), String::new()]);
    }

    #[test]
    fn the_long_flag_is_recognised() {
        let p = parse_ok(&["qp_tunnel", "--get-input", "ezfio"]);
        assert!(p.get_input);
        assert_eq!(p.anonymous, vec!["ezfio".to_owned()]);
        let p = parse_ok(&["qp_tunnel", "--help", "ezfio"]);
        assert!(p.help);
    }
}
