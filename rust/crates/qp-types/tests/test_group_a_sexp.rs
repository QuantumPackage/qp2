//! Byte-level tests of [`qp_types::sexplib_compat`], the module that
//! defines the md5 contract.
//!
//! The reference values are not free choices:
//!
//! * the printer cases come from `sexplib0/sexp.ml:57-119, 253-276`
//!   (`must_escape`, `escaped`, `to_buffer_mach`);
//! * the float cases come from `sexplib0/sexp_conv.ml:19-23`
//!   (`%.15G`, then `%.17G` when the short form does not read back) and
//!   from C's definition of `%e`, `%f` and `%G`, computed independently
//!   with `printf`/Python on the same doubles.
//!
//! Every expected string is a literal, so a change to the printer shows
//! up as a diff rather than as a mysterious md5 change.

use qp_types::sexplib_compat::{
    Sexp, escape_atom, float_to_sexp, float_to_string, format_e, format_f, format_g, must_escape,
};

// ---------------------------------------------------------------------
// Printer spacing
// ---------------------------------------------------------------------

#[test]
fn sibling_atoms_are_separated_but_lists_are_not() {
    let a = || Sexp::atom("a");
    let b = || Sexp::atom("b");

    // Two unescaped sibling atoms: one space.
    assert_eq!(Sexp::list(vec![a(), b()]).to_string(), "(a b)");
    assert_eq!(
        Sexp::list(vec![a(), b(), Sexp::atom("c")]).to_string(),
        "(a b c)"
    );

    // A list never takes a space and clears the flag on both sides.
    assert_eq!(
        Sexp::list(vec![a(), Sexp::list(vec![b()])]).to_string(),
        "(a(b))"
    );
    assert_eq!(
        Sexp::list(vec![Sexp::list(vec![a()]), b()]).to_string(),
        "((a)b)"
    );
    assert_eq!(
        Sexp::list(vec![a(), Sexp::list(vec![]), b()]).to_string(),
        "(a()b)"
    );
    assert_eq!(
        Sexp::list(vec![Sexp::list(vec![a(), b()]), Sexp::list(vec![a()])]).to_string(),
        "((a b)(a))"
    );

    // An escaped atom sets the flag to false, on both sides.
    assert_eq!(
        Sexp::list(vec![a(), Sexp::atom("b c"), Sexp::atom("d")]).to_string(),
        "(a\"b c\"d)"
    );
    assert_eq!(
        Sexp::list(vec![Sexp::atom("a b"), b()]).to_string(),
        "(\"a b\"b)"
    );

    // Top-level nodes: no leading space, whatever the previous call did.
    assert_eq!(Sexp::atom("a").to_string(), "a");
    assert_eq!(Sexp::list(vec![]).to_string(), "()");
    assert_eq!(Sexp::atom("").to_string(), "\"\"");
}

#[test]
fn the_derived_pair_shape_matches_a_shell() {
    // `[@@deriving sexp]` writes a pair as a two-element list with no
    // extra level: the shell of `Gto.ml` prints `(((sym S)(lc …))1)`.
    let sym = Sexp::pair(Sexp::atom("sym"), Sexp::atom("S"));
    let lc = Sexp::pair(Sexp::atom("lc"), Sexp::list(vec![]));
    let shell = Sexp::list(vec![Sexp::pair(sym, lc), Sexp::atom("1")]);
    assert_eq!(shell.to_string(), "(((sym S)(lc()))1)");

    // And the record field order is the declaration order.
    let point = Sexp::list(vec![
        Sexp::pair(Sexp::atom("x"), float_to_sexp(1.0)),
        Sexp::pair(Sexp::atom("y"), float_to_sexp(2.0)),
        Sexp::pair(Sexp::atom("z"), float_to_sexp(3.0)),
    ]);
    assert_eq!(point.to_string(), "((x 1)(y 2)(z 3))");
}

// ---------------------------------------------------------------------
// Escaping
// ---------------------------------------------------------------------

#[test]
fn must_escape_follows_the_sexplib_rule() {
    // Empty, and the five always-escaped bytes.
    assert!(must_escape(""));
    for s in [
        "\"", "(", ")", ";", "\\", "a(b", "a)b", "a;b", "a\"b", "a\\b",
    ] {
        assert!(must_escape(s), "{s:?} must be escaped");
    }

    // Ordinary atoms, including the lone `|` and `#`.
    for s in ["a", "abc", "H", "6-31g_star", "a|b", "a#b", "|", "#"] {
        assert!(!must_escape(s), "{s:?} must not be escaped");
    }

    // `|` and `#` only escape as an **adjacent** pair; the backward scan
    // in sexplib does not look further than the immediate predecessor.
    // Checked against sexplib0: `"a#b|c"` prints as `a#b|c`.
    assert!(must_escape("#|"));
    assert!(must_escape("|#"));
    assert!(must_escape("a#|b"));
    assert!(must_escape("a|#b"));
    assert!(!must_escape("a#b|c"));
    assert!(!must_escape("a|b#c"));

    // Space and everything outside printable ASCII.
    assert!(must_escape("a b"));
    assert!(must_escape("a\tb"));
    assert!(must_escape("a\nb"));
    assert!(must_escape("\u{7f}"));
    assert!(must_escape("é"));
    assert!(!must_escape("~"));
}

#[test]
fn escapes_use_the_sexplib_spelling() {
    assert_eq!(escape_atom(""), "\"\"");
    assert_eq!(escape_atom("abc"), "\"abc\"");
    assert_eq!(escape_atom("a b"), "\"a b\"");
    assert_eq!(escape_atom("a\"b"), "\"a\\\"b\"");
    assert_eq!(escape_atom("a\\b"), "\"a\\\\b\"");
    assert_eq!(escape_atom("a\nb"), "\"a\\nb\"");
    assert_eq!(escape_atom("a\tb"), "\"a\\tb\"");
    assert_eq!(escape_atom("a\rb"), "\"a\\rb\"");
    assert_eq!(escape_atom("a\u{8}b"), "\"a\\bb\"");
    // Everything else is `\` followed by three **decimal** digits, and a
    // UTF-8 character is escaped byte by byte: é is 195, 169.
    assert_eq!(escape_atom("é"), "\"\\195\\169\"");
    assert_eq!(escape_atom("\u{1}"), "\"\\001\"");
    assert_eq!(escape_atom("\u{7f}"), "\"\\127\"");
    assert_eq!(escape_atom("\u{80}"), "\"\\194\\128\"");
}

// ---------------------------------------------------------------------
// C printf number formats
// ---------------------------------------------------------------------

#[test]
fn e_has_a_signed_two_digit_exponent() {
    assert_eq!(format_e(1.0, 15, false), "1.000000000000000e+00");
    assert_eq!(format_e(-1.5e-7, 2, false), "-1.50e-07");
    assert_eq!(format_e(1.0, 0, true), "1E+00");
    assert_eq!(format_e(1e-15, 15, false), "1.000000000000000e-15");
    assert_eq!(format_e(1e100, 1, true), "1.0E+100");
    assert_eq!(format_e(0.0, 3, false), "0.000e+00");
}

#[test]
fn f_is_fixed_with_the_requested_decimals() {
    assert_eq!(format_f(14710.0, 6), "14710.000000");
    assert_eq!(format_f(-0.002332, 6), "-0.002332");
    assert_eq!(format_f(1.0 / 3.0, 8), "0.33333333");
    assert_eq!(format_f(0.0, 0), "0");
    assert_eq!(format_f(-0.0, 3), "-0.000");
}

#[test]
fn g_rounds_once_and_drops_trailing_zeros() {
    // One significant digit, then the trailing-zero rule.
    assert_eq!(format_g(14710.0, 15, true), "14710");
    assert_eq!(format_g(0.0, 15, true), "0");
    assert_eq!(format_g(-0.0, 15, true), "-0");
    assert_eq!(format_g(1.0, 17, false), "1");

    // The fixed/exponential switch is `-4 <= X < precision`.
    assert_eq!(format_g(1e-4, 15, true), "0.0001");
    assert_eq!(format_g(1e-5, 15, true), "1E-05");
    assert_eq!(format_g(1e-5, 17, false), "1.0000000000000001e-05");
    assert_eq!(format_g(1e15, 17, false), "1000000000000000");
    // At 17 digits the exponent 16 is still below the precision.
    assert_eq!(format_g(1e16, 17, false), "10000000000000000");
    assert_eq!(format_g(1e16, 15, false), "1e+16");

    // Rounding at the requested precision (Python's `%.2g`).
    assert_eq!(format_g(1234.5, 2, false), "1.2e+03");
    assert_eq!(format_g(0.0001234, 2, false), "0.00012");
    assert_eq!(format_g(0.00001234, 2, false), "1.2e-05");

    // A precision of zero means one significant digit.
    assert_eq!(format_g(1234.5, 0, false), "1e+03");
}

#[test]
fn non_finite_values_use_the_c_spelling() {
    assert_eq!(format_g(f64::INFINITY, 15, true), "INF");
    assert_eq!(format_g(f64::NEG_INFINITY, 15, true), "-INF");
    assert_eq!(format_g(f64::NAN, 15, true), "NAN");
    assert_eq!(format_e(f64::NEG_INFINITY, 2, false), "-inf");
    assert_eq!(format_f(f64::NAN, 3), "nan");
}

// ---------------------------------------------------------------------
// The sexplib float rule
// ---------------------------------------------------------------------

#[test]
fn float_atoms_use_15g_unless_it_does_not_round_trip() {
    // Short forms that read back exactly.
    assert_eq!(float_to_string(0.0), "0");
    assert_eq!(float_to_string(-0.0), "-0");
    assert_eq!(float_to_string(0.1), "0.1");
    assert_eq!(float_to_string(2.0), "2");
    assert_eq!(float_to_string(14710.0), "14710");
    assert_eq!(float_to_string(-0.002332), "-0.002332");
    assert_eq!(float_to_string(1e-4), "0.0001");
    assert_eq!(float_to_string(1e-5), "1E-05");
    assert_eq!(float_to_string(1e15), "1E+15");
    assert_eq!(float_to_string(1e16), "1E+16");
    // Only one double lives near 5e-324 (the smallest subnormal, not
    // `f64::MIN_POSITIVE`, which is the smallest normal), so 15 digits
    // already suffice.
    assert_eq!(float_to_string(f64::from_bits(1)), "4.94065645841247E-324");
    // The smallest normal needs the 17-digit fallback.
    assert_eq!(
        float_to_string(f64::MIN_POSITIVE),
        "2.2250738585072014E-308"
    );

    // Fallbacks: 15 significant digits lose the value.
    assert_eq!(float_to_string(0.1 + 0.2), "0.30000000000000004");
    assert_eq!(float_to_string(1.0 / 3.0), "0.33333333333333331");
    assert_eq!(float_to_string(f64::MAX), "1.7976931348623157E+308");

    // Every corpus-style continuation is exact at 15 digits; check that
    // a sample of the values in `basis-md5.jsonl` keeps its short form.
    for x in [
        14710.0, 2207.0, 502.8, 142.6, 46.47, 16.7, 6.356, 1.316, 0.3897,
    ] {
        let short = format_g(x, 15, true);
        assert_eq!(float_to_string(x), short, "{x}");
        assert_eq!(short.parse::<f64>().expect("parses"), x);
    }
}

#[test]
fn float_to_sexp_is_the_rule_as_an_atom() {
    assert_eq!(float_to_sexp(1.0).to_string(), "1");
    assert_eq!(float_to_sexp(0.1 + 0.2).to_string(), "0.30000000000000004");
    assert_eq!(float_to_sexp(-0.0).to_string(), "-0");
    assert_eq!(
        Sexp::list(vec![float_to_sexp(0.5), Sexp::atom("x")]).to_string(),
        "(0.5 x)"
    );
}
