//! Nuclear charges and electron counts (`ocaml/Charge.ml`).
//!
//! `Charge.t` is a bare `float` in the reference: ghost atoms carry
//! charge 0 and pseudo-potentials can make the charge fractional, so the
//! type is not an integer. The derived S-expression is a bare float atom,
//! which is why [`Charge::to_sexp`] exists next to the `Display` impl of
//! [`Charge`]: the printer rules differ (`%.17g` for the tool layer's text
//! files, `%.15G`-with-fallback for the md5 contract).

use anyhow::{Context, Result};

use crate::sexplib_compat::{Sexp, float_to_sexp, format_g};

/// A charge, in units of the elementary charge.
///
/// The inner value is public because OCaml's `Charge.t` is transparent
/// (`type t = float`); there is no invariant to protect (a fractional or
/// negative charge is meaningful: ghosts and pseudo-potentials).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Charge(pub f64);

impl Charge {
    /// `Charge.of_float`.
    pub fn of_float(x: f64) -> Charge {
        Charge(x)
    }

    /// `Charge.of_int`.
    pub fn of_int(i: i32) -> Charge {
        Charge(f64::from(i))
    }

    /// `Charge.of_string`: `float_of_string s`.
    ///
    /// Errors instead of raising `Failure`; the context names the token.
    /// OCaml's `float_of_string` accepts `_` separators and `0x` hex
    /// floats, which Rust's parser does not — nothing in the tool layer
    /// produces those, and a file containing one is an input error.
    pub fn of_string(s: &str) -> Result<Charge> {
        s.parse::<f64>()
            .map(Charge)
            .with_context(|| format!("Charge: not a float: {s:?}"))
    }

    /// `Charge.to_float`.
    pub fn to_float(self) -> f64 {
        self.0
    }

    /// `Charge.to_int = int_of_float`: truncation toward zero.
    ///
    /// OCaml leaves out-of-range conversions unspecified; Rust's `as`
    /// saturates and maps NaN to 0. Only used for display (`Molecule.name`)
    /// and for counting electrons, where the value is a small integer.
    pub fn to_int(self) -> i32 {
        self.0 as i32
    }

    /// The derived `sexp_of_t` of `type t = float`: a bare float atom,
    /// printed by the sexplib float rule.
    pub fn to_sexp(self) -> Sexp {
        float_to_sexp(self.0)
    }
}

impl std::fmt::Display for Charge {
    /// `Charge.to_string`: `%.17g` (post-fix `Charge.ml:12-15`).
    ///
    /// `%.17g`, not `%f`, because the charge is a float and `qp_edit`
    /// writes it back through this function: 17 significant digits always
    /// round-trip.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format_g(self.0, 17, false))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_string_is_seventeen_significant_digits() {
        assert_eq!(Charge::of_int(6).to_string(), "6");
        assert_eq!(Charge::of_float(0.1).to_string(), "0.10000000000000001");
        assert_eq!(Charge::of_int(-1).to_string(), "-1");
        assert_eq!(
            Charge::of_float(1.0 / 3.0).to_string(),
            "0.33333333333333331"
        );
    }

    #[test]
    fn to_int_truncates_toward_zero() {
        assert_eq!(Charge::of_float(6.9).to_int(), 6);
        assert_eq!(Charge::of_float(-6.9).to_int(), -6);
        assert_eq!(Charge::of_int(0).to_int(), 0);
    }

    #[test]
    fn sexp_uses_the_sexplib_float_rule_not_printf_17g() {
        // 0.1: %.15G round-trips, so the short form is used.
        assert_eq!(Charge::of_float(0.1).to_sexp().to_string(), "0.1");
        // 1/3 does not round-trip at 15 digits: fall back to %.17G.
        assert_eq!(
            Charge::of_float(1.0 / 3.0).to_sexp().to_string(),
            "0.33333333333333331"
        );
    }

    #[test]
    fn of_string_rejects_garbage() {
        assert!(Charge::of_string("+6.0").is_ok());
        assert!(Charge::of_string(" 6").is_err());
        assert!(Charge::of_string("six").is_err());
    }
}
