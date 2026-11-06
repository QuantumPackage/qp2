//! Cartesian coordinates in bohr (`ocaml/Point3d.ml`), with the unit
//! conversions of `ocaml/Units.ml`.
//!
//! # Units
//!
//! Every [`Point3d`] stores **bohr**, always. [`Units`] only selects the
//! conversion applied at the boundary ([`Point3d::of_tuple`],
//! [`Point3d::of_string`], [`Point3d::to_string`]). The
//! Angstrom → bohr conversion is a **division** by
//! `bohr_to_angstrom = 0.52917721067121` ([`units::BOHR_TO_ANGSTROM`],
//! `Units.ml:8`), as in the reference `Point3d.to_bohr` added by commit
//! `15e725aa` (`Point3d.ml:15-19`: `x /. Units.bohr_to_angstrom`).
//!
//! Dividing, rather than multiplying by the pre-computed reciprocal
//! `angstrom_to_bohr = 1. /. 0.52917721067121`, is what makes the
//! Angstrom round trip of [`Point3d::to_string`]/[`Point3d::of_string`]
//! a fixed point: `fl(fl(x*b)/b)` returns to `x`, while
//! `fl(fl(x*b)*(1/b))` drifts by one ulp at each round trip on some
//! inputs. That is the rationale of the reference commit, mirrored here.
//!
//! # Printing
//!
//! [`Point3d::to_string`] is `%.17g` per component (`Point3d.ml:52-61`,
//! post-fix): 17 significant digits is the shortest form that always
//! reads back as the same `double`, which the former `%16.8f` did not.
//! The conversion to the requested unit happens *before* printing, so a
//! Bohr coordinate is printed unscaled and an Angstrom one is multiplied
//! by `bohr_to_angstrom`.

use anyhow::{Context, Result, bail};

use crate::sexplib_compat::{Sexp, float_to_sexp, format_g};

/// Unit systems (`ocaml/Units.ml`).
pub mod units {
    /// `Units.bohr_to_angstrom`, the CODATA value hard-coded in
    /// `Units.ml:8`.
    ///
    /// This is the **divisor** of the Angstrom → bohr **coordinate**
    /// conversion: [`Point3d::of_tuple`] and [`Point3d::of_string`] compute
    /// `x / BOHR_TO_ANGSTROM` (reference commit `15e725aa`,
    /// `Point3d.to_bohr`). It is also the multiplier of the reverse
    /// bohr → Angstrom rendering in [`Point3d::to_string`].
    pub const BOHR_TO_ANGSTROM: f64 = 0.52917721067121;

    /// `Units.angstrom_to_bohr = 1. /. bohr_to_angstrom`
    /// (`Units.ml:7`): the correctly-rounded reciprocal of
    /// [`BOHR_TO_ANGSTROM`].
    ///
    /// Used by [`crate::element`] for the Angstrom → bohr **radius**
    /// conversions (covalent and van der Waals), which still multiply by
    /// it because commit `15e725aa` did not change `Element.ml`.
    /// [`Point3d`] coordinates no longer use it: they divide by
    /// `BOHR_TO_ANGSTROM` instead.
    pub const ANGSTROM_TO_BOHR: f64 = 1.0 / BOHR_TO_ANGSTROM;

    /// Which unit a caller is speaking in. `Default` is `Angstrom`, the
    /// default of every optional `?units` argument in the OCaml tool
    /// layer.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
    pub enum Units {
        /// Atomic units: lengths in bohr, no conversion.
        Bohr,
        /// Lengths in angstrom, converted to bohr by **dividing** by
        /// [`BOHR_TO_ANGSTROM`] (reference commit `15e725aa`).
        #[default]
        Angstrom,
    }
}

pub use units::Units;

/// A point in space, in **bohr**.
///
/// The fields are public because OCaml's `Point3d.t` is a record with
/// public fields and callers pattern-match on it; no invariant relates
/// them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3d {
    /// `x / a0`.
    pub x: f64,
    /// `y / a0`.
    pub y: f64,
    /// `z / a0`.
    pub z: f64,
}

/// `Point3d.to_bohr`'s per-coordinate conversion — added by the reference
/// commit `15e725aa` (`Point3d.ml:15-19`) and called from `of_tuple`
/// (`:21-24`) and `of_string` (`:27-36`): apply the unit conversion of
/// `units` to one coordinate.
///
/// Angstrom is a **division** by [`units::BOHR_TO_ANGSTROM`]
/// (`0.52917721067121`, `Units.ml:8`), exactly the reference's
/// `x /. Units.bohr_to_angstrom`. Dividing, rather than multiplying by
/// the pre-computed reciprocal `1. /. bohr_to_angstrom`, is what makes
/// `fl(fl(x*b)/b)` a fixed point of the [`Point3d::to_string`] /
/// [`Point3d::of_string`] round trip; `fl(fl(x*b)*(1/b))` drifts by one
/// ulp at each round trip on some inputs.
fn to_bohr(units: Units, x: f64) -> f64 {
    match units {
        Units::Bohr => x,
        Units::Angstrom => x / units::BOHR_TO_ANGSTROM,
    }
}

impl Point3d {
    /// `Point3d.of_tuple ~units (x, y, z)`: convert each coordinate from
    /// `units` to bohr.
    pub fn of_tuple(units: Units, (x, y, z): (f64, f64, f64)) -> Point3d {
        Point3d {
            x: to_bohr(units, x),
            y: to_bohr(units, y),
            z: to_bohr(units, z),
        }
    }

    /// `Point3d.of_string ~units s`: read the first three whitespace-
    /// separated tokens of `s` as `x y z` in `units`.
    ///
    /// The split is on the single character `' '` followed by a filter of
    /// empty tokens (`Point3d.ml:29-30`), so runs of spaces and leading
    /// spaces are fine but a tab is *not* a separator. Extra tokens after
    /// the third are ignored, exactly as in the reference (`l.(0..2)`).
    ///
    /// Errors (never panics) when a token does not parse or when fewer
    /// than three tokens are present. OCaml raises `Failure` /
    /// `Invalid_argument` there; the message context names the token.
    pub fn of_string(units: Units, s: &str) -> Result<Point3d> {
        let tokens: Vec<&str> = s.split(' ').filter(|t| !t.is_empty()).collect();
        if tokens.len() < 3 {
            bail!(
                "Point3d: expected 3 coordinates, got {} in {s:?}",
                tokens.len()
            );
        }
        let mut c = [0.0_f64; 3];
        for (i, slot) in c.iter_mut().enumerate() {
            *slot = tokens[i]
                .parse::<f64>()
                .with_context(|| format!("Point3d: bad coordinate {:?}", tokens[i]))?;
        }
        Ok(Point3d::of_tuple(units, (c[0], c[1], c[2])))
    }

    /// `Point3d.distance2`: the squared distance between two points, in
    /// bohr².
    ///
    /// The sum is `((dx*dx) + (dy*dy)) + (dz*dz)`, left to right, as in
    /// `Point3d.ml:44`; the order is part of the reproducibility contract.
    ///
    /// OCaml wraps the result in `Positive_float`, which raises on a
    /// negative or NaN value. A sum of squares cannot be negative, so this
    /// returns the plain `f64` (NaN only if an input is NaN, which
    /// `Positive_float` would also accept).
    pub fn distance2(p1: &Point3d, p2: &Point3d) -> f64 {
        let dx = p2.x - p1.x;
        let dy = p2.y - p1.y;
        let dz = p2.z - p1.z;
        dx * dx + dy * dy + dz * dz
    }

    /// `Point3d.distance`: `sqrt` of [`Point3d::distance2`], in bohr.
    pub fn distance(p1: &Point3d, p2: &Point3d) -> f64 {
        Point3d::distance2(p1, p2).sqrt()
    }

    /// `Point3d.to_string ~units p`: `"x y z"`, each component printed
    /// with `%.17g` *after* conversion to `units`.
    ///
    /// One space between components, no padding. See
    /// [`crate::sexplib_compat::format_g`] for the `%g` semantics
    /// (17 significant digits, trailing zeros removed, exponent form
    /// outside `1e-4 .. 1e17`).
    pub fn to_string(units: Units, p: &Point3d) -> String {
        let f = match units {
            Units::Bohr => 1.0,
            Units::Angstrom => units::BOHR_TO_ANGSTROM,
        };
        format!(
            "{} {} {}",
            format_g(p.x * f, 17, false),
            format_g(p.y * f, 17, false),
            format_g(p.z * f, 17, false)
        )
    }

    /// `to_string` with each component right-aligned in a fixed-width field
    /// of `width` characters (no padding below the content: a component
    /// longer than `width` is printed verbatim, never truncated).
    ///
    /// The component *content* is byte-identical to [`Point3d::to_string`]
    /// (the same `format_g(.., 17, false)` per component), so the
    /// `to_string`/`of_string` round-trip property is unaffected; only
    /// whitespace is added, and the parser is whitespace-tolerant.
    ///
    /// Used for the fixed-width Molecule block (w10, DEVIATIONS C4);
    /// [`Point3d::to_string`] keeps the reference layout for the parity call
    /// sites.
    pub fn to_string_width(units: Units, p: &Point3d, width: usize) -> String {
        let f = match units {
            Units::Bohr => 1.0,
            Units::Angstrom => units::BOHR_TO_ANGSTROM,
        };
        [p.x, p.y, p.z]
            .iter()
            .map(|v| format_g(*v * f, 17, false))
            .map(|s| format!("{:>width$}", s, width = width))
            .collect::<Vec<_>>()
            .join(" ")
    }

    /// The derived `sexp_of_t`: `((x 1)(y 2)(z 3))` in bohr.
    ///
    /// Component order is `x`, `y`, `z` because that is the field
    /// declaration order in `Point3d.ml:5-9`; `[@@deriving sexp]` emits
    /// fields in that order and the md5 contract depends on it.
    pub fn to_sexp(p: &Point3d) -> Sexp {
        Sexp::list(vec![
            Sexp::pair(Sexp::atom("x"), float_to_sexp(p.x)),
            Sexp::pair(Sexp::atom("y"), float_to_sexp(p.y)),
            Sexp::pair(Sexp::atom("z"), float_to_sexp(p.z)),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn angstrom_round_trip_is_a_fixed_point() {
        // x -> of_string(Angstrom) -> to_string(Angstrom) must be exact
        // for these values: with the division of `15e725aa`,
        // `fl(fl(x*b)/b)` is a fixed point of the round trip, so a second
        // trip returns the same bohr point (measured 2026-11-06).
        let p = Point3d::of_string(Units::Angstrom, "0.54386314 0.00000000 -3.78645152")
            .expect("fixture parses");
        assert_eq!(
            Point3d::to_string(Units::Angstrom, &p),
            "0.54386314000000002 0 -3.78645152"
        );
        // And a second round trip is identical, not one ulp off.
        let q = Point3d::of_string(Units::Angstrom, &Point3d::to_string(Units::Angstrom, &p))
            .expect("own output parses");
        assert_eq!(q, p);
    }

    #[test]
    fn to_string_width_pads_components_without_truncating() {
        // Zero, short and long components (bohr, no conversion): the
        // padded fields carry the `to_string` content, each with length
        // `max(20, content)`.
        let p = Point3d::of_tuple(Units::Bohr, (0.0, 1.5, -1234567.0));
        let plain = Point3d::to_string(Units::Bohr, &p);
        let components = plain.split(' ').collect::<Vec<_>>();
        assert_eq!(components.len(), 3);
        for component in &components {
            assert_eq!(format!("{:>20}", *component).len(), 20.max(component.len()));
        }
        assert_eq!(
            Point3d::to_string_width(Units::Bohr, &p, 20),
            format!(
                "{:>20} {:>20} {:>20}",
                components[0], components[1], components[2]
            )
        );
    }

    #[test]
    fn to_string_width_prints_long_components_verbatim() {
        // An exponent form whose `%.17g` string is longer than the width:
        // printed verbatim, never truncated.  The expected string is built
        // from `format_g` itself (the content source of `to_string`).
        let v = 1e-5 * 1.234567890123456;
        let cx = format_g(v, 17, false);
        let cy = format_g(0.0, 17, false);
        let cz = format_g(-v, 17, false);
        assert!(
            cx.len() > 20 && cz.len() > 20,
            "fixture must exceed the width"
        );
        let p = Point3d::of_tuple(Units::Bohr, (v, 0.0, -v));
        let s = Point3d::to_string_width(Units::Bohr, &p, 20);
        assert_eq!(s, format!("{:>20} {:>20} {:>20}", cx, cy, cz));
        for content in [cx.as_str(), cy.as_str(), cz.as_str()] {
            assert!(s.contains(content), "content {content:?} truncated");
        }
    }

    #[test]
    fn to_string_width_keeps_the_sign_inside_the_field() {
        // A negative component keeps its sign inside the right-aligned
        // field: the padding sits to the left of the minus.
        let p = Point3d::of_tuple(Units::Bohr, (-2.5, 0.0, -7.0));
        let s = Point3d::to_string_width(Units::Bohr, &p, 20);
        let first = &s[..20];
        assert_eq!(first, format!("{:>20}", "-2.5"));
        assert!(first.starts_with(' ') && first.ends_with("-2.5"));
        assert_eq!(&s[21..], format!("{:>20} {:>20}", "0", "-7"));
    }

    #[test]
    fn to_string_width_round_trips_through_of_string() {
        // The parser ignores the padding, so the round trip of the Angstrom
        // fixture value is unaffected by the fixed width.
        let p = Point3d::of_string(Units::Angstrom, "0.54386314 0.00000000 -3.78645152")
            .expect("fixture parses");
        let s = Point3d::to_string_width(Units::Angstrom, &p, 20);
        assert_eq!(
            Point3d::of_string(Units::Angstrom, &s).expect("padded form parses"),
            p
        );
    }

    #[test]
    fn bohr_conversion_divides_by_bohr_to_angstrom() {
        let p = Point3d::of_tuple(Units::Bohr, (1.0, 2.0, 3.0));
        assert_eq!(
            p,
            Point3d {
                x: 1.0,
                y: 2.0,
                z: 3.0
            }
        );
        assert_eq!(Point3d::to_string(Units::Bohr, &p), "1 2 3");
        // 0.52917721067121 bohr is exactly 1 angstrom.
        let q = Point3d::of_tuple(Units::Angstrom, (1.0, 1.0, 1.0));
        assert_eq!(q.x, 1.0 / units::BOHR_TO_ANGSTROM);
        assert_eq!(Point3d::to_string(Units::Angstrom, &q), "1 1 1");
        // The c2h2 corpus value: `fl(x / b) != fl(x * (1 / b))` here
        // (they differ by exactly one ulp), so this pins the divide form
        // the reference uses since `15e725aa`.
        let h = 1.6644 / units::BOHR_TO_ANGSTROM;
        assert_eq!(Point3d::of_tuple(Units::Angstrom, (1.6644, 0.0, 0.0)).x, h);
        assert_ne!(h, 1.6644 * units::ANGSTROM_TO_BOHR);
    }

    #[test]
    fn distance_is_the_3_4_5_triangle() {
        let a = Point3d::of_tuple(Units::Bohr, (0.0, 0.0, 0.0));
        let b = Point3d::of_tuple(Units::Bohr, (3.0, 4.0, 0.0));
        assert_eq!(Point3d::distance2(&a, &b), 25.0);
        assert_eq!(Point3d::distance(&a, &b), 5.0);
        assert_eq!(Point3d::distance(&a, &a), 0.0);
    }

    #[test]
    fn of_string_needs_three_tokens() {
        assert!(Point3d::of_string(Units::Bohr, "1 2").is_err());
        assert!(Point3d::of_string(Units::Bohr, "").is_err());
        // A tab is not a separator in the reference (`String_ext.split
        // ~on:' '`), so this is one unparsable token, not three.
        assert!(Point3d::of_string(Units::Bohr, "1\t2\t3").is_err());
        // Extra tokens are ignored.
        let p = Point3d::of_string(Units::Bohr, "1 2 3 4").expect("three tokens present");
        assert_eq!(
            p,
            Point3d {
                x: 1.0,
                y: 2.0,
                z: 3.0
            }
        );
    }

    #[test]
    fn sexp_of_a_point_uses_the_field_order() {
        let p = Point3d::of_tuple(Units::Bohr, (1.0, 2.0, 3.0));
        assert_eq!(Point3d::to_sexp(&p).to_string(), "((x 1)(y 2)(z 3))");
    }
}
