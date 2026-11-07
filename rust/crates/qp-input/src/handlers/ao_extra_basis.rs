//! `Ao_extra_basis`: the auxiliary AO basis, edited like the AO basis.
//!
//! Port of `ocaml/Input_ao_extra_basis.ml`.  The group is `ao_extra_basis`
//! and every key is prefixed `ao_extra_`:
//!
//! | field | EZFIO key | element type |
//! |---|---|---|
//! | `ao_extra_basis` | `ao_extra_basis/ao_extra_basis` | `character*(256)` scalar |
//! | `ao_extra_num` | `ao_extra_basis/ao_extra_num` | integer scalar |
//! | `ao_extra_prim_num` | `ao_extra_basis/ao_extra_prim_num` | integer array, length `ao_extra_num` |
//! | `ao_extra_prim_num_max` | derived | integer scalar |
//! | `ao_extra_nucl` | `ao_extra_basis/ao_extra_nucl` | integer array, length `ao_extra_num` |
//! | `ao_extra_power` | `ao_extra_basis/ao_extra_power` | integer array, `(ao_extra_num,3)` |
//! | `ao_extra_coef` | `ao_extra_basis/ao_extra_coef` | real array, `(ao_extra_num,ao_extra_prim_num_max)` |
//! | `ao_extra_expo` | `ao_extra_basis/ao_extra_expo` | real array, same shape |
//! | `ao_extra_cartesian` | `ao_extra_basis/ao_extra_cartesian` | logical scalar, default `false` |
//! | `ao_extra_normalized` | `ao_extra_basis/ao_extra_normalized` | logical scalar, default `true` |
//! | `primitives_normalized_extra` | `ao_extra_basis/primitives_normalized_extra` | logical scalar, default `true` |
//!
//! # The generated TOML is wrong; this file is not
//!
//! `scripts/ezfio_interface` derives the module config from `EZFIO.cfg`, and
//! the `ao_extra_basis` config it produces points at `ezfio_dir =
//! ["ao_basis"]` with `source_cfg = "src/ao_basis/EZFIO.cfg"`, listing only
//! the `ao_*` fields.  The reference handler does **not** read `ao_basis`:
//! it uses its own group and its own `ao_extra_*` keys (see the table).  The
//! group and key names below are therefore hard-coded from the OCaml source,
//! not taken from the config.  See `NOTES.md`.
//!
//! # Not ported
//!
//! `of_long_basis` builds a record from a `Long_basis.t`; `qp_edit` never
//! calls it (no arm of `get`/`set`/`reorder` uses it, and it is not exported
//! to another handler), so it is documented here and left out.  Its
//! `AO_extra_basis` wrapper would also need `AO_coef.of_float` on values that
//! the digest path already validates.
//!
//! # The absent-section back-fill
//!
//! `read`'s `try ... with _ -> None` arm writes `MD5("None")` and, because
//! OCaml builds records right-to-left, the three logicals have already been
//! materialised when `read_ao_extra_basis ()` raises `Failure "No basis"`.
//! On a tree that never had the section the reference therefore leaves
//! `ao_extra_cartesian`, `ao_extra_normalized`, `primitives_normalized_extra`
//! and `ao_extra_md5` behind (plus the group's `.version`).  This port
//! reproduces that set exactly — the values come from
//! `$QP_ROOT/data/ezfio_defaults/ao_extra_basis.ezfio_interface_default` at
//! run time, and the stamp is written only when it is absent, so a directory
//! that already has it is never modified (contract C2).
//!
//! Units: exponents are in bohr⁻², coefficients are dimensionless.

use anyhow::{Context, Result};
use qp_ezfio::EzFio;
use qp_types::basis;
use qp_types::gto::{AngmomXyz, Fmt};
use qp_types::md5::MD5_OF_NONE;

use super::ao_basis::{self, Shells};
use crate::engine::make_header;

/// The EZFIO group this handler owns.
const GROUP: &str = "ao_extra_basis";

/// The section of `$QP_ROOT/data/ezfio_defaults/` the defaulted fields are
/// read from (`Qpackage.get_ezfio_default "ao_extra_basis"` in
/// `Input_ao_extra_basis.ml`).  It equals [`GROUP`] today, but it is a
/// different concept: the section is the module id derived from `EZFIO.cfg`,
/// the group is the EZFIO directory.
const MODULE: &str = "ao_extra_basis";

/// One auxiliary basis, in the field order of the OCaml record.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AoExtraBasis {
    /// Basis-set name stored in `ao_extra_basis/ao_extra_basis`.
    pub(crate) ao_extra_basis: String,
    /// Number of Cartesian functions.
    pub(crate) ao_extra_num: i64,
    /// Number of primitives of each Cartesian function.
    pub(crate) ao_extra_prim_num: Vec<i64>,
    /// `max(ao_extra_prim_num)`, or 0 when the array is absent.
    pub(crate) ao_extra_prim_num_max: i64,
    /// Nucleus index of each Cartesian function (1-based).
    pub(crate) ao_extra_nucl: Vec<i64>,
    /// Cartesian monomial of each function.
    pub(crate) ao_extra_power: Vec<AngmomXyz>,
    /// Contraction coefficients, column major.
    pub(crate) ao_extra_coef: Vec<f64>,
    /// Exponents, same layout as `ao_extra_coef`.
    pub(crate) ao_extra_expo: Vec<f64>,
    /// Whether the basis uses Cartesian functions.
    pub(crate) ao_extra_cartesian: bool,
    /// Whether the contracted basis functions are normalised.
    pub(crate) ao_extra_normalized: bool,
    /// Whether the primitive functions are normalised.
    pub(crate) primitives_normalized_extra: bool,
}

impl Shells for AoExtraBasis {
    fn num(&self) -> i64 {
        self.ao_extra_num
    }
    fn prim_num(&self) -> &[i64] {
        &self.ao_extra_prim_num
    }
    fn prim_num_max(&self) -> i64 {
        self.ao_extra_prim_num_max
    }
    fn power(&self) -> &[AngmomXyz] {
        &self.ao_extra_power
    }
    fn nucl(&self) -> &[i64] {
        &self.ao_extra_nucl
    }
    fn coef(&self) -> &[f64] {
        &self.ao_extra_coef
    }
    fn expo(&self) -> &[f64] {
        &self.ao_extra_expo
    }
}

/// The record fields, without the `try ... with _ -> None` wrapper.
///
/// **Field order is load-bearing here.**  OCaml evaluates record fields
/// right-to-left, so `read_primitives_normalized_extra`, `read_ao_extra_normalized`
/// and `read_ao_extra_cartesian` run — and materialise their defaults —
/// *before* `read_ao_extra_basis ()` can raise `Failure "No basis"`.  On a
/// tree that never had this section the reference therefore leaves the three
/// logicals behind when the basis is absent; a Rust struct literal (evaluated
/// in source order) would not.  The three `let`s below restore that order.
///
/// The logicals take their value from `$QP_ROOT/data/ezfio_defaults/` at run
/// time (`get_default "<key>"`), falling back to the `EZFIO.cfg` constants
/// (`false`, `true`, `true`) only when no installed value exists.
fn read_record(ezfio: &EzFio) -> Result<AoExtraBasis> {
    let primitives_normalized_extra =
        ao_basis::read_logical(ezfio, MODULE, GROUP, "primitives_normalized_extra", true)?;
    let ao_extra_normalized =
        ao_basis::read_logical(ezfio, MODULE, GROUP, "ao_extra_normalized", true)?;
    let ao_extra_cartesian =
        ao_basis::read_logical(ezfio, MODULE, GROUP, "ao_extra_cartesian", false)?;
    Ok(AoExtraBasis {
        ao_extra_basis: ao_basis::read_name(ezfio, GROUP, "ao_extra_basis")?,
        ao_extra_num: ao_basis::read_num(ezfio, GROUP, "ao_extra_num")?,
        ao_extra_prim_num: ao_basis::read_prim_num(ezfio, GROUP, "ao_extra_prim_num")?,
        ao_extra_prim_num_max: ao_basis::read_prim_num_max(ezfio, GROUP, "ao_extra_prim_num")?,
        ao_extra_nucl: ao_basis::read_nucl(ezfio, GROUP, "ao_extra_nucl")?,
        ao_extra_power: ao_basis::read_power(ezfio, GROUP, "ao_extra_power")?,
        ao_extra_coef: ao_basis::read_coef(ezfio, GROUP, "ao_extra_coef")?,
        ao_extra_expo: ao_basis::read_expo(ezfio, GROUP, "ao_extra_expo")?,
        ao_extra_cartesian,
        ao_extra_normalized,
        primitives_normalized_extra,
    })
}

/// `read`: `None` on any failure of the record or of the digest.
///
/// The OCaml body is one `try ... with _ -> None`; on failure it writes
/// `MD5("None")` to `ao_extra_md5`.  Decision D4 forbids rewriting the stamp
/// of a calculation that already has one, but contract C3 requires the key on
/// a tree that never had the section (the reference's writer also creates the
/// group and its `.version` there).  The compromise is the guarded write in
/// [`stamp_md5_of_none_when_absent`]: absent key → create it, present key →
/// leave it byte-identical.  The digest is still computed, because the
/// reference's `try` covers it: a record whose digest cannot be built (a
/// short coefficient matrix, a missing `x` component) reads as `None`.
pub(crate) fn read(ezfio: &EzFio) -> Result<Option<AoExtraBasis>> {
    let record = match read_record(ezfio) {
        Ok(record) => record,
        Err(_) => {
            stamp_md5_of_none_when_absent(ezfio)?;
            return Ok(None);
        }
    };
    match ao_basis::to_md5(&record) {
        Ok(_) => Ok(Some(record)),
        Err(_) => {
            stamp_md5_of_none_when_absent(ezfio)?;
            Ok(None)
        }
    }
}

/// The `with _ ->` arm's write: `MD5("None")`, only when the key is absent.
///
/// Writing an absent key creates the group and its `.version`, exactly like
/// the reference; rewriting an existing stamp would violate contract C2.
fn stamp_md5_of_none_when_absent(ezfio: &EzFio) -> Result<()> {
    if !ezfio.has(GROUP, "ao_extra_md5")? {
        ezfio.set::<String>(GROUP, "ao_extra_md5", &MD5_OF_NONE.to_owned())?;
    }
    Ok(())
}

/// `to_basis`: the collapsed basis, one entry per shell.
pub(crate) fn to_basis(b: &AoExtraBasis) -> Result<basis::Basis> {
    ao_basis::to_basis(b)
}

/// `to_md5`: `Basis.to_md5` of the collapsed basis.
pub(crate) fn to_md5(b: &AoExtraBasis) -> Result<String> {
    ao_basis::to_md5(b)
}

/// `ordering`: see [`ao_basis::ordering`]; the algorithm is identical.
pub(crate) fn ordering(b: &AoExtraBasis) -> Result<Vec<i64>> {
    ao_basis::ordering(b)
}

/// `reorder`: see [`ao_basis::reorder`]; the algorithm is identical.
pub(crate) fn reorder(b: &AoExtraBasis) -> Result<AoExtraBasis> {
    let parts = ao_basis::reorder_parts(b)?;
    Ok(AoExtraBasis {
        ao_extra_prim_num: parts.prim_num,
        ao_extra_nucl: parts.nucl,
        ao_extra_power: parts.power,
        ao_extra_coef: parts.coef,
        ao_extra_expo: parts.expo,
        ..b.clone()
    })
}

/// `write`: the eleven fields, in record order, MD5 first.
///
/// The three logicals are written before the matrices, between `power` and
/// `coef`, exactly as in the reference.
pub(crate) fn write(ezfio: &EzFio, b: &AoExtraBasis) -> Result<()> {
    let num = ao_basis::usize_of(b.ao_extra_num)?;
    let prim_num_max = ao_basis::usize_of(b.ao_extra_prim_num_max)?;

    ezfio.set::<String>(GROUP, "ao_extra_md5", &to_md5(b)?)?;
    ezfio.set::<String>(GROUP, "ao_extra_basis", &b.ao_extra_basis)?;

    let prim_num = super::ezfio_array_of_list(vec![num], b.ao_extra_prim_num.clone())?;
    ezfio.set_array::<i64>(GROUP, "ao_extra_prim_num", &prim_num)?;

    let nucl = super::ezfio_array_of_list(vec![num], b.ao_extra_nucl.clone())?;
    ezfio.set_array::<i64>(GROUP, "ao_extra_nucl", &nucl)?;

    let mut power = Vec::with_capacity(3 * b.ao_extra_power.len());
    power.extend(b.ao_extra_power.iter().map(|xyz| i64::from(xyz.x)));
    power.extend(b.ao_extra_power.iter().map(|xyz| i64::from(xyz.y)));
    power.extend(b.ao_extra_power.iter().map(|xyz| i64::from(xyz.z)));
    let power = super::ezfio_array_of_list(vec![num, 3], power)?;
    ezfio.set_array::<i64>(GROUP, "ao_extra_power", &power)?;

    ezfio.set::<bool>(GROUP, "ao_extra_cartesian", &b.ao_extra_cartesian)?;
    ezfio.set::<bool>(GROUP, "ao_extra_normalized", &b.ao_extra_normalized)?;
    ezfio.set::<bool>(
        GROUP,
        "primitives_normalized_extra",
        &b.primitives_normalized_extra,
    )?;

    let coef = super::ezfio_array_of_list(vec![num, prim_num_max], b.ao_extra_coef.clone())?;
    ezfio.set_array::<f64>(GROUP, "ao_extra_coef", &coef)?;

    let expo = super::ezfio_array_of_list(vec![num, prim_num_max], b.ao_extra_expo.clone())?;
    ezfio.set_array::<f64>(GROUP, "ao_extra_expo", &expo)?;
    Ok(())
}

/// `to_rst`: the editable block, without the section header.
///
/// Same body as [`ao_basis::to_rst`], with the two normalisation questions
/// inserted between the Cartesian flag and the basis set.  The order of the
/// questions is the order of the `Printf.sprintf` arguments:
/// `primitives_normalized_extra` first, then `ao_extra_normalized`.
pub(crate) fn to_rst(b: &AoExtraBasis) -> Result<String> {
    let short_basis = to_basis(b)?;
    let basis_text = basis::to_string(Fmt::Gamess, &short_basis, None);
    let basis_text = basis_text
        .split('\n')
        .map(|line| format!("  {line}"))
        .collect::<Vec<_>>()
        .join("\n");

    let mut print_sym = String::new();
    for (i, xyz) in b.ao_extra_power.iter().enumerate() {
        let nucl = b
            .ao_extra_nucl
            .get(i)
            .context("Invalid_argument: ao_extra_nucl index out of bounds")?;
        // `AngmomXyz`'s `Display` ignores the formatter width (it writes
        // through `write_str`), so the symbol string is padded, not the value.
        let xyz = xyz.to_symbol_monomial();
        print_sym.push_str(&format!(" {:>5}  {:>6}     {:<9}\n", i + 1, nucl, xyz));
    }

    Ok(format!(
        "\nName of the AO basis ::\n\n  ao_extra_basis = {}\n\n\
         Cartesian coordinates (6d,10f,...) ::\n\n  ao_extra_cartesian = {}\n\n\
         Use normalized primitive functions ::\n\n  primitives_normalized_extra = {}\n\n\
         Use normalized basis functions ::\n\n  ao_extra_normalized = {}\n\n\
         Basis set (read-only) ::\n\n{}\n\n\n\
         ======= ========= ===================\n Basis   Nucleus   Cartesian function\n\
         ======= ========= ===================\n{}\n\
         ======= ========= ===================\n\n",
        b.ao_extra_basis,
        b.ao_extra_cartesian,
        b.primitives_normalized_extra,
        b.ao_extra_normalized,
        basis_text,
        print_sym
    ))
}

/// `get` for this section: `None` renders as `""`.
pub(crate) fn get(ezfio: &EzFio, title: &str) -> Result<String> {
    match read(ezfio)? {
        Some(record) => Ok(make_header(title) + &to_rst(&record)?),
        None => Ok(String::new()),
    }
}

/// A minimal, already ordered single-shell basis used by the tests.
#[cfg(test)]
fn one_s_shell() -> AoExtraBasis {
    AoExtraBasis {
        ao_extra_basis: "test".to_owned(),
        ao_extra_num: 1,
        ao_extra_prim_num: vec![1],
        ao_extra_prim_num_max: 1,
        ao_extra_nucl: vec![1],
        ao_extra_power: vec![AngmomXyz::new(0, 0, 0)],
        ao_extra_coef: vec![1.0],
        ao_extra_expo: vec![2.0],
        ao_extra_cartesian: false,
        ao_extra_normalized: true,
        primitives_normalized_extra: true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::handlers::test_support::hcn_ezfio;
    use qp_types::md5::MD5_OF_NONE;

    /// The `hcn` fixture has only the three logicals and a stamp: no name,
    /// so the section is absent.  The reference rewrites the stamp with
    /// `MD5("None")` on that path; the port must leave it alone (D4).
    #[test]
    fn read_returns_none_and_does_not_restamp() {
        let (_dir, ezfio) = hcn_ezfio();
        assert!(read(&ezfio).expect("read").is_none());
        assert_eq!(
            ezfio.get::<String>(GROUP, "ao_extra_md5").expect("md5"),
            MD5_OF_NONE
        );
        assert!(!ezfio.has(GROUP, "ao_extra_basis").expect("has"));
        assert!(!ezfio.has(GROUP, "ao_extra_num").expect("has"));
    }

    /// On a tree that never had the section, `read` leaves exactly the files
    /// the reference leaves: the three logicals (materialised before the
    /// missing name fails the record), the `MD5("None")` stamp, and the
    /// group's `.version` created by the first write.  Contract C3.
    #[test]
    fn read_on_an_absent_section_creates_the_reference_file_set() {
        use std::fs;

        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        assert!(read(&ezfio).expect("read").is_none());

        let mut names: Vec<String> = fs::read_dir(ezfio.root().join(GROUP))
            .expect("group directory")
            .map(|entry| {
                entry
                    .expect("entry")
                    .file_name()
                    .to_string_lossy()
                    .into_owned()
            })
            .collect();
        names.sort();
        assert_eq!(
            names,
            [
                ".version",
                "ao_extra_cartesian",
                "ao_extra_md5",
                "ao_extra_normalized",
                "primitives_normalized_extra",
            ]
        );
        assert_eq!(
            ezfio.get::<String>(GROUP, "ao_extra_md5").expect("md5"),
            MD5_OF_NONE
        );
        assert!(
            !ezfio
                .get::<bool>(GROUP, "ao_extra_cartesian")
                .expect("cartesian")
        );
        assert!(
            ezfio
                .get::<bool>(GROUP, "ao_extra_normalized")
                .expect("normalized")
        );
        assert!(
            ezfio
                .get::<bool>(GROUP, "primitives_normalized_extra")
                .expect("primitives")
        );
    }

    #[test]
    fn read_materialises_the_three_logicals_with_their_defaults() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        ezfio
            .set::<String>(GROUP, "ao_extra_basis", &"test".to_owned())
            .expect("name");
        ezfio.set::<i64>(GROUP, "ao_extra_num", &1).expect("num");
        // `AO_prim_number.of_int 0` would fail the whole read before the
        // logicals are reached, so the primitive array has to be there.
        ezfio
            .set_array::<i64>(
                GROUP,
                "ao_extra_prim_num",
                &super::super::ezfio_array_of_list(vec![1], vec![1]).expect("array"),
            )
            .expect("prim_num");
        // `read` still fails (no power, so the digest cannot be built), but
        // the logicals are materialised on the way through the fields.
        assert!(read(&ezfio).expect("read").is_none());
        assert!(
            !ezfio
                .get::<bool>(GROUP, "ao_extra_cartesian")
                .expect("cartesian")
        );
        assert!(
            ezfio
                .get::<bool>(GROUP, "ao_extra_normalized")
                .expect("normalized")
        );
        assert!(
            ezfio
                .get::<bool>(GROUP, "primitives_normalized_extra")
                .expect("primitives")
        );
    }

    /// The extra handler uses its own group: a value stored under `ao_basis`
    /// must not be picked up.
    #[test]
    fn read_ignores_the_ao_basis_group() {
        let (_dir, ezfio) = hcn_ezfio();
        assert!(ezfio.has("ao_basis", "ao_basis").expect("has"));
        let record = read(&ezfio).expect("read");
        assert!(record.is_none(), "read fell back to the ao_basis group");
    }

    /// `write` does **not** stamp `ao_extra_num`: like `ao_num`, the
    /// reference only reads it (the Fortran side derives it from
    /// `ao_extra_prim_num`).  On a directory that already has the key the
    /// round trip is complete.
    #[test]
    fn write_round_trips_through_read() {
        let (_dir, ezfio) = hcn_ezfio();
        let record = one_s_shell();
        ezfio.set::<i64>(GROUP, "ao_extra_num", &1).expect("num");
        write(&ezfio, &record).expect("write");
        let again = read(&ezfio).expect("read").expect("some");
        assert_eq!(record, again);
        assert_eq!(
            ezfio.get::<String>(GROUP, "ao_extra_md5").expect("md5"),
            to_md5(&record).expect("digest")
        );
    }

    #[test]
    fn to_rst_has_the_two_normalisation_questions() {
        let record = one_s_shell();
        let rendered = to_rst(&record).expect("to_rst");
        assert!(rendered.starts_with(
            "\nName of the AO basis ::\n\n  ao_extra_basis = test\n\n\
             Cartesian coordinates (6d,10f,...) ::\n\n  ao_extra_cartesian = false\n\n\
             Use normalized primitive functions ::\n\n  primitives_normalized_extra = true\n\n\
             Use normalized basis functions ::\n\n  ao_extra_normalized = true\n\n"
        ));
        assert!(rendered.ends_with("======= ========= ===================\n\n"));
        assert!(rendered.contains("     1       1     s        \n"));
    }

    #[test]
    fn ordering_and_reorder_fix_a_misordered_shell() {
        // A P shell stored as z, x, y, with the shared contraction of a
        // real shell, is reordered to x, y, z.
        let record = AoExtraBasis {
            ao_extra_power: vec![
                AngmomXyz::new(0, 0, 1),
                AngmomXyz::new(1, 0, 0),
                AngmomXyz::new(0, 1, 0),
            ],
            ao_extra_num: 3,
            ao_extra_prim_num: vec![1, 1, 1],
            ao_extra_nucl: vec![1, 1, 1],
            ao_extra_coef: vec![1.0, 1.0, 1.0],
            ao_extra_expo: vec![2.0, 2.0, 2.0],
            ..one_s_shell()
        };
        assert_eq!(ordering(&record).expect("ordering"), vec![1, 2, 0]);
        let reordered = reorder(&record).expect("reorder");
        assert_eq!(
            reordered.ao_extra_power,
            vec![
                AngmomXyz::new(1, 0, 0),
                AngmomXyz::new(0, 1, 0),
                AngmomXyz::new(0, 0, 1)
            ]
        );
        assert_eq!(ordering(&reordered).expect("ordering"), vec![0, 1, 2]);
    }

    /// `do_work` ends on `assert false` when the three lists have different
    /// lengths.  A *shorter* `power` fails earlier, at `ao_power.(i)`.
    #[test]
    fn to_long_basis_reports_the_length_assertion() {
        let mut record = one_s_shell();
        record.ao_extra_power.push(AngmomXyz::new(0, 0, 0));
        let error = ao_basis::to_long_basis(&record).expect_err("mismatch");
        assert!(crate::ocaml_compat::is_assert_failure(&error), "{error}");

        // Shorter: the OCaml `Array.init` index error, not an assertion.
        let mut short = one_s_shell();
        short.ao_extra_num = 2;
        let error = ao_basis::to_long_basis(&short).expect_err("short");
        assert_eq!(
            error.to_string(),
            "Invalid_argument: ao_power index out of bounds"
        );
    }
}
