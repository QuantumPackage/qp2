//! The `Qptypes` numeric wrappers the handlers need, reproduced exactly.
//!
//! `ocaml/Qptypes.ml` wraps every integer and float the tool layer exchanges
//! with EZFIO in a one-field module whose `of_int`/`of_float` validates the
//! value.  The wrappers are transparent — `to_int x = x` — but their *checks*
//! are observable: they print warnings to stdout (`Qptypes.warning` is
//! `print_string`), raise `Invalid_argument` with a pinned message, or trip an
//! `assert` (which `qp_edit` maps to exit code 2, not 1).
//!
//! The port keeps plain `i64`/`f64` in the handler records and routes every
//! value through one of the functions below, mirroring the call sites of the
//! OCaml `read_*`/`write_*` functions.  Nothing here touches EZFIO except the
//! `*_get_max` readers, which are exactly the OCaml `get_max` functions: they
//! look at the *current* stored maximum, so they must be called at the same
//! point in the sequence as the reference does.
//!
//! # Error classes
//!
//! * A violated `assert` returns [`crate::ocaml_compat::AssertFailure`]
//!   (driver: banner, exit 2).
//! * A violated `Invalid_argument` returns a plain `anyhow` message (driver:
//!   banner, exit 1).
//! * [`warning`] is not an error: it writes to stdout, without a newline,
//!   exactly like `Qptypes.warning`.
//!
//! Units: none.  The float checks are unit-free range checks; the `%f` in
//! their messages uses six decimals, as OCaml's `Printf.sprintf "%f"`.

use anyhow::Result;
use qp_ezfio::EzFio;
use qp_types::sexplib_compat::format_f;

use crate::ocaml_compat::assertion;

/// `Qpackage.bit_kind`: bits per Fortran integer of the reference build.
///
/// Hardcoded from `src/bitmask/bitmasks_module.f90`
/// (`integer, parameter :: bit_kind_size = 64` and the matching `bit_kind`
/// of 8 bytes); the OCaml build reads the same values from the module file.
/// `Qptypes.Bit_kind` accepts `1 | 2 | 4 | 8` at parse time.
pub(crate) const BIT_KIND: i64 = 8;

/// `Qpackage.bit_kind_size`: bits per machine word, 64.
///
/// Same source as [`BIT_KIND`]; re-exported so the handlers never spell the
/// literal themselves.
pub(crate) use qp_types::bitlist::BIT_KIND_SIZE;

/// `Qptypes.warning`: `print_string s`, on **stdout**, with no newline.
///
/// The OCaml definition is literally `let warning = print_string`.  Callers
/// pass the message without a terminator, as the reference does.
pub(crate) fn warning(message: &str) {
    print!("{message}");
}

// ---------------------------------------------------------------------
// Integer wrappers with an EZFIO-backed maximum
// ---------------------------------------------------------------------

/// `N_int_number.get_max`: `determinants/n_int`, or 128 when absent.
pub(crate) fn n_int_get_max(ezfio: &EzFio) -> Result<i64> {
    if ezfio.has("determinants", "n_int")? {
        Ok(ezfio.get::<i64>("determinants", "n_int")?)
    } else {
        Ok(128)
    }
}

/// `N_int_number.of_int`: `min = 1`, `max = `[`n_int_get_max`].
///
/// The `x > 128` warning fires *before* the maximum check, as in OCaml.
pub(crate) fn n_int_of_int(ezfio: &EzFio, x: i64) -> Result<i64> {
    n_int_of_int_max(x, n_int_get_max(ezfio)?)
}

/// `N_int_number.of_int ~max`.
///
/// `max == 1` disables the upper check, the OCaml `match max with 1 -> ()`
/// quirk: a value of exactly one is the only legal one, so the check cannot
/// fail anyway.
pub(crate) fn n_int_of_int_max(x: i64, max: i64) -> Result<i64> {
    if x < 1 {
        return Err(assertion(format!(
            "N_int_number.of_int: assertion (x >= 1) failed for x = {x}"
        )));
    }
    if x > 128 {
        warning("N_int > 128");
    }
    if max != 1 && x > max {
        anyhow::bail!("N_int_number: {x}");
    }
    Ok(x)
}

/// `Bitlist.n_int_of_mo_num`: `(mo_num - 1) / 64 + 1`, through
/// `N_int_number.of_int` with the *current* maximum.
///
/// The division is OCaml's truncating integer division, so `mo_num = 0`
/// yields 1 (Rust's `/` on `i64` truncates the same way).
pub(crate) fn n_int_of_mo_num(ezfio: &EzFio, mo_num: i64) -> Result<i64> {
    n_int_of_int(ezfio, (mo_num - 1) / BIT_KIND_SIZE as i64 + 1)
}

/// `Det_number.get_max`: `determinants/n_det`, or 50 billion when absent.
pub(crate) fn det_number_get_max(ezfio: &EzFio) -> Result<i64> {
    if ezfio.has("determinants", "n_det")? {
        Ok(ezfio.get::<i64>("determinants", "n_det")?)
    } else {
        Ok(50_000_000_000)
    }
}

/// `Det_number.of_int`: `min = 1`, `max = `[`det_number_get_max`].
pub(crate) fn det_number_of_int(ezfio: &EzFio, x: i64) -> Result<i64> {
    det_number_of_int_max(x, det_number_get_max(ezfio)?)
}

/// `Det_number.of_int ~max`.
pub(crate) fn det_number_of_int_max(x: i64, max: i64) -> Result<i64> {
    if x < 1 {
        return Err(assertion(format!(
            "Det_number.of_int: assertion (x >= 1) failed for x = {x}"
        )));
    }
    if x > 50_000_000_000 {
        warning("More than 50 billion determinants");
    }
    if max != 1 && x > max {
        anyhow::bail!("Det_number: {x}");
    }
    Ok(x)
}

/// `States_number.of_int`: asserts `x > 0` and warns above 1000.
///
/// There is no maximum: `States_number` has no `~max` argument.
pub(crate) fn states_of_int(x: i64) -> Result<i64> {
    if x < 1 {
        return Err(assertion(format!(
            "States_number.of_int: assertion (x > 0) failed for x = {x}"
        )));
    }
    if x > 1000 {
        warning("More than 1000 states");
    }
    Ok(x)
}

/// `MO_number.get_max`: `mo_basis/mo_num`, or 10 000 when absent.
pub(crate) fn mo_number_get_max(ezfio: &EzFio) -> Result<i64> {
    if ezfio.has("mo_basis", "mo_num")? {
        Ok(ezfio.get::<i64>("mo_basis", "mo_num")?)
    } else {
        Ok(10_000)
    }
}

/// `MO_number.of_int ~max`.
pub(crate) fn mo_number_of_int_max(x: i64, max: i64) -> Result<i64> {
    if x < 1 {
        return Err(assertion(format!(
            "MO_number.of_int: assertion (x >= 1) failed for x = {x}"
        )));
    }
    if x > 10_000 {
        warning("More than 10_000 MOs");
    }
    if max != 1 && x > max {
        anyhow::bail!("MO_number: {x}");
    }
    Ok(x)
}

/// `MO_number.of_int` with the default maximum.
pub(crate) fn mo_number_of_int(ezfio: &EzFio, x: i64) -> Result<i64> {
    mo_number_of_int_max(x, mo_number_get_max(ezfio)?)
}

/// `AO_number.get_max`: `ao_basis/ao_num`, or 10 000 when absent.
pub(crate) fn ao_number_get_max(ezfio: &EzFio) -> Result<i64> {
    if ezfio.has("ao_basis", "ao_num")? {
        Ok(ezfio.get::<i64>("ao_basis", "ao_num")?)
    } else {
        Ok(10_000)
    }
}

/// `AO_number.of_int ~max`.
pub(crate) fn ao_number_of_int_max(x: i64, max: i64) -> Result<i64> {
    if x < 1 {
        return Err(assertion(format!(
            "AO_number.of_int: assertion (x >= 1) failed for x = {x}"
        )));
    }
    if x > 10_000 {
        warning("More than 10_000 AOs");
    }
    if max != 1 && x > max {
        anyhow::bail!("AO_number: {x}");
    }
    Ok(x)
}

/// `AO_number.of_int` with the default maximum.
pub(crate) fn ao_number_of_int(ezfio: &EzFio, x: i64) -> Result<i64> {
    ao_number_of_int_max(x, ao_number_get_max(ezfio)?)
}

/// `Nucl_number.get_max`: `nuclei/nucl_num`, or 10 000 when absent.
pub(crate) fn nucl_number_get_max(ezfio: &EzFio) -> Result<i64> {
    if ezfio.has("nuclei", "nucl_num")? {
        Ok(ezfio.get::<i64>("nuclei", "nucl_num")?)
    } else {
        Ok(10_000)
    }
}

/// `Nucl_number.of_int ~max`.
pub(crate) fn nucl_number_of_int_max(x: i64, max: i64) -> Result<i64> {
    if x < 1 {
        return Err(assertion(format!(
            "Nucl_number.of_int: assertion (x >= 1) failed for x = {x}"
        )));
    }
    if x > 10_000 {
        warning("More than 10_000 nuclei");
    }
    if max != 1 && x > max {
        anyhow::bail!("Nucl_number: {x}");
    }
    Ok(x)
}

/// `AO_prim_number.of_int`: asserts `x > 0`; no warning, no maximum.
pub(crate) fn ao_prim_number_of_int(x: i64) -> Result<i64> {
    if x < 1 {
        return Err(assertion(format!(
            "AO_prim_number.of_int: assertion (x > 0) failed for x = {x}"
        )));
    }
    Ok(x)
}

/// `Elec_alpha_number.of_int`: asserts `x > 0`.
pub(crate) fn elec_alpha_of_int(x: i64) -> Result<i64> {
    if x < 1 {
        return Err(assertion(format!(
            "Elec_alpha_number.of_int: assertion (x > 0) failed for x = {x}"
        )));
    }
    Ok(x)
}

/// `Elec_beta_number.of_int`: asserts `x >= 0`.
pub(crate) fn elec_beta_of_int(x: i64) -> Result<i64> {
    if x < 0 {
        return Err(assertion(format!(
            "Elec_beta_number.of_int: assertion (x >= 0) failed for x = {x}"
        )));
    }
    Ok(x)
}

/// `Bit_kind.of_int`: accepts `1`, `2`, `4` and `8`.
///
/// The rejection is an `Invalid_argument` with a fixed message and no
/// mention of the value, exactly as in `Qptypes.ml`.
pub(crate) fn bit_kind_of_int(x: i64) -> Result<i64> {
    match x {
        1 | 2 | 4 | 8 => Ok(x),
        _ => anyhow::bail!("Bit_kind should be (1|2|4|8)."),
    }
}

/// `AO_basis_name.of_string`: asserts the name is not empty.
///
/// The reference is `assert (x <> "") ; x`, so an empty string is an
/// `Assert_failure` (exit 2), not a `Failure`.  When does that happen?  Both
/// `read_ao_basis` and `read_ao_extra_basis` call it on the stored scalar
/// after checking it is not the literal `"None"`, so only a zero-byte
/// `ao_basis` file reaches it.
pub(crate) fn ao_basis_name_of_string(name: &str) -> Result<String> {
    if name.is_empty() {
        return Err(assertion(
            "AO_basis_name.of_string: assertion (x <> \"\") failed",
        ));
    }
    Ok(name.to_string())
}

/// `MD5.of_string`: asserts 32 bytes, each below `'g'`.
///
/// The reference checks `String.length x = 32` then folds `x.[i] < 'g'` over
/// the bytes; it does **not** recompute the digest.  The message below names
/// the wrapper and the rejected value, since the OCaml `Assert_failure`
/// payload is a source position the port does not have.
pub(crate) fn md5_of_string(digest: &str) -> Result<String> {
    if digest.len() != 32 || !digest.bytes().all(|byte| byte < b'g') {
        return Err(assertion(format!(
            "MD5.of_string: assertion (32 characters below 'g') failed for x = {digest:?}"
        )));
    }
    Ok(digest.to_string())
}

// ---------------------------------------------------------------------
// Float wrappers
// ---------------------------------------------------------------------

/// `Det_coef.of_float`: `-1 <= x <= 1`.
///
/// The message is `"Det_coef : (-1. <= x <= 1.) : x=%f"` with six decimals.
pub(crate) fn det_coef_of_float(x: f64) -> Result<f64> {
    // `!(x >= -1. && x <= 1.)` in the OCaml, NaN included.
    if x.is_nan() || !(-1.0..=1.0).contains(&x) {
        anyhow::bail!("Det_coef : (-1. <= x <= 1.) : x={}", format_f(x, 6));
    }
    Ok(x)
}

/// `AO_expo.of_float`: `x >= 0`.
pub(crate) fn ao_expo_of_float(x: f64) -> Result<f64> {
    if x.is_nan() || x < 0.0 {
        anyhow::bail!("AO_expo : (x >= 0.) : x={}", format_f(x, 6));
    }
    Ok(x)
}

/// `Positive_float.of_float`: `x >= 0` (NaN is rejected, as `not (x >= 0.)`).
pub(crate) fn positive_float_of_float(x: f64) -> Result<f64> {
    if x.is_nan() || x < 0.0 {
        anyhow::bail!("Positive_float : (x >= 0.) : x={}", format_f(x, 6));
    }
    Ok(x)
}

/// `MO_occ.of_float`: clamp to `[0, 2]`.
///
/// The reference clamps rather than rejects; a NaN falls through all three
/// branches and is returned unchanged (`if x < 0. then 0. else if x > 2.
/// then 2. else x`).  `f64::clamp` compares the same way, so it propagates
/// NaN unchanged too.
pub(crate) fn mo_occ_of_float(x: f64) -> f64 {
    x.clamp(0.0, 2.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch() -> (tempfile::TempDir, EzFio) {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        (dir, ezfio)
    }

    #[test]
    fn maxima_fall_back_to_the_documented_defaults() {
        let (_dir, ezfio) = scratch();
        assert_eq!(n_int_get_max(&ezfio).expect("n_int"), 128);
        assert_eq!(det_number_get_max(&ezfio).expect("n_det"), 50_000_000_000);
        assert_eq!(mo_number_get_max(&ezfio).expect("mo_num"), 10_000);
        assert_eq!(ao_number_get_max(&ezfio).expect("ao_num"), 10_000);
        assert_eq!(nucl_number_get_max(&ezfio).expect("nucl_num"), 10_000);
    }

    #[test]
    fn maxima_read_the_stored_value() {
        let (_dir, ezfio) = scratch();
        ezfio.set::<i64>("determinants", "n_int", &2).expect("set");
        ezfio.set::<i64>("determinants", "n_det", &3).expect("set");
        ezfio.set::<i64>("mo_basis", "mo_num", &55).expect("set");
        ezfio.set::<i64>("ao_basis", "ao_num", &59).expect("set");
        ezfio.set::<i64>("nuclei", "nucl_num", &3).expect("set");
        assert_eq!(n_int_get_max(&ezfio).expect("n_int"), 2);
        assert_eq!(det_number_get_max(&ezfio).expect("n_det"), 3);
        assert_eq!(mo_number_get_max(&ezfio).expect("mo_num"), 55);
        assert_eq!(ao_number_get_max(&ezfio).expect("ao_num"), 59);
        assert_eq!(nucl_number_get_max(&ezfio).expect("nucl_num"), 3);
    }

    /// The `max` check is `x > max` with `max == 1` meaning "unchecked", and
    /// the message is the OCaml `Invalid_argument` payload.
    #[test]
    fn integer_wrappers_report_the_ocaml_messages() {
        let error = n_int_of_int_max(3, 2).expect_err("3 > 2");
        assert_eq!(error.to_string(), "N_int_number: 3");
        assert_eq!(
            det_number_of_int_max(2, 1).expect("max 1 disables the check"),
            2
        );
        assert_eq!(n_int_of_int_max(129, 200).expect("warned"), 129);
        assert_eq!(mo_number_of_int_max(5, 1).expect("unchecked"), 5);
        let error = ao_number_of_int_max(0, 10).expect_err("assert");
        assert!(error.to_string().contains("x >= 1"), "{error}");
        assert!(crate::ocaml_compat::is_assert_failure(&error));
    }

    #[test]
    fn bit_kind_accepts_only_the_four_documented_widths() {
        for good in [1, 2, 4, 8] {
            assert_eq!(bit_kind_of_int(good).expect("accepted"), good);
        }
        let error = bit_kind_of_int(3).expect_err("rejected");
        assert_eq!(error.to_string(), "Bit_kind should be (1|2|4|8).");
    }

    #[test]
    fn float_wrappers_use_the_ocaml_messages() {
        assert_eq!(det_coef_of_float(1.0).expect("bound"), 1.0);
        assert_eq!(det_coef_of_float(-1.0).expect("bound"), -1.0);
        let error = det_coef_of_float(1.5).expect_err("out of range");
        assert_eq!(
            error.to_string(),
            "Det_coef : (-1. <= x <= 1.) : x=1.500000"
        );
        let error = ao_expo_of_float(-0.5).expect_err("negative");
        assert_eq!(error.to_string(), "AO_expo : (x >= 0.) : x=-0.500000");
        let error = positive_float_of_float(f64::NAN).expect_err("NaN");
        assert_eq!(error.to_string(), "Positive_float : (x >= 0.) : x=nan");
    }

    #[test]
    fn mo_occupation_clamps_and_lets_nan_through() {
        assert_eq!(mo_occ_of_float(-1.0), 0.0);
        assert_eq!(mo_occ_of_float(3.0), 2.0);
        assert_eq!(mo_occ_of_float(1.5), 1.5);
        assert!(mo_occ_of_float(f64::NAN).is_nan());
    }

    /// `n_int_of_mo_num` is `(mo_num - 1) / 64 + 1`, with truncating
    /// division: the reference returns 1 for `mo_num = 0`.
    #[test]
    fn n_int_of_mo_num_uses_the_reference_division() {
        let (_dir, ezfio) = scratch();
        assert_eq!(n_int_of_mo_num(&ezfio, 0).expect("zero"), 1);
        assert_eq!(n_int_of_mo_num(&ezfio, 1).expect("one"), 1);
        assert_eq!(n_int_of_mo_num(&ezfio, 64).expect("64"), 1);
        assert_eq!(n_int_of_mo_num(&ezfio, 65).expect("65"), 2);
        assert_eq!(n_int_of_mo_num(&ezfio, 55).expect("55"), 1);
    }

    #[test]
    fn ao_basis_name_rejects_only_the_empty_string() {
        assert_eq!(
            ao_basis_name_of_string("aug-cc-pvdz").expect("name"),
            "aug-cc-pvdz"
        );
        let error = ao_basis_name_of_string("").expect_err("empty");
        assert!(crate::ocaml_compat::is_assert_failure(&error));
    }

    #[test]
    fn md5_of_string_checks_the_length_and_the_alphabet() {
        let digest = "c4997d5a5cac61bf6d88d46396d6c608";
        assert_eq!(md5_of_string(digest).expect("digest"), digest);
        // 31 characters: the length assertion fires.
        let error = md5_of_string(&digest[..31]).expect_err("short");
        assert!(crate::ocaml_compat::is_assert_failure(&error));
        // 32 characters, but `'g'` is not below `'g'`.
        let error = md5_of_string(&format!("g{}", &digest[1..])).expect_err("alphabet");
        assert!(crate::ocaml_compat::is_assert_failure(&error));
        // Uppercase hex is accepted: `'A' < 'g'`.
        assert!(md5_of_string(&digest.to_uppercase()).is_ok());
    }
}
