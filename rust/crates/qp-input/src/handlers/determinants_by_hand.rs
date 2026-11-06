//! `Determinants_by_hand`: the determinant list, the CI coefficients and the
//! state bookkeeping.
//!
//! Port of `ocaml/Input_determinants_by_hand.ml` (commit `15e725aa`).  This is
//! by far the largest handler: it reads and writes ten EZFIO keys, renders four
//! of them to RST, and parses the edited text back through a hand-built
//! S-expression because the determinant lines are bit strings, not numbers.
//!
//! # The RST section
//!
//! `to_rst` prints `expected_s2`, `n_det` and `state_average_weight` as
//! `key = value` lines, then one block per determinant of the *edited* list
//! (`n_det_qp_edit`):
//!
//! ```text
//!   <coefficients separated by tabs>
//!     <alpha bit string>
//!     <beta bit string>
//! ```
//!
//! `of_rst` rebuilds the record by concatenating two sources:
//!
//! * every `key = value` line of the **whole** edited text whose first
//!   character is a space (this keeps user edits to the three header values,
//!   including a duplicated line, which the record reader then rejects);
//! * the fields the OCaml reads from EZFIO instead of the RST: `mo_label`,
//!   `bit_kind`, `n_int`, `n_states`, `n_det_qp_edit` and the two generated
//!   arrays `psi_coef`/`psi_det`.
//!
//! The field order of the concatenation is fixed and is *not* the record
//! order; `ppx_sexp_conv` accepts any order.
//!
//! # Preserved quirks
//!
//! * `psi_coef` is state-major: state `j` occupies `[j*n_det, (j+1)*n_det)`.
//!   `to_rst` guards the column start (`ishift < len`) but not the row index
//!   (`i + ishift < len`); a short array therefore raises `Invalid_argument`
//!   in the reference.  The port returns an error with the same exit class.
//! * `read` catches **every** exception of the record fields and turns it into
//!   `None` — including a failed materialisation write and an `assert`
//!   violation.  Only `read_n_det_qp_edit` and `read_n_det`, evaluated before
//!   the `try`, can fail loudly.
//! * `write_n_states` overwrites `state_average_weight` with `1.0`s; the
//!   caller's weights are written afterwards and win.
//! * `write_n_det_qp_edit` clamps to the *stored* `n_det`, not to the record's.
//! * `of_rst` ignores the parsed `n_int` and `bit_kind`: it writes
//!   `N_int_number.get_max ()` and `Qpackage.bit_kind` (8).
//! * A coefficient line is found by the two-line pattern `"" :: c :: _`: the
//!   parser needs the blank line that `to_rst` emits before every coefficient
//!   line (and before every determinant).
//!
//! Units: none.  Coefficients are dimensionless; `expected_s2` is `S(S+1)`.

use anyhow::{Context, Result, bail};
use qp_ezfio::EzFio;
use qp_types::bitlist::Bitlist;
use qp_types::determinant::Determinant;

use super::qptypes;
use super::{ParseResult, apply_parsed, ezfio_array_of_list};
use crate::engine::make_header;
use crate::ocaml_compat::{assertion, string_of_float};
use crate::sexp::{self, Sexp, SexpError};

/// The EZFIO group this handler owns.
const GROUP: &str = "determinants";

/// The ten record fields, in declaration order.
const FIELDS: [&str; 10] = [
    "n_int",
    "bit_kind",
    "n_det",
    "n_det_qp_edit",
    "n_states",
    "expected_s2",
    "psi_coef",
    "psi_det",
    "state_average_weight",
    "mo_label",
];

/// The generated OCaml record (`Input_determinants_by_hand.t`).
///
/// The wrappers of `Qptypes` are transparent, so the port stores the plain
/// values; every validation is applied at the same place as the reference,
/// through [`super::qptypes`].  `psi_det` keeps the parsed determinant objects
/// rather than their words so that `write` and the sexp builder share one
/// representation.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct DeterminantsByHand {
    /// Number of 64-bit words per spin string (`N_int_number`).
    pub n_int: i64,
    /// Bytes per Fortran integer, always 8 (`Bit_kind`).
    pub bit_kind: i64,
    /// Number of determinants in the wave function (`Det_number`).
    pub n_det: i64,
    /// Number of determinants shown by the editor (`Det_number`).
    pub n_det_qp_edit: i64,
    /// Number of states (`States_number`).
    pub n_states: i64,
    /// Target `S(S+1)` (`Positive_float`).
    pub expected_s2: f64,
    /// State-major CI coefficients (`Det_coef array`).
    pub psi_coef: Vec<f64>,
    /// Determinants, in file order (`Determinant array`).
    pub psi_det: Vec<Determinant>,
    /// State-average weights (`Positive_float array`).
    pub state_average_weight: Vec<f64>,
    /// MO label spelling (`MO_label`).
    pub mo_label: String,
}

// ---------------------------------------------------------------------
// Materialising readers
// ---------------------------------------------------------------------

/// `MO_label.of_string`: case-insensitive and trimmed, with the reference's
/// rejection path.
///
/// On failure the reference prints the rejected string on **stdout**
/// (`print_endline`) before raising the `Failure`; the message is the exact
/// multi-line text of `ocaml/MO_label.ml`.
fn mo_label_of_string(s: &str) -> Result<String> {
    let canonical = match s.trim().to_ascii_lowercase().as_str() {
        "guess" => Some("Guess"),
        "canonical" => Some("Canonical"),
        "natural" => Some("Natural"),
        "localized" => Some("Localized"),
        "orthonormalized" => Some("Orthonormalized"),
        "mcscf" => Some("MCSCF"),
        "mom" => Some("MOM"),
        "none" => Some("None"),
        _ => None,
    };
    match canonical {
        Some(label) => Ok(label.to_string()),
        None => {
            println!("{s}");
            bail!(
                "MO_label should be one of:\nGuess | Orthonormalized | Canonical | Natural | Localized | MCSCF | MOM | None."
            )
        }
    }
}

/// `read_mo_label`: copy the MO basis label when the determinant label is
/// absent, then convert.
fn read_mo_label(ezfio: &EzFio) -> Result<String> {
    if !ezfio.has(GROUP, "mo_label")? && ezfio.has("mo_basis", "mo_label")? {
        let label = ezfio.get::<String>("mo_basis", "mo_label")?;
        ezfio.set::<String>(GROUP, "mo_label", &label)?;
    }
    let raw = ezfio.get::<String>(GROUP, "mo_label")?;
    mo_label_of_string(&raw)
}

/// `read_n_int`: derive `(mo_num - 1) / 64 + 1` when absent, then validate.
fn read_n_int(ezfio: &EzFio) -> Result<i64> {
    if !ezfio.has(GROUP, "n_int")? {
        let mo_num = ezfio.get::<i64>("mo_basis", "mo_num")?;
        let derived = qptypes::n_int_of_mo_num(ezfio, mo_num)?;
        ezfio.set::<i64>(GROUP, "n_int", &derived)?;
    }
    let value = ezfio.get::<i64>(GROUP, "n_int")?;
    qptypes::n_int_of_int(ezfio, value)
}

/// `read_bit_kind`: materialise `Qpackage.bit_kind` (8) when absent.
fn read_bit_kind(ezfio: &EzFio) -> Result<i64> {
    if !ezfio.has(GROUP, "bit_kind")? {
        ezfio.set::<i64>(GROUP, "bit_kind", &qptypes::BIT_KIND)?;
    }
    let value = ezfio.get::<i64>(GROUP, "bit_kind")?;
    qptypes::bit_kind_of_int(value)
}

/// `read_n_det`: materialise 1 when absent, then validate against itself.
fn read_n_det(ezfio: &EzFio) -> Result<i64> {
    if !ezfio.has(GROUP, "n_det")? {
        ezfio.set::<i64>(GROUP, "n_det", &1)?;
    }
    let value = ezfio.get::<i64>(GROUP, "n_det")?;
    qptypes::det_number_of_int(ezfio, value)
}

/// `read_n_det_qp_edit`: materialise from `n_det` when absent.
fn read_n_det_qp_edit(ezfio: &EzFio) -> Result<i64> {
    if !ezfio.has(GROUP, "n_det_qp_edit")? {
        let n_det = read_n_det(ezfio)?;
        ezfio.set::<i64>(GROUP, "n_det_qp_edit", &n_det)?;
    }
    let value = ezfio.get::<i64>(GROUP, "n_det_qp_edit")?;
    qptypes::det_number_of_int(ezfio, value)
}

/// `read_n_states`: materialise 1 when absent.
fn read_n_states(ezfio: &EzFio) -> Result<i64> {
    if !ezfio.has(GROUP, "n_states")? {
        ezfio.set::<i64>(GROUP, "n_states", &1)?;
    }
    let value = ezfio.get::<i64>(GROUP, "n_states")?;
    qptypes::states_of_int(value)
}

/// `read_expected_s2`: materialise `S(S+1)` from the electron counts.
///
/// The counts are read with the low-level getters, **without** the
/// `Elec_*_number.of_int` checks — the reference does not validate them here.
fn read_expected_s2(ezfio: &EzFio) -> Result<f64> {
    if !ezfio.has(GROUP, "expected_s2")? {
        let na = ezfio.get::<i64>("electrons", "elec_alpha_num")?;
        let nb = ezfio.get::<i64>("electrons", "elec_beta_num")?;
        let s = 0.5 * (na - nb) as f64;
        ezfio.set::<f64>(GROUP, "expected_s2", &(s * (s + 1.0)))?;
    }
    let value = ezfio.get::<f64>(GROUP, "expected_s2")?;
    qptypes::positive_float_of_float(value)
}

/// `read_psi_coef`: materialise the `(1, n_states)` unit vector when absent,
/// then read the flat coefficient array.
fn read_psi_coef(ezfio: &EzFio, read_only: bool) -> Result<Vec<f64>> {
    if !ezfio.has(GROUP, "psi_coef")? {
        let n_states = read_n_states(ezfio)?;
        let data: Vec<f64> = (0..n_states)
            .map(|i| if i == 0 { 1.0 } else { 0.0 })
            .collect();
        let array = ezfio_array_of_list(vec![1, usize_of(n_states)?], data)?;
        ezfio.set_array::<f64>(GROUP, "psi_coef", &array)?;
    }
    let key = if read_only {
        "psi_coef_qp_edit"
    } else {
        "psi_coef"
    };
    let array = ezfio.get_array::<f64>(GROUP, key)?;
    array
        .data()
        .iter()
        .map(|value| qptypes::det_coef_of_float(*value))
        .collect()
}

/// `read_psi_det`: materialise the Hartree–Fock determinant when absent, then
/// read and validate the `(n_int, 2, n_det)` word array.
fn read_psi_det(ezfio: &EzFio, read_only: bool) -> Result<Vec<Determinant>> {
    let n_int = read_n_int(ezfio)?;
    let n_alpha = ezfio.get::<i64>("electrons", "elec_alpha_num")?;
    let n_beta = ezfio.get::<i64>("electrons", "elec_beta_num")?;
    qptypes::elec_alpha_of_int(n_alpha)?;
    qptypes::elec_beta_of_int(n_beta)?;

    if !ezfio.has(GROUP, "psi_det")? {
        let mo_num = qptypes::mo_number_get_max(ezfio)?;
        let mut alpha_numbers = Vec::with_capacity(usize_of(n_alpha)?);
        for n in (1..=n_alpha).rev() {
            qptypes::mo_number_of_int_max(n, mo_num)?;
            alpha_numbers.push(usize_of(n)?);
        }
        let alpha = Bitlist::of_mo_number_list(usize_of(n_int)?, &alpha_numbers)?;
        let beta_numbers: Vec<usize> = (1..=n_beta).rev().map(usize_of).collect::<Result<_>>()?;
        let beta = Bitlist::of_mo_number_list(usize_of(n_int)?, &beta_numbers)?;
        let data: Vec<i64> = alpha
            .words()
            .iter()
            .chain(beta.words())
            .map(|word| *word as i64)
            .collect();
        let array = ezfio_array_of_list(vec![usize_of(n_int)?, 2, 1], data)?;
        ezfio.set_array::<i64>(GROUP, "psi_det", &array)?;
    }

    let n_int = usize_of(n_int)?;
    let key = if read_only {
        "psi_det_qp_edit"
    } else {
        "psi_det"
    };
    let array = ezfio.get_array::<i64>(GROUP, key)?;
    let dims = array.dims();
    // OCaml `assert`s; `read` maps the `Assert_failure` to `None`.
    if dims.len() != 3 || dims[0] != n_int {
        return Err(assertion(format!(
            "Input_determinants_by_hand.read_psi_det: assertion (n_int = dim.(0)) failed for \
             n_int = {n_int}, dim = {dims:?}"
        )));
    }
    if dims[1] != 2 {
        return Err(assertion(format!(
            "Input_determinants_by_hand.read_psi_det: assertion (dim.(1) = 2) failed for \
             dim = {dims:?}"
        )));
    }
    let expected_dets = usize_of(if read_only {
        read_n_det_qp_edit(ezfio)?
    } else {
        read_n_det(ezfio)?
    })?;
    if dims[2] != expected_dets {
        return Err(assertion(format!(
            "Input_determinants_by_hand.read_psi_det: assertion (dim.(2) = n_det) failed for \
             dim = {dims:?}, n_det = {expected_dets}"
        )));
    }
    let data = array.data();
    let mut determinants = Vec::with_capacity(dims[2]);
    for i in 0..dims[2] {
        let start = 2 * n_int * i;
        let words = &data[start..start + 2 * n_int];
        let words: Vec<u64> = words.iter().map(|word| *word as u64).collect();
        let determinant =
            Determinant::of_words(n_int, usize_of(n_alpha)?, usize_of(n_beta)?, &words)
                .with_context(|| format!("determinant {i}"))?;
        determinants.push(determinant);
    }
    Ok(determinants)
}

/// `read_state_average_weight`: materialise uniform weights when absent or
/// when the stored length differs from `n_states`.
fn read_state_average_weight(ezfio: &EzFio) -> Result<Vec<f64>> {
    let n_states = read_n_states(ezfio)?;
    if !ezfio.has(GROUP, "state_average_weight")? {
        write_state_average_weight(ezfio, &uniform_weights(n_states)?)?;
    }
    let array = ezfio.get_array::<f64>(GROUP, "state_average_weight")?;
    let result: Vec<f64> = array
        .data()
        .iter()
        .map(|value| qptypes::positive_float_of_float(*value))
        .collect::<Result<_>>()?;
    if result.len() == usize_of(n_states)? {
        Ok(result)
    } else {
        let data = uniform_weights(n_states)?;
        write_state_average_weight(ezfio, &data)?;
        Ok(data)
    }
}

/// `Array.init n_states (fun _ -> 1. / float_of_int n_states)`.
fn uniform_weights(n_states: i64) -> Result<Vec<f64>> {
    let n_states = n_states as f64;
    let mut weights = Vec::new();
    for _ in 0..usize_of(n_states as i64)? {
        weights.push(qptypes::positive_float_of_float(1.0 / n_states)?);
    }
    Ok(weights)
}

/// `read ?full ()`: `None` when there are no MOs, or when any record field
/// fails.
///
/// The two determinant counters are read **before** the guard and before the
/// `try`: a `n_det_qp_edit` above `n_det` is an `Invalid_argument` that escapes
/// `read` instead of turning into `None`, exactly as in the reference.
pub(crate) fn read(ezfio: &EzFio, full: bool) -> Result<Option<DeterminantsByHand>> {
    let n_det_qp_edit = read_n_det_qp_edit(ezfio)?;
    let n_det = read_n_det(ezfio)?;
    let read_only = if full { false } else { n_det_qp_edit != n_det };

    if !ezfio.has("mo_basis", "mo_num")? {
        // No molecular orbitals, so no determinants.
        return Ok(None);
    }
    match read_record(ezfio, read_only) {
        Ok(record) => Ok(Some(record)),
        Err(_) => Ok(None),
    }
}

/// The record fields, in declaration order.
///
/// OCaml evaluates record fields in an unspecified order; the port uses the
/// declaration order and the resulting set of materialised files is the same
/// (each writer only creates its own key).
fn read_record(ezfio: &EzFio, read_only: bool) -> Result<DeterminantsByHand> {
    Ok(DeterminantsByHand {
        n_int: read_n_int(ezfio)?,
        bit_kind: read_bit_kind(ezfio)?,
        n_det: read_n_det(ezfio)?,
        n_det_qp_edit: read_n_det_qp_edit(ezfio)?,
        expected_s2: read_expected_s2(ezfio)?,
        psi_coef: read_psi_coef(ezfio, read_only)?,
        psi_det: read_psi_det(ezfio, read_only)?,
        n_states: read_n_states(ezfio)?,
        state_average_weight: read_state_average_weight(ezfio)?,
        mo_label: read_mo_label(ezfio)?,
    })
}

// ---------------------------------------------------------------------
// Writers
// ---------------------------------------------------------------------

/// `write_n_states`: set `n_states` and overwrite the weights with `1.0`.
fn write_n_states(ezfio: &EzFio, n_states: i64) -> Result<()> {
    ezfio.set::<i64>(GROUP, "n_states", &n_states)?;
    if n_states < 0 {
        // OCaml `Array.make` raises `Invalid_argument "Array.make"`.
        bail!("Array.make");
    }
    let array = ezfio_array_of_list(vec![usize_of(n_states)?], vec![1.0; usize_of(n_states)?])?;
    ezfio.set_array::<f64>(GROUP, "state_average_weight", &array)?;
    Ok(())
}

/// `write_state_average_weight`: truncate to the **current** `n_states`.
fn write_state_average_weight(ezfio: &EzFio, weights: &[f64]) -> Result<()> {
    let n_states = read_n_states(ezfio)?;
    let array = ezfio_array_of_list(vec![usize_of(n_states)?], weights.to_vec())?;
    ezfio.set_array::<f64>(GROUP, "state_average_weight", &array)?;
    Ok(())
}

/// `write_n_det_qp_edit`: clamp to the stored `n_det`.
fn write_n_det_qp_edit(ezfio: &EzFio, n: i64) -> Result<()> {
    let n_det = read_n_det(ezfio)?;
    ezfio.set::<i64>(GROUP, "n_det_qp_edit", &n.min(n_det))?;
    Ok(())
}

/// `write_psi_coef`: write both the working and the editor copy.
fn write_psi_coef(ezfio: &EzFio, n_det: i64, n_states: i64, coefficients: &[f64]) -> Result<()> {
    let array = ezfio_array_of_list(
        vec![usize_of(n_det)?, usize_of(n_states)?],
        coefficients.to_vec(),
    )?;
    ezfio.set_array::<f64>(GROUP, "psi_coef", &array)?;
    ezfio.set_array::<f64>(GROUP, "psi_coef_qp_edit", &array)?;
    Ok(())
}

/// `write_psi_det`: flatten the determinants and write both copies.
fn write_psi_det(
    ezfio: &EzFio,
    n_int: i64,
    n_det: i64,
    determinants: &[Determinant],
) -> Result<()> {
    let mut data = Vec::new();
    for determinant in determinants {
        data.extend(determinant.words().iter().map(|word| *word as i64));
    }
    let array = ezfio_array_of_list(vec![usize_of(n_int)?, 2, usize_of(n_det)?], data)?;
    ezfio.set_array::<i64>(GROUP, "psi_det", &array)?;
    ezfio.set_array::<i64>(GROUP, "psi_det_qp_edit", &array)?;
    Ok(())
}

/// `write ?force`: the eight writes of the reference, in its order.
///
/// `psi_coef`/`psi_det` (and their `*_qp_edit` copies) are skipped when
/// `n_det > n_det_qp_edit` and `force` is false: the editor must not resurrect
/// determinants the user has just truncated away.
pub(crate) fn write(ezfio: &EzFio, record: &DeterminantsByHand, force: bool) -> Result<()> {
    ezfio.set::<i64>(GROUP, "n_int", &record.n_int)?;
    ezfio.set::<i64>(GROUP, "bit_kind", &record.bit_kind)?;
    ezfio.set::<i64>(GROUP, "n_det", &record.n_det)?;
    write_n_states(ezfio, record.n_states)?;
    ezfio.set::<f64>(GROUP, "expected_s2", &record.expected_s2)?;
    if force || record.n_det <= record.n_det_qp_edit {
        write_n_det_qp_edit(ezfio, record.n_det)?;
        write_psi_coef(ezfio, record.n_det, record.n_states, &record.psi_coef)?;
        write_psi_det(ezfio, record.n_int, record.n_det, &record.psi_det)?;
    }
    write_state_average_weight(ezfio, &record.state_average_weight)?;
    ezfio.set::<String>(GROUP, "mo_label", &record.mo_label)?;
    Ok(())
}

// ---------------------------------------------------------------------
// RST rendering
// ---------------------------------------------------------------------

/// `to_rst`: the rendered body, without the section header.
///
/// Fails (rather than panicking) when the record is shorter than the
/// `n_det_qp_edit`/`n_states` it announces, which is the reference's
/// `Invalid_argument` from `Array.init`/`Array.get`.
pub(crate) fn to_rst(ezfio: &EzFio, record: &DeterminantsByHand) -> Result<String> {
    let max = ezfio.get::<i64>("mo_basis", "mo_num")?;
    // The reference validates `mo_num` (warning above 10 000) and then never
    // uses it: `Determinant.to_string` ignores `~mo_num`.
    qptypes::mo_number_of_int_max(max, max)?;
    let n_states = read_n_states(ezfio)?;
    let n_det_qp_edit = record.n_det_qp_edit;
    let mut det_text = String::new();
    for i in 0..n_det_qp_edit {
        if i > 0 {
            det_text.push('\n');
        }
        let mut coefficients = String::new();
        for j in 0..n_states {
            if j > 0 {
                coefficients.push('\t');
            }
            let ishift = j * n_det_qp_edit;
            if ishift < record.psi_coef.len() as i64 {
                let index = usize_of(i + ishift).map_err(|_| {
                    anyhow::anyhow!("Input_determinants_by_hand.to_rst: negative index")
                })?;
                let value = record.psi_coef.get(index).ok_or_else(|| {
                    anyhow::anyhow!(
                        "Invalid_argument: index out of bounds in psi_coef (index {index}, length {})",
                        record.psi_coef.len()
                    )
                })?;
                coefficients.push_str(&string_of_float(*value));
            } else {
                coefficients.push_str("0.");
            }
        }
        let index = usize_of(i)?;
        let determinant = record.psi_det.get(index).ok_or_else(|| {
            anyhow::anyhow!(
                "Invalid_argument: index out of bounds in psi_det (index {index}, length {})",
                record.psi_det.len()
            )
        })?;
        let printed = determinant.to_string();
        let indented = printed
            .split('\n')
            .map(|line| format!("  {line}"))
            .collect::<Vec<_>>()
            .join("\n");
        det_text.push_str(&format!("  {coefficients}\n{indented}\n"));
    }

    let weights = record
        .state_average_weight
        .iter()
        .map(|weight| string_of_float(*weight))
        .collect::<Vec<_>>()
        .join("\t");
    Ok(format!(
        "\nForce the selected wave function to be an eigenfunction of S^2.\nIf true, input the expected value of S^2 ::\n\n  expected_s2 = {}\n\nNumber of determinants ::\n\n  n_det = {}\n\nState average weights ::\n\n  state_average_weight = ({weights})\n\nDeterminants ::\n\n{det_text}\n",
        string_of_float(record.expected_s2),
        record.n_det,
    ))
}

// ---------------------------------------------------------------------
// RST parsing
// ---------------------------------------------------------------------

/// `of_rst`: rebuild the record from the edited section.
///
/// Returns `Ok(None)` after printing the sexplib `Error:  (...)` when the
/// generated S-expression does not match the record; returns `Err` when an
/// EZFIO read or a `Qptypes.of_*` validation raises, which the driver answers
/// with `Info: Read error in <title>` and one retry.
pub(crate) fn of_rst(ezfio: &EzFio, body: &str) -> ParseResult<DeterminantsByHand> {
    let Some(idx) = body.find("\nDeterminants") else {
        // OCaml `assert false`.
        return Err(assertion(
            "Input_determinants_by_hand.of_rst: assertion false failed (no \"\\nDeterminants\")",
        ));
    };
    let determinants_text = &body[idx..];

    // The header is filtered from the *whole* body, not from the part before
    // "\nDeterminants".
    let mut header = String::new();
    for line in body.split('\n') {
        if line.is_empty() || !line.contains('=') || !line.starts_with(' ') {
            continue;
        }
        header.push('(');
        header.push_str(line.replace('=', " ").trim());
        header.push(')');
    }

    let mut determinant_lines: Vec<&str> = determinants_text.split('\n').map(str::trim).collect();
    if determinant_lines.is_empty() {
        bail!("Error in determinants");
    }
    determinant_lines.remove(0);

    let coefficients = read_coefficients(&determinant_lines)?;
    let psi_coef = coefficients_to_sexp(&coefficients)?;
    let psi_det = determinants_to_sexp(ezfio, &determinant_lines)?;

    let n_int = qptypes::n_int_get_max(ezfio)?;
    let source = format!(
        "{header}(mo_label {})(bit_kind {})(n_int {n_int})(n_states {})(psi_coef {psi_coef})(psi_det {psi_det})(n_det_qp_edit {})",
        read_mo_label(ezfio)?,
        qptypes::BIT_KIND,
        read_n_states(ezfio)?,
        read_n_det_qp_edit(ezfio)?,
    );
    // `Generic_input_of_rst.evaluate_sexp` builds `"(" ^ s ^ ")"` and uses
    // that *wrapped* string both for the parse and for `fail_msg`, so the
    // `##` snippet offsets are relative to it.
    let sexp_source = format!("({source})");

    match convert(&sexp_source) {
        Ok(record) => Ok(Some(record)),
        Err(error) => {
            sexp::fail_msg(&sexp_source, &error);
            Ok(None)
        }
    }
}

/// The `read_coefs` half of `of_rst`: one `Vec<f64>` per determinant.
///
/// The pattern needs the blank line before the coefficient line; a stray blank
/// line is skipped in pairs (`"" :: ""`) or reached through the catch-all.
fn read_coefficients(lines: &[&str]) -> Result<Vec<Vec<f64>>> {
    let mut coefficients = Vec::new();
    let mut index = 0;
    while index < lines.len() {
        if lines[index].is_empty() {
            if index + 1 < lines.len() && lines[index + 1].is_empty() {
                index += 2;
                continue;
            }
            if index + 1 < lines.len() {
                let values = lines[index + 1]
                    .split('\t')
                    .map(|text| {
                        let value: f64 = text
                            .parse()
                            .with_context(|| format!("Float.of_string {text:?}"))?;
                        qptypes::det_coef_of_float(value)
                    })
                    .collect::<Result<Vec<_>>>()?;
                coefficients.push(values);
                index += 2;
                continue;
            }
        }
        index += 1;
    }
    Ok(coefficients)
}

/// The `build_result` half of `of_rst`: state-major, space-separated floats in
/// the `((...))` that the `psi_coef` field expects.
fn coefficients_to_sexp(coefficients: &[Vec<f64>]) -> Result<String> {
    let first = coefficients
        .first()
        .ok_or_else(|| anyhow::anyhow!("Failure(\"hd\")"))?;
    let n_states = first.len();
    let mut values = Vec::with_capacity(n_states * coefficients.len());
    for state in 0..n_states {
        for row in coefficients {
            let value = row.get(state).ok_or_else(|| {
                anyhow::anyhow!(
                    "Invalid_argument: index out of bounds in a coefficient row (index {state}, \
                     length {})",
                    row.len()
                )
            })?;
            values.push(string_of_float(*value));
        }
    }
    Ok(format!("({})", values.join(" ")))
}

/// The `read_dets` half of `of_rst`: one `(w1 w2 ...)` per determinant.
fn determinants_to_sexp(ezfio: &EzFio, lines: &[&str]) -> Result<String> {
    let n_int = qptypes::n_int_of_int(
        ezfio,
        (qptypes::mo_number_get_max(ezfio)? - 1) / qptypes::BIT_KIND_SIZE as i64 + 1,
    )?;
    let n_alpha = ezfio.get::<i64>("electrons", "elec_alpha_num")?;
    let n_beta = ezfio.get::<i64>("electrons", "elec_beta_num")?;
    qptypes::elec_alpha_of_int(n_alpha)?;
    qptypes::elec_beta_of_int(n_beta)?;

    let mut determinants = String::new();
    let mut index = 0;
    while index < lines.len() {
        if lines[index].is_empty() && index + 3 < lines.len() {
            let alpha = Bitlist::of_string(lines[index + 2], '-', '+')?;
            let beta = Bitlist::of_string(lines[index + 3], '-', '+')?;
            let determinant = Determinant::of_bitlist_couple(
                usize_of(n_int)?,
                usize_of(n_alpha)?,
                usize_of(n_beta)?,
                &alpha,
                &beta,
            )?;
            determinants.push('(');
            determinants.push_str(
                &determinant
                    .words()
                    .iter()
                    .map(|word| (*word as i64).to_string())
                    .collect::<Vec<_>>()
                    .join(" "),
            );
            determinants.push(')');
            index += 4;
            continue;
        }
        index += 1;
    }
    Ok(format!("({determinants})"))
}

/// `t_of_sexp` with the exact record schema.
fn convert(source: &str) -> std::result::Result<DeterminantsByHand, SexpError> {
    let node = sexp::parse(source)?;
    let pairs = sexp::record_fields(&node, &FIELDS)?;
    let field = |name: &str| -> &Sexp {
        // `record_fields` guarantees every declared field is present.
        pairs
            .iter()
            .find(|(seen, _, _)| *seen == name)
            .map(|(_, _, value)| *value)
            .expect("record_fields checks completeness")
    };
    let pair = |name: &str| -> &Sexp {
        pairs
            .iter()
            .find(|(seen, _, _)| *seen == name)
            .map(|(_, pair, _)| *pair)
            .expect("record_fields checks completeness")
    };

    Ok(DeterminantsByHand {
        n_int: int_field("n_int", pair("n_int"), field("n_int"))?,
        bit_kind: int_field("bit_kind", pair("bit_kind"), field("bit_kind"))?,
        n_det: int_field("n_det", pair("n_det"), field("n_det"))?,
        n_det_qp_edit: int_field(
            "n_det_qp_edit",
            pair("n_det_qp_edit"),
            field("n_det_qp_edit"),
        )?,
        n_states: int_field("n_states", pair("n_states"), field("n_states"))?,
        expected_s2: float_field("expected_s2", pair("expected_s2"), field("expected_s2"))?,
        psi_coef: float_list_field("psi_coef", pair("psi_coef"), field("psi_coef"))?,
        psi_det: determinant_list_field("psi_det", pair("psi_det"), field("psi_det"))?,
        state_average_weight: float_list_field(
            "state_average_weight",
            pair("state_average_weight"),
            field("state_average_weight"),
        )?,
        mo_label: mo_label_field(pair("mo_label"), field("mo_label"))?,
    })
}

/// Report a value-conversion error on the whole `(name value)` pair.
fn value_error(message: impl Into<String>, pair: &Sexp) -> SexpError {
    let (start, end) = pair.span();
    SexpError {
        message: message.into(),
        start,
        end,
    }
}

/// `int_of_sexp`, including the OCaml literal syntax (`0x`, `0o`, `0b`, `_`).
fn int_field(name: &str, pair: &Sexp, node: &Sexp) -> std::result::Result<i64, SexpError> {
    let Some(text) = node.atom() else {
        return Err(value_error(
            format!("{name}_of_sexp: expected an atom"),
            pair,
        ));
    };
    ocaml_int_of_string(text)
        .ok_or_else(|| value_error(format!("{name}_of_sexp: invalid integer {text:?}"), pair))
}

/// `float_of_sexp`.  Rust's parser covers the decimal forms the editor
/// produces; hexadecimal floats and `_` separators are not accepted (a
/// documented deviation).
fn float_field(name: &str, pair: &Sexp, node: &Sexp) -> std::result::Result<f64, SexpError> {
    let Some(text) = node.atom() else {
        return Err(value_error(
            format!("{name}_of_sexp: expected an atom"),
            pair,
        ));
    };
    parse_ocaml_float(text)
        .ok_or_else(|| value_error(format!("{name}_of_sexp: invalid float {text:?}"), pair))
}

/// `float array_of_sexp`.
fn float_list_field(
    name: &str,
    pair: &Sexp,
    node: &Sexp,
) -> std::result::Result<Vec<f64>, SexpError> {
    let Sexp::List { items, .. } = node else {
        return Err(value_error(
            format!("{name}_of_sexp: expected a list"),
            pair,
        ));
    };
    items
        .iter()
        .map(|item| {
            let Some(text) = item.atom() else {
                return Err(value_error(
                    format!("{name}_of_sexp: expected a float atom"),
                    pair,
                ));
            };
            parse_ocaml_float(text)
                .ok_or_else(|| value_error(format!("{name}_of_sexp: invalid float {text:?}"), pair))
        })
        .collect()
}

/// `Determinant.t array_of_sexp`: a list of lists of int64 words.
fn determinant_list_field(
    name: &str,
    pair: &Sexp,
    node: &Sexp,
) -> std::result::Result<Vec<Determinant>, SexpError> {
    let Sexp::List { items, .. } = node else {
        return Err(value_error(
            format!("{name}_of_sexp: expected a list"),
            pair,
        ));
    };
    items
        .iter()
        .map(|item| {
            let Sexp::List { items: words, .. } = item else {
                return Err(value_error(
                    format!("{name}_of_sexp: expected a list of words"),
                    pair,
                ));
            };
            let mut parsed = Vec::with_capacity(words.len());
            for word in words {
                let Some(text) = word.atom() else {
                    return Err(value_error(
                        format!("{name}_of_sexp: expected an int64 atom"),
                        pair,
                    ));
                };
                let value = ocaml_int_of_string(text).ok_or_else(|| {
                    value_error(format!("{name}_of_sexp: invalid int64 {text:?}"), pair)
                })?;
                parsed.push(value as u64);
            }
            Ok(Determinant::from_words(&parsed))
        })
        .collect()
}

/// `MO_label.t_of_sexp`: the constructor atoms (or their one-element list
/// form, which `ppx_sexp_conv` also accepts).
fn mo_label_field(pair: &Sexp, node: &Sexp) -> std::result::Result<String, SexpError> {
    let atom = match node {
        Sexp::Atom { value, .. } => value.as_str(),
        Sexp::List { items, .. } if items.len() == 1 => match &items[0] {
            Sexp::Atom { value, .. } => value.as_str(),
            Sexp::List { .. } => {
                return Err(value_error("mo_label_of_sexp: unknown constructor", pair));
            }
        },
        Sexp::List { .. } => {
            return Err(value_error("mo_label_of_sexp: unknown constructor", pair));
        }
    };
    match atom {
        "Guess" | "Canonical" | "Natural" | "Localized" | "Orthonormalized" | "MCSCF" | "MOM"
        | "None" => Ok(atom.to_string()),
        other => Err(value_error(
            format!("mo_label_of_sexp: unknown constructor {other:?}"),
            pair,
        )),
    }
}

/// `int_of_string` as OCaml defines it.
fn ocaml_int_of_string(text: &str) -> Option<i64> {
    let (sign, rest) = match text.as_bytes().first() {
        Some(b'-') => (-1i64, &text[1..]),
        Some(b'+') => (1, &text[1..]),
        _ => (1, text),
    };
    let cleaned: String = rest.chars().filter(|c| *c != '_').collect();
    let (radix, digits) = if let Some(digits) = cleaned
        .strip_prefix("0x")
        .or_else(|| cleaned.strip_prefix("0X"))
    {
        (16, digits)
    } else if let Some(digits) = cleaned
        .strip_prefix("0o")
        .or_else(|| cleaned.strip_prefix("0O"))
    {
        (8, digits)
    } else if let Some(digits) = cleaned
        .strip_prefix("0b")
        .or_else(|| cleaned.strip_prefix("0B"))
    {
        (2, digits)
    } else {
        (10, cleaned.as_str())
    };
    if digits.is_empty() {
        return None;
    }
    i64::from_str_radix(digits, radix)
        .ok()
        .map(|value| sign * value)
}

/// `float_of_string` as OCaml defines it for the forms the editor emits.
fn parse_ocaml_float(text: &str) -> Option<f64> {
    text.parse::<f64>().ok()
}

/// `ezfio_array_of_list` needs a dimension; a negative one is a bug in the
/// data, not a panic.
fn usize_of(value: i64) -> Result<usize> {
    usize::try_from(value).map_err(|_| anyhow::anyhow!("negative dimension or index: {value}"))
}

// ---------------------------------------------------------------------
// Section hooks
// ---------------------------------------------------------------------

/// `get` for this section: `None` renders as `""`.
pub(crate) fn get(ezfio: &EzFio, title: &str) -> Result<String> {
    match read(ezfio, false)? {
        Some(record) => Ok(make_header(title) + &to_rst(ezfio, &record)?),
        None => Ok(String::new()),
    }
}

/// `set` for this section: parse, and apply unless `write_back` is false.
pub(crate) fn set(ezfio: &EzFio, title: &str, body: &str, write_back: bool) -> Result<()> {
    apply_parsed(
        title,
        body,
        write_back,
        |text| of_rst(ezfio, text),
        |record| write(ezfio, record, false),
    )
}

// ---------------------------------------------------------------------
// Driver actions (`-n` and `-s`)
// ---------------------------------------------------------------------

/// `update_ndet`: keep the first `n_det` determinants of every state and
/// rewrite the file with `force:true`.
///
/// The reference prints `Reducing n_det to <n>` on stdout even when the value
/// does not change, and its guard `n_det_new > n_det_new` is always false (a
/// dead check that the port keeps out).
pub(crate) fn update_ndet(ezfio: &EzFio, n_det_new: i64) -> Result<()> {
    println!("Reducing n_det to {n_det_new}");
    let mut record = match read(ezfio, true)? {
        Some(record) => record,
        None => bail!("No determinants in file"),
    };
    let n_det_old = record.n_det;
    let n_states = record.n_states;
    if n_det_new == n_det_old {
        // The reference falls through to the copy loop, which is a no-op.
    }
    for state in 0..n_states {
        let ishift_old = state * n_det_old;
        let ishift_new = state * n_det_new;
        for i in 0..n_det_new {
            let source = usize_of(i + ishift_old)?;
            let target = usize_of(i + ishift_new)?;
            let value = *record.psi_coef.get(source).ok_or_else(|| {
                anyhow::anyhow!(
                    "Invalid_argument: psi_coef index {source} out of bounds (length {})",
                    record.psi_coef.len()
                )
            })?;
            if let Some(slot) = record.psi_coef.get_mut(target) {
                *slot = value;
            } else {
                return Err(anyhow::anyhow!(
                    "Invalid_argument: psi_coef index {target} out of bounds (length {})",
                    record.psi_coef.len()
                ));
            }
        }
    }
    record.n_det = qptypes::det_number_of_int(ezfio, n_det_new)?;
    write(ezfio, &record, true)
}

/// `extract_states`: copy the requested states to the front and keep only
/// them.
///
/// The range is already sorted and deduplicated by `EditConfig` (OCaml
/// `Range.to_int_list` sorts with `sort_uniq`).  The reference prints
/// `Extracting states <range>` and one `OK` per shifted state.
pub(crate) fn extract_states(ezfio: &EzFio, states: &[i64]) -> Result<()> {
    println!("Extracting states {}", format_range(states));
    let mut record = match read(ezfio, true)? {
        Some(record) => record,
        None => bail!("No determinants in file"),
    };
    let n_det = record.n_det;
    let n_states = record.n_states;
    for state in states {
        if *state > n_states {
            bail!("State to extract should not be greater than n_states");
        }
    }
    let mut sorted = states.to_vec();
    sorted.sort_unstable();
    let mut state_shift = 0i64;
    for state in &sorted {
        let j = state - 1;
        if j > 0 {
            let ishift = j * n_det;
            for i in 0..n_det {
                let source = usize_of(i + ishift)?;
                let target = usize_of(state_shift + i)?;
                let value = *record.psi_coef.get(source).ok_or_else(|| {
                    anyhow::anyhow!(
                        "Invalid_argument: psi_coef index {source} out of bounds (length {})",
                        record.psi_coef.len()
                    )
                })?;
                if let Some(slot) = record.psi_coef.get_mut(target) {
                    *slot = value;
                } else {
                    return Err(anyhow::anyhow!(
                        "Invalid_argument: psi_coef index {target} out of bounds (length {})",
                        record.psi_coef.len()
                    ));
                }
            }
            println!("OK");
            // `Printf.printf "OK\n%!"` flushes; `println!` on stdout is line
            // buffered when attached to a terminal and flushed at exit
            // otherwise, which is equivalent here.
        }
        state_shift += n_det;
    }
    record.n_states = qptypes::states_of_int(sorted.len() as i64)?;
    write(ezfio, &record, true)
}

/// `Range.to_string`: `[a,b]` (the `EditConfig` spelling).
fn format_range(states: &[i64]) -> String {
    format!(
        "[{}]",
        states
            .iter()
            .map(i64::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The body captured from the reference for the `hcn` fixture (388 bytes).
    const HCN_BODY: &str = "\nForce the selected wave function to be an eigenfunction of S^2.\nIf true, input the expected value of S^2 ::\n\n  expected_s2 = 0.\n\nNumber of determinants ::\n\n  n_det = 1\n\nState average weights ::\n\n  state_average_weight = (1.)\n\nDeterminants ::\n\n  1.\n    +++++++---------------------------------------------------------\n    +++++++---------------------------------------------------------\n\n";

    /// A minimal EZFIO directory that reproduces the `hcn` fixture's
    /// determinant inputs: `mo_num = 64` (one word), 7 alpha and 7 beta
    /// electrons, one determinant, one state, `Canonical` MOs.
    fn hcn_like() -> (tempfile::TempDir, EzFio) {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        ezfio.set::<i64>("mo_basis", "mo_num", &64).expect("mo_num");
        ezfio
            .set::<String>("mo_basis", "mo_label", &"Canonical".to_string())
            .expect("mo_label");
        ezfio
            .set::<i64>("electrons", "elec_alpha_num", &7)
            .expect("alpha");
        ezfio
            .set::<i64>("electrons", "elec_beta_num", &7)
            .expect("beta");
        (dir, ezfio)
    }

    /// Exact decompressed payloads captured from the reference
    /// (`tests/golden/transcripts/hcn/added/determinants/`).
    const REFERENCE: [(&str, &str); 7] = [
        ("bit_kind", "                   8\n"),
        ("expected_s2", "   0.000000000000000e+00\n"),
        ("mo_label", "Canonical\n"),
        ("n_int", "                   1\n"),
        (
            "psi_coef",
            "  2\n                   1                    1 \n   1.000000000000000e+00\n",
        ),
        (
            "psi_det",
            "  3\n                   1                    2                    1 \n                 127\n                 127\n",
        ),
        (
            "state_average_weight",
            "  1\n                   1 \n   1.000000000000000e+00\n",
        ),
    ];

    fn read_raw(ezfio: &EzFio, key: &str) -> String {
        use std::io::Read;

        let path = ezfio.root().join(GROUP).join(format!("{key}.gz"));
        if path.exists() {
            let mut text = String::new();
            flate2::read::GzDecoder::new(std::fs::File::open(&path).expect("open payload"))
                .read_to_string(&mut text)
                .expect("decompress payload");
            text
        } else {
            std::fs::read_to_string(ezfio.root().join(GROUP).join(key)).expect("read payload")
        }
    }

    #[test]
    fn read_materialises_the_seven_contract_files() {
        let (_dir, ezfio) = hcn_like();

        let record = read(&ezfio, false).expect("read").expect("some");

        assert_eq!(record.n_int, 1);
        assert_eq!(record.bit_kind, 8);
        assert_eq!(record.n_det, 1);
        assert_eq!(record.n_det_qp_edit, 1);
        assert_eq!(record.n_states, 1);
        assert_eq!(record.expected_s2, 0.0);
        assert_eq!(record.psi_coef, vec![1.0]);
        assert_eq!(record.psi_det.len(), 1);
        assert_eq!(record.psi_det[0].words(), &[127, 127]);
        assert_eq!(record.state_average_weight, vec![1.0]);
        assert_eq!(record.mo_label, "Canonical");

        // Byte-for-byte against the reference payloads (C3).
        for (key, expected) in REFERENCE {
            assert_eq!(read_raw(&ezfio, key), expected, "payload {key}");
        }
        // C4: the `*_qp_edit` scratch pair is not materialised by `read`.
        assert!(!ezfio.has(GROUP, "psi_coef_qp_edit").expect("has"));
        assert!(!ezfio.has(GROUP, "psi_det_qp_edit").expect("has"));
    }

    #[test]
    fn to_rst_matches_the_reference_body_byte_for_byte() {
        let (_dir, ezfio) = hcn_like();
        let record = read(&ezfio, false).expect("read").expect("some");
        assert_eq!(to_rst(&ezfio, &record).expect("render"), HCN_BODY);
    }

    #[test]
    fn of_rst_round_trips_the_reference_body() {
        let (_dir, ezfio) = hcn_like();
        let record = read(&ezfio, false).expect("read").expect("some");
        // `set` hands `of_rst` the text with one trailing newline removed.
        let body = HCN_BODY.strip_suffix('\n').expect("one newline");
        let parsed = of_rst(&ezfio, body).expect("parse").expect("some");
        assert_eq!(parsed.n_int, record.n_int);
        assert_eq!(parsed.bit_kind, record.bit_kind);
        assert_eq!(parsed.n_det, record.n_det);
        assert_eq!(parsed.n_det_qp_edit, record.n_det_qp_edit);
        assert_eq!(parsed.n_states, record.n_states);
        assert_eq!(parsed.psi_coef, record.psi_coef);
        assert_eq!(parsed.psi_det, record.psi_det);
        assert_eq!(parsed.state_average_weight, record.state_average_weight);
        assert_eq!(parsed.mo_label, record.mo_label);
    }

    #[test]
    fn of_rst_takes_the_header_values_from_the_edited_text() {
        let (_dir, ezfio) = hcn_like();
        let body = HCN_BODY.replace("  expected_s2 = 0.", "  expected_s2 = 2.5");
        let parsed = of_rst(&ezfio, &body).expect("parse").expect("some");
        assert_eq!(parsed.expected_s2, 2.5);
    }

    #[test]
    fn of_rst_reports_a_missing_determinants_header_and_returns_none() {
        let (_dir, ezfio) = hcn_like();
        let error = of_rst(&ezfio, "expected_s2 = 0.\n").expect_err("assertion escapes");
        assert!(crate::ocaml_compat::is_assert_failure(&error));
    }

    #[test]
    fn of_rst_reports_a_duplicate_field_and_writes_nothing() {
        let (_dir, ezfio) = hcn_like();
        // The header filter keeps the second `expected_s2` line too, and the
        // record reader rejects the duplicate.
        let body = HCN_BODY.replace(
            "  expected_s2 = 0.",
            "  expected_s2 = 0.\n  expected_s2 = 1.",
        );
        assert!(of_rst(&ezfio, &body).expect("reported").is_none());
    }

    /// `write(force:false)` with `n_det <= n_det_qp_edit` writes the scratch
    /// pair too, and the seven contract files are unchanged.
    #[test]
    fn write_force_false_writes_the_working_and_editor_copies() {
        let (_dir, ezfio) = hcn_like();
        let record = read(&ezfio, false).expect("read").expect("some");
        write(&ezfio, &record, false).expect("write");
        assert_eq!(
            read_raw(&ezfio, "psi_coef_qp_edit"),
            REFERENCE[4].1,
            "psi_coef_qp_edit"
        );
        assert_eq!(
            read_raw(&ezfio, "psi_det_qp_edit"),
            REFERENCE[5].1,
            "psi_det_qp_edit"
        );
        for (key, expected) in REFERENCE {
            assert_eq!(read_raw(&ezfio, key), expected, "payload {key}");
        }
    }

    /// A record whose `n_det` exceeds `n_det_qp_edit` must not resurrect the
    /// truncated determinants when `force` is false.
    #[test]
    fn write_force_false_skips_the_arrays_when_det_count_shrinks() {
        let (_dir, ezfio) = hcn_like();
        let mut record = read(&ezfio, false).expect("read").expect("some");
        ezfio.set::<i64>(GROUP, "n_det", &2).expect("n_det");
        record.n_det = 2;
        record.psi_coef = vec![1.0, 0.5];
        write(&ezfio, &record, false).expect("write");
        // `n_det` is written, the arrays are not: the stored payload is still
        // the one-determinant, one-state array from `read`.
        assert_eq!(read_raw(&ezfio, "psi_coef"), REFERENCE[4].1);
    }

    #[test]
    fn update_ndet_keeps_the_first_determinants_of_every_state() {
        let (_dir, ezfio) = hcn_like();
        let mut record = read(&ezfio, false).expect("read").expect("some");
        record.n_det = 3;
        record.n_states = 2;
        // `state_average_weight` always has one entry per state; `write`
        // truncates it to the (unchanged) `n_states`.
        record.state_average_weight = vec![0.5, 0.5];
        record.psi_coef = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        record.psi_det = vec![
            Determinant::from_words(&[127, 127]),
            Determinant::from_words(&[191, 127]),
            Determinant::from_words(&[223, 127]),
        ];
        write(&ezfio, &record, true).expect("write");

        update_ndet(&ezfio, 2).expect("update_ndet");

        let after = read(&ezfio, true).expect("read").expect("some");
        assert_eq!(after.n_det, 2);
        assert_eq!(after.psi_coef, vec![0.1, 0.2, 0.4, 0.5]);
    }

    #[test]
    fn extract_states_copies_the_requested_states_to_the_front() {
        let (_dir, ezfio) = hcn_like();
        let mut record = read(&ezfio, false).expect("read").expect("some");
        record.n_det = 2;
        record.n_states = 3;
        record.state_average_weight = vec![1.0 / 3.0; 3];
        record.psi_coef = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        record.psi_det = vec![
            Determinant::from_words(&[127, 127]),
            Determinant::from_words(&[191, 127]),
        ];
        write(&ezfio, &record, true).expect("write");

        extract_states(&ezfio, &[2, 3]).expect("extract");

        let after = read(&ezfio, true).expect("read").expect("some");
        assert_eq!(after.n_states, 2);
        assert_eq!(after.psi_coef, vec![0.3, 0.4, 0.5, 0.6]);
    }

    #[test]
    fn extract_states_rejects_a_state_above_n_states() {
        let (_dir, ezfio) = hcn_like();
        let error = extract_states(&ezfio, &[2]).expect_err("state 2 > n_states 1");
        assert_eq!(
            error.to_string(),
            "State to extract should not be greater than n_states"
        );
    }
}
