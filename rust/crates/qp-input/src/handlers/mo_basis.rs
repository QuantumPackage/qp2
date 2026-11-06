//! `Mo_basis`: the molecular orbitals, their label, classes and occupations.
//!
//! Port of `ocaml/Input_mo_basis.ml`.  The EZFIO group is `mo_basis`:
//!
//! | field | EZFIO key | element type |
//! |---|---|---|
//! | `mo_num` | `mo_basis/mo_num` | integer scalar |
//! | `mo_label` | `mo_basis/mo_label` | `character*(64)` scalar |
//! | `mo_class` | `mo_basis/mo_class` | `character*(32)` array, length `mo_num` |
//! | `mo_occ` | `mo_basis/mo_occ` | real array, length `mo_num` |
//! | `mo_coef` | `mo_basis/mo_coef` | real array, `(ao_num,mo_num)` |
//! | `mo_coef_imag` | `mo_basis/mo_coef_imag` | real array, same shape, optional |
//! | `ao_md5` | `mo_basis/ao_md5` | `character*(32)` scalar |
//!
//! `mo_coef` is stored **column major**: the OCaml reader divides the flat
//! payload into `mo_num` consecutive blocks of `ao_num = len / mo_num`
//! values, so block `j` is the coefficient vector of MO `j` over all AOs.
//! [`write`] rebuilds the same flat order by concatenating the rows.
//!
//! # This handler has no `of_rst`
//!
//! The reference module exposes `read`, `write`, `reorder`, `to_string` and
//! `to_rst` only.  `qp_edit.set` therefore skips the section (there is no
//! parser to run), and the RST block is documentation: editing it does
//! nothing.  `set` in [`super::Handler`] is the corresponding no-op.
//!
//! # Read paths do not restamp the MD5
//!
//! Two stamps are involved.  `Input_ao_basis.read` rewrites
//! `ao_basis/ao_md5`; decision D4 removes that write (see
//! [`super::ao_basis`]).  `read_ao_md5` here rewrites `mo_basis/ao_md5` when
//! the key is absent; the port does not write it either.  An absent key means
//! "produced by this same AO basis" — the reference writes the computed value
//! and then reads it back, so the comparison can only succeed — hence the
//! port treats absence as equal and returns the computed digest.  A present
//! key is still validated and compared, and a mismatch fails with the
//! reference message.
//!
//! Units: occupations are electrons per MO, in `[0, 2]`.

use anyhow::{Context, Result, anyhow, bail};
use qp_ezfio::EzFio;
use qp_types::md5::MD5_OF_NONE;

use super::ao_basis;
use super::qptypes;
use crate::engine::make_header;
use crate::ocaml_compat::{MoLabel, assertion};

/// The EZFIO group this handler owns.
const GROUP: &str = "mo_basis";

/// One `MO_class.t`: a keyword and, when it was built from a range, the list
/// of MO indices it names.
///
/// The OCaml type carries a `MO_number.t list` in every constructor; EZFIO
/// only ever stores the bare keyword (the reader accepts nothing else), but
/// [`MoClass::to_string`] has to reproduce the padded `"Core     : (1, 2)"`
/// form because `write` uses it and the round trip must be stable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoClass {
    /// Frozen core orbitals.
    Core(Vec<i64>),
    /// Inactive (doubly occupied, not correlated) orbitals.
    Inactive(Vec<i64>),
    /// Active orbitals.
    Active(Vec<i64>),
    /// Virtual orbitals.
    Virtual(Vec<i64>),
    /// Deleted orbitals.
    Deleted(Vec<i64>),
}

/// The printed form (`MO_class.to_string` in the reference): the bare keyword,
/// or the padded `"Core     : (1, 2)"` when the class carries indices.
///
/// The padding is part of the format: the five spellings have different
/// lengths and the reference pads them to the width of `"Inactive : "`.
/// `MO_number.to_int` is the identity on the stored indices.
impl std::fmt::Display for MoClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (name, indices) = match self {
            MoClass::Core(indices) => ("Core", indices),
            MoClass::Inactive(indices) => ("Inactive", indices),
            MoClass::Active(indices) => ("Active", indices),
            MoClass::Virtual(indices) => ("Virtual", indices),
            MoClass::Deleted(indices) => ("Deleted", indices),
        };
        if indices.is_empty() {
            return f.write_str(name);
        }
        let list = indices
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        write!(f, "{name:<8} : ({list})")
    }
}

impl MoClass {
    /// `MO_class.of_string`: lower-case the whole string, accept exactly the
    /// five keywords, and keep an empty index list.
    ///
    /// There is no `trim`: EZFIO character arrays are padded, but the crate's
    /// `Element for String` already strips the padding when it reads a field,
    /// so a direct call with a padded string fails here, as it does in OCaml.
    pub(crate) fn of_string(s: &str) -> Result<MoClass> {
        match s.to_lowercase().as_str() {
            "core" => Ok(MoClass::Core(Vec::new())),
            "inactive" => Ok(MoClass::Inactive(Vec::new())),
            "active" => Ok(MoClass::Active(Vec::new())),
            "virtual" => Ok(MoClass::Virtual(Vec::new())),
            "deleted" => Ok(MoClass::Deleted(Vec::new())),
            _ => bail!("MO_class should be (Core|Inactive|Active|Virtual|Deleted)"),
        }
    }
}

/// One molecular-orbital set, in the field order of the OCaml record.
#[derive(Debug, Clone, PartialEq)]
pub struct MoBasis {
    /// Number of molecular orbitals.
    pub mo_num: i64,
    /// How the orbitals were produced.
    pub(crate) mo_label: MoLabel,
    /// One class per MO.
    pub mo_class: Vec<MoClass>,
    /// One occupation per MO.
    pub(crate) mo_occ: Vec<f64>,
    /// `mo_num` rows of `ao_num` coefficients, in MO order.
    pub(crate) mo_coef: Vec<Vec<f64>>,
    /// Imaginary part of the coefficients (complex wave functions), if any.
    pub(crate) mo_coef_imag: Option<Vec<Vec<f64>>>,
    /// Digest of the AO basis these MOs were computed with.
    pub(crate) ao_md5: String,
}

/// `read_mo_num`: the MO count, checked against the alpha-electron count and
/// against `MO_number`'s own bounds.
///
/// The electron check comes **before** the wrapper: `MO_number.of_int` would
/// also reject 0, but the reference reports the missing-orbital case with its
/// own message and the CLI prints that one verbatim.
fn read_mo_num(ezfio: &EzFio) -> Result<i64> {
    let elec_alpha_num = ezfio.get::<i64>("electrons", "elec_alpha_num")?;
    let result = ezfio.get::<i64>(GROUP, "mo_num")?;
    if result < elec_alpha_num {
        bail!("More alpha electrons than MOs");
    }
    qptypes::mo_number_of_int(ezfio, result)
}

/// `read_mo_label`: materialise `"None"` when the key is absent, then parse.
///
/// The write is contract C3 (a defaulted attribute is created); it is not an
/// MD5 restamp.  `MO_label.of_string` echoes a rejected spelling to stdout.
fn read_mo_label(ezfio: &EzFio) -> Result<MoLabel> {
    if !ezfio.has(GROUP, "mo_label")? {
        ezfio.set::<String>(GROUP, "mo_label", &"None".to_owned())?;
    }
    MoLabel::of_string(&ezfio.get::<String>(GROUP, "mo_label")?)
}

/// `read_mo_class`: materialise one `"Active"` per MO when the key is absent,
/// then parse every entry.
fn read_mo_class(ezfio: &EzFio) -> Result<Vec<MoClass>> {
    if !ezfio.has(GROUP, "mo_class")? {
        let mo_num = usize_of(read_mo_num(ezfio)?)?;
        let data = vec![MoClass::Active(Vec::new()).to_string(); mo_num];
        let array = super::ezfio_array_of_list(vec![mo_num], data)?;
        ezfio.set_array::<String>(GROUP, "mo_class", &array)?;
    }
    let array = ezfio.get_array::<String>(GROUP, "mo_class")?;
    array.data().iter().map(|s| MoClass::of_string(s)).collect()
}

/// `read_mo_occ`: one occupation per MO, clamped to `[0, 2]`.
///
/// The reference guards the materialisation with `not (has mo_label)`, which
/// is always false once `read_mo_label` has run; the guard is reproduced so
/// that the field order stays explicit, but it never fires inside [`read`].
fn read_mo_occ(ezfio: &EzFio) -> Result<Vec<f64>> {
    if !ezfio.has(GROUP, "mo_label")? {
        let elec_alpha_num = ezfio.get::<i64>("electrons", "elec_alpha_num")?;
        let elec_beta_num = ezfio.get::<i64>("electrons", "elec_beta_num")?;
        let mo_num = usize_of(read_mo_num(ezfio)?)?;
        let data = (0..mo_num)
            .map(|i| {
                let i = i as i64;
                if i < elec_beta_num {
                    2.0
                } else if i < elec_alpha_num {
                    1.0
                } else {
                    0.0
                }
            })
            .collect::<Vec<_>>();
        let array = super::ezfio_array_of_list(vec![mo_num], data)?;
        ezfio.set_array::<f64>(GROUP, "mo_occ", &array)?;
    }
    let array = ezfio.get_array::<f64>(GROUP, "mo_occ")?;
    Ok(array
        .data()
        .iter()
        .copied()
        .map(qptypes::mo_occ_of_float)
        .collect())
}

/// Split a flat coefficient payload into `mo_num` rows of `len / mo_num`.
///
/// This is the OCaml `let ao_num = Array.length a / mo_num in Array.sub`.
/// The integer division is the reference's own way of ignoring a trailing
/// partial column; it can never index past the end because
/// `mo_num * ao_num <= len`.
fn split_rows(data: &[f64], mo_num: usize) -> Result<Vec<Vec<f64>>> {
    if mo_num == 0 {
        // `MO_number.of_int` already rejects 0; a direct call is a bug.
        return Err(assertion(
            "split_rows: assertion (mo_num > 0) failed for mo_num = 0",
        ));
    }
    let ao_num = data.len() / mo_num;
    Ok((0..mo_num)
        .map(|j| data[j * ao_num..(j + 1) * ao_num].to_vec())
        .collect())
}

/// `read_mo_coef`: `mo_num` rows of `ao_num` coefficients.
fn read_mo_coef(ezfio: &EzFio) -> Result<Vec<Vec<f64>>> {
    let data = ezfio.get_array::<f64>(GROUP, "mo_coef")?.into_data();
    split_rows(&data, usize_of(read_mo_num(ezfio)?)?)
}

/// `read_mo_coef_imag`: the same shape, or `None` when the key is absent.
fn read_mo_coef_imag(ezfio: &EzFio) -> Result<Option<Vec<Vec<f64>>>> {
    if !ezfio.has(GROUP, "mo_coef_imag")? {
        return Ok(None);
    }
    let data = ezfio.get_array::<f64>(GROUP, "mo_coef_imag")?.into_data();
    Ok(Some(split_rows(&data, usize_of(read_mo_num(ezfio)?)?)?))
}

/// `read_ao_md5`: the stored digest must describe the current AO basis.
///
/// *No write* when the key is absent (decision D4); the computed digest is
/// used as the comparison value instead, which is what the reference compares
/// after stamping it.
fn read_ao_md5(ezfio: &EzFio) -> Result<String> {
    let ao_md5 = match ao_basis::read(ezfio)? {
        None => MD5_OF_NONE.to_owned(),
        Some(result) => ao_basis::to_md5(&result)?,
    };
    let result = if !ezfio.has(GROUP, "ao_md5")? {
        ao_md5.clone()
    } else {
        qptypes::md5_of_string(&ezfio.get::<String>(GROUP, "ao_md5")?)?
    };
    if ao_md5 != result {
        bail!("The current MOs don't correspond to the current AOs.");
    }
    Ok(result)
}

/// `read`: `None` when `mo_basis/mo_num` is absent, every other failure
/// propagates.
///
/// Unlike the AO handlers there is no `try ... with _ -> None`: a missing
/// `mo_occ`, a mismatch with the AO basis or an unparsable label aborts the
/// editor with the `=-` banner.
pub fn read(ezfio: &EzFio) -> Result<Option<MoBasis>> {
    if !ezfio.has(GROUP, "mo_num")? {
        return Ok(None);
    }
    // Field order inside the record literal = declaration order, the same
    // order the reference uses.  The order matters: `read_mo_label` writes
    // the default that disables `read_mo_occ`'s materialisation.
    Ok(Some(MoBasis {
        mo_num: read_mo_num(ezfio)?,
        mo_label: read_mo_label(ezfio)?,
        mo_class: read_mo_class(ezfio)?,
        mo_occ: read_mo_occ(ezfio)?,
        mo_coef: read_mo_coef(ezfio)?,
        mo_coef_imag: read_mo_coef_imag(ezfio)?,
        ao_md5: read_ao_md5(ezfio)?,
    }))
}

/// `reorder b ordering`: rebuild every row as `row[ordering[i]]`.
///
/// The output row has `ordering.len()` entries, not the input length: the
/// reference uses `Array.init (Array.length mo)`, so a short `ordering`
/// silently truncates a row and a long one reads out of bounds.  `mo_num`,
/// the label, the classes, the occupations and the digest are untouched.
///
/// # Errors
///
/// `"index out of bounds"` when `ordering` names a column the row does not
/// have (the OCaml `Invalid_argument` from `mo.(ordering.(i))`).
pub(crate) fn reorder(b: &MoBasis, ordering: &[i64]) -> Result<MoBasis> {
    fn reorder_rows(rows: &[Vec<f64>], ordering: &[i64]) -> Result<Vec<Vec<f64>>> {
        rows.iter()
            .map(|row| {
                ordering
                    .iter()
                    .map(|&index| {
                        let index =
                            usize::try_from(index).map_err(|_| anyhow!("index out of bounds"))?;
                        row.get(index)
                            .copied()
                            .ok_or_else(|| anyhow!("index out of bounds"))
                    })
                    .collect()
            })
            .collect()
    }
    Ok(MoBasis {
        mo_coef: reorder_rows(&b.mo_coef, ordering)?,
        mo_coef_imag: b
            .mo_coef_imag
            .as_deref()
            .map(|rows| reorder_rows(rows, ordering))
            .transpose()?,
        ..b.clone()
    })
}

/// `write`: the seven fields, in record order.
///
/// `write_mo_num` is the only stamp of a derived scalar here: `ao_num` comes
/// from the first row of `mo_coef`, as in the reference, which then indexes
/// `a.(0)` and raises when the matrix is empty.
pub(crate) fn write(ezfio: &EzFio, b: &MoBasis) -> Result<()> {
    ezfio.set::<i64>(GROUP, "mo_num", &b.mo_num)?;
    ezfio.set::<String>(GROUP, "mo_label", &b.mo_label.to_string())?;

    let class = b
        .mo_class
        .iter()
        .map(MoClass::to_string)
        .collect::<Vec<_>>();
    let class = super::ezfio_array_of_list(vec![b.mo_class.len()], class)?;
    ezfio.set_array::<String>(GROUP, "mo_class", &class)?;

    let occ = super::ezfio_array_of_list(vec![b.mo_occ.len()], b.mo_occ.clone())?;
    ezfio.set_array::<f64>(GROUP, "mo_occ", &occ)?;

    let ao_num = b
        .mo_coef
        .first()
        .context("Invalid_argument: index out of bounds")?
        .len();
    let coef = b.mo_coef.iter().flatten().copied().collect::<Vec<_>>();
    let coef = super::ezfio_array_of_list(vec![ao_num, b.mo_coef.len()], coef)?;
    ezfio.set_array::<f64>(GROUP, "mo_coef", &coef)?;

    if let Some(imag) = &b.mo_coef_imag {
        let ao_num = imag
            .first()
            .context("Invalid_argument: index out of bounds")?
            .len();
        let data = imag.iter().flatten().copied().collect::<Vec<_>>();
        let imag = super::ezfio_array_of_list(vec![ao_num, imag.len()], data)?;
        ezfio.set_array::<f64>(GROUP, "mo_coef_imag", &imag)?;
    }

    ezfio.set::<String>(GROUP, "ao_md5", &b.ao_md5)?;
    Ok(())
}

/// `print_five` for one block of one to four MOs.
///
/// `%15.10f` is fixed notation with ten decimals, right justified in fifteen
/// columns; Rust's `{:15.10}` on `f64` produces the same text.  The
/// single-MO case keeps the reference's trailing space — it is part of the
/// golden RST (`ref_mo.rst` ends with one).
fn print_block(mo_coef: &[Vec<f64>], ao_num: usize, imin: usize, cols: usize) -> Result<String> {
    let mut lines = Vec::with_capacity(1 + ao_num);
    let mut header = format!("  #{:>15}", imin + 1);
    for k in 1..cols {
        header.push_str(&format!(" {:>15}", imin + k + 1));
    }
    lines.push(header);
    for row in 0..ao_num {
        let mut line = format!("  {:>3}", row + 1);
        for k in 0..cols {
            let value = mo_coef[imin + k]
                .get(row)
                .context("Invalid_argument: index out of bounds")?;
            line.push_str(&format!(" {value:15.10}"));
        }
        if cols == 1 {
            line.push(' ');
        }
        lines.push(line);
    }
    Ok(lines.join("\n"))
}

/// `mo_coef_to_string`: blocks of four MOs, separated by a blank line.
///
/// The block loop is `if i + 4 < mo_num then i..i+3 else i..mo_num-1`, so the
/// last block holds one to four columns (a multiple of four leaves the last
/// block full).  The reference indexes `mo_coef.(0)` before the loop, which
/// is where an empty matrix fails.
fn mo_coef_to_string(mo_coef: &[Vec<f64>]) -> Result<String> {
    let ao_num = mo_coef
        .first()
        .context("Invalid_argument: index out of bounds")?
        .len();
    let mo_num = mo_coef.len();
    let mut blocks = Vec::new();
    let mut i = 0;
    loop {
        let cols = if i + 4 < mo_num { 4 } else { mo_num - i };
        blocks.push(print_block(mo_coef, ao_num, i, cols)?);
        i += cols;
        if i >= mo_num {
            break;
        }
    }
    Ok(blocks.join("\n\n"))
}

/// `to_rst`: the editable block, without the section header.
///
/// The final `%s\n` is the only newline after the last coefficient line.
pub fn to_rst(b: &MoBasis) -> Result<String> {
    Ok(format!(
        "\nLabel of the molecular orbitals ::\n\n  mo_label = {}\n\n\
         Total number of MOs ::\n\n  mo_num = {}\n\n\
         MO coefficients ::\n\n{}\n",
        b.mo_label,
        b.mo_num,
        mo_coef_to_string(&b.mo_coef)?
    ))
}

/// `get` for this section: `None` renders as `""`.
pub(crate) fn get(ezfio: &EzFio, title: &str) -> Result<String> {
    match read(ezfio)? {
        Some(record) => Ok(make_header(title) + &to_rst(&record)?),
        None => Ok(String::new()),
    }
}

/// A negative count is a corrupt file, not a panic.
fn usize_of(value: i64) -> Result<usize> {
    usize::try_from(value).map_err(|_| anyhow!("negative dimension or index: {value}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::handlers::test_support::hcn_ezfio;

    /// The whole MO section, byte for byte against the reference binary's
    /// output.
    #[test]
    fn get_renders_the_reference_section_byte_for_byte() {
        let (_dir, ezfio) = hcn_ezfio();
        let rendered = get(&ezfio, "MO basis").expect("get");
        assert_eq!(rendered.len(), 57_952, "section length");
        assert_eq!(
            qp_types::md5::md5_hex(&rendered),
            "07f436627be74e44e7ba9127a2d1d633",
            "the MO basis section no longer matches the reference"
        );
    }

    /// The read path adds no `mo_basis/ao_md5` when the key is missing and
    /// does not touch the one that is there.
    #[test]
    fn read_does_not_restamp_the_ao_md5() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        // Minimal, consistent basis: one S function.
        let ao = ao_basis::AoBasis {
            ao_basis: "test".to_owned(),
            ao_num: 1,
            ao_prim_num: vec![1],
            ao_prim_num_max: 1,
            ao_nucl: vec![1],
            ao_power: vec![qp_types::gto::AngmomXyz::new(0, 0, 0)],
            ao_coef: vec![1.0],
            ao_expo: vec![2.0],
            ao_cartesian: false,
        };
        for (group, key, value) in [
            ("nuclei", "nucl_num", 1i64),
            ("electrons", "elec_alpha_num", 1),
            ("electrons", "elec_beta_num", 1),
        ] {
            ezfio.set::<i64>(group, key, &value).expect("scalar");
        }
        ao_basis::write(&ezfio, &ao).expect("write AO");
        // `ao_basis::write` does not stamp the derived `ao_num`; the
        // reference reads it.
        ezfio.set::<i64>("ao_basis", "ao_num", &1).expect("ao_num");
        ezfio.set::<i64>(GROUP, "mo_num", &1).expect("mo_num");
        ezfio
            .set::<String>(GROUP, "mo_label", &"Canonical".to_owned())
            .expect("label");
        let class =
            super::super::ezfio_array_of_list(vec![1], vec!["Active".to_owned()]).expect("class");
        ezfio
            .set_array::<String>(GROUP, "mo_class", &class)
            .expect("class write");
        let occ = super::super::ezfio_array_of_list(vec![1], vec![2.0]).expect("occ");
        ezfio
            .set_array::<f64>(GROUP, "mo_occ", &occ)
            .expect("occ write");
        let coef = super::super::ezfio_array_of_list(vec![1, 1], vec![1.0]).expect("coef");
        ezfio
            .set_array::<f64>(GROUP, "mo_coef", &coef)
            .expect("coef write");
        assert!(!ezfio.has(GROUP, "ao_md5").expect("has"));

        let record = read(&ezfio).expect("read").expect("some");
        assert_eq!(record.ao_md5, ao_basis::to_md5(&ao).expect("md5"));
        assert!(
            !ezfio.has(GROUP, "ao_md5").expect("has"),
            "read stamped mo_basis/ao_md5"
        );
    }

    #[test]
    fn read_returns_none_without_mo_num() {
        let (_dir, ezfio) = hcn_ezfio();
        ezfio.remove(GROUP, "mo_num").expect("remove");
        assert!(read(&ezfio).expect("read").is_none());
    }

    /// A digest that does not describe the current AO basis is fatal: the MO
    /// handler has no `try ... with _ -> None`.
    #[test]
    fn read_rejects_a_mismatched_ao_md5() {
        let (_dir, ezfio) = hcn_ezfio();
        ezfio
            .set::<String>(GROUP, "ao_md5", &MD5_OF_NONE.to_owned())
            .expect("set");
        let error = read(&ezfio).expect_err("mismatch");
        assert_eq!(
            error.to_string(),
            "The current MOs don't correspond to the current AOs."
        );
    }

    /// `read_mo_class` creates one `Active` per MO and `read_mo_label` one
    /// `None`; `read_mo_occ` never creates anything (its guard is already
    /// false).
    #[test]
    fn read_materialises_the_label_and_the_classes() {
        let (_dir, ezfio) = hcn_ezfio();
        ezfio.remove(GROUP, "mo_label").expect("remove label");
        ezfio.remove(GROUP, "mo_class").expect("remove class");
        let record = read(&ezfio).expect("read").expect("some");
        assert_eq!(record.mo_label, MoLabel::None);
        assert_eq!(record.mo_class.len(), 55);
        assert!(
            record
                .mo_class
                .iter()
                .all(|c| c == &MoClass::Active(Vec::new()))
        );
        assert_eq!(
            ezfio.get::<String>(GROUP, "mo_label").expect("label"),
            "None"
        );
    }

    /// `read_mo_occ`'s guard is `not (has mo_label)`, and `read_mo_label`
    /// has just created the key, so the occupation array is never
    /// materialised: a missing `mo_occ` is an error, not a default.
    #[test]
    fn read_does_not_materialise_mo_occ() {
        let (_dir, ezfio) = hcn_ezfio();
        ezfio.remove(GROUP, "mo_occ").expect("remove occ");
        let error = read(&ezfio).expect_err("missing occ");
        assert!(error.to_string().contains("mo_basis/mo_occ"), "{error}");
        assert!(
            !ezfio.has(GROUP, "mo_occ").expect("has"),
            "mo_occ was created"
        );
    }

    #[test]
    fn mo_class_round_trips_through_its_printed_form() {
        assert_eq!(MoClass::Active(Vec::new()).to_string(), "Active");
        assert_eq!(
            MoClass::Core(vec![1, 2, 3]).to_string(),
            "Core     : (1, 2, 3)"
        );
        assert_eq!(MoClass::Inactive(vec![4]).to_string(), "Inactive : (4)");
        assert_eq!(
            MoClass::of_string("aCtIvE").expect("parse").to_string(),
            "Active"
        );
        assert_eq!(
            MoClass::of_string(" active")
                .expect_err("no trim")
                .to_string(),
            "MO_class should be (Core|Inactive|Active|Virtual|Deleted)"
        );
    }

    /// Blocks of four: 5 MOs give a full block and a single-MO block with
    /// the reference's trailing space, 4 give one full block, 1 one padded
    /// line.
    #[test]
    fn mo_coef_to_string_uses_blocks_of_four() {
        // `mo_coef` is one row per MO, so `[vec![1.0, 1.5]]` is one MO over
        // two AOs: a header, then two AO lines.
        let one = mo_coef_to_string(&[vec![1.0, 1.5]]).expect("1 MO");
        assert_eq!(
            one,
            "  #              1\n    1    1.0000000000 \n    2    1.5000000000 "
        );

        // Four MOs over one AO: a single full block, no blank line.
        let four = mo_coef_to_string(&[vec![1.0], vec![2.0], vec![3.0], vec![4.0]]).expect("4 MOs");
        assert_eq!(four.lines().count(), 2);
        assert_eq!(
            four.lines().next().expect("header"),
            "  #              1               2               3               4"
        );
        assert_eq!(
            four.lines().nth(1).expect("row"),
            "    1    1.0000000000    2.0000000000    3.0000000000    4.0000000000"
        );

        // Five: a full block, a blank line, then the one-MO block with the
        // reference's trailing space.
        let five = mo_coef_to_string(&[vec![1.0], vec![2.0], vec![3.0], vec![4.0], vec![5.0]])
            .expect("5 MOs");
        // `lines()` yields the empty line between the two blocks too.
        assert_eq!(five.lines().count(), 5, "header+row, blank, header+row");
        assert_eq!(
            five.lines().last().expect("last line"),
            "    1    5.0000000000 "
        );
        assert!(five.contains("\n\n  #              5\n"));
    }

    #[test]
    fn reorder_permutes_every_coefficient_row() {
        let record = MoBasis {
            mo_num: 2,
            mo_label: MoLabel::Canonical,
            mo_class: vec![MoClass::Active(Vec::new()); 2],
            mo_occ: vec![2.0, 0.0],
            mo_coef: vec![vec![1.0, 2.0], vec![3.0, 4.0]],
            mo_coef_imag: Some(vec![vec![5.0, 6.0], vec![7.0, 8.0]]),
            ao_md5: MD5_OF_NONE.to_owned(),
        };
        let reordered = reorder(&record, &[1, 0]).expect("reorder");
        assert_eq!(reordered.mo_coef, vec![vec![2.0, 1.0], vec![4.0, 3.0]]);
        assert_eq!(
            reordered.mo_coef_imag,
            Some(vec![vec![6.0, 5.0], vec![8.0, 7.0]])
        );
        // Everything else is untouched.
        assert_eq!(reordered.mo_occ, record.mo_occ);
        assert_eq!(reordered.mo_num, record.mo_num);
        let error = reorder(&record, &[0, 2]).expect_err("out of bounds");
        assert_eq!(error.to_string(), "index out of bounds");
    }

    #[test]
    fn write_round_trips_through_read() {
        let (_dir, ezfio) = hcn_ezfio();
        let record = read(&ezfio).expect("read").expect("some");
        write(&ezfio, &record).expect("write");
        let again = read(&ezfio).expect("read").expect("some");
        assert_eq!(record, again);
    }

    /// `write_mo_coef` reads `a.(0)`: an empty matrix is the OCaml
    /// `Invalid_argument`, not a panic.
    #[test]
    fn write_rejects_an_empty_coefficient_matrix() {
        let (_dir, ezfio) = hcn_ezfio();
        let mut record = read(&ezfio).expect("read").expect("some");
        record.mo_coef.clear();
        let error = write(&ezfio, &record).expect_err("empty");
        assert_eq!(error.to_string(), "Invalid_argument: index out of bounds");
    }
}
