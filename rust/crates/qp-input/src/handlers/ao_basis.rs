//! `Ao_basis`: the AO basis set, its MD5, its ordering and its RST block.
//!
//! Port of `ocaml/Input_ao_basis.ml`.  The EZFIO group is `ao_basis`:
//!
//! | field | EZFIO key | element type |
//! |---|---|---|
//! | `ao_basis` | `ao_basis/ao_basis` | `character*(256)` scalar |
//! | `ao_num` | `ao_basis/ao_num` | integer scalar |
//! | `ao_prim_num` | `ao_basis/ao_prim_num` | integer array, length `ao_num` |
//! | `ao_prim_num_max` | derived | integer scalar, `max(ao_prim_num)` |
//! | `ao_nucl` | `ao_basis/ao_nucl` | integer array, length `ao_num` |
//! | `ao_power` | `ao_basis/ao_power` | integer array, `(ao_num,3)` |
//! | `ao_coef` | `ao_basis/ao_coef` | real array, `(ao_num,ao_prim_num_max)` |
//! | `ao_expo` | `ao_basis/ao_expo` | real array, `(ao_num,ao_prim_num_max)` |
//! | `ao_cartesian` | `ao_basis/ao_cartesian` | logical scalar |
//!
//! `ao_prim_num_max` has **no EZFIO key of its own** in this handler: the
//! reference recomputes it from `ao_prim_num` on every read, although
//! `EZFIO.cfg` declares it with `default: =maxval(ao_basis.ao_prim_num)`.
//! The `hcn` fixture does not contain the key, and `read` still works.
//!
//! # One entry per Cartesian function, not per shell
//!
//! `ao_num` counts *Cartesian functions*: a `P` shell occupies three entries
//! of `ao_power`, each with its own `ao_prim_num`, `ao_coef` and `ao_expo`
//! row.  [`to_long_basis`] therefore builds a one-primitive-list `Gto` per
//! entry, and [`LongBasis::to_basis`] keeps only the entry whose monomial is
//! the shell's first (`x` for `P`, `xx` for `D`, ...).  That collapse is what
//! makes the MD5 independent of the number of Cartesian components.
//!
//! # Read paths do not restamp the MD5
//!
//! The OCaml `read` writes `ao_basis/ao_md5` twice: with the computed digest
//! on success, and with `MD5("None")` on failure.  Decision D4 removes both:
//! reading a calculation must not modify it, and check mode (contract C2)
//! compares the directory byte for byte.  `read` still *computes* the digest,
//! because the reference's `try` covers that computation: a record whose
//! digest cannot be built reads as `None`.  The stamp is written only by
//! [`write`].
//!
//! Units: exponents are in bohr⁻², coefficients are dimensionless.

use anyhow::{Context, Result, anyhow, bail};
use qp_ezfio::EzFio;
use qp_types::basis::{self, LongBasis};
use qp_types::gto::{AngmomXyz, Fmt, GaussianPrimitive, Gto, GtoError};

use super::qptypes;
use crate::engine::make_header;
use crate::ocaml_compat::assertion;

/// The EZFIO group this handler owns.
const GROUP: &str = "ao_basis";

/// One AO basis set, in the field order of the OCaml record.
///
/// The vectors are the flattened EZFIO arrays, in Fortran (column-major)
/// order: `ao_coef[i*ao_num + j]` is primitive `i` of function `j`.  The only
/// invariant a constructor in this module establishes is that `ao_num`,
/// `ao_prim_num`, `ao_nucl` and `ao_power` describe the same functions;
/// `read` can return a record whose coefficient matrix is short, because the
/// reference only notices when the digest is built or the section rendered.
#[derive(Debug, Clone, PartialEq)]
pub struct AoBasis {
    /// Basis-set name stored in `ao_basis/ao_basis`.
    pub(crate) ao_basis: String,
    /// Number of Cartesian functions.
    pub(crate) ao_num: i64,
    /// Number of primitives of each Cartesian function.
    pub(crate) ao_prim_num: Vec<i64>,
    /// `max(ao_prim_num)`, or 0 when `ao_prim_num` is absent.
    pub(crate) ao_prim_num_max: i64,
    /// Nucleus index of each Cartesian function (1-based).
    pub(crate) ao_nucl: Vec<i64>,
    /// Cartesian monomial of each function (`x^a y^b z^c`).
    pub(crate) ao_power: Vec<AngmomXyz>,
    /// Contraction coefficients, `ao_prim_num_max` rows of `ao_num`.
    pub(crate) ao_coef: Vec<f64>,
    /// Exponents, same layout as `ao_coef`.
    pub(crate) ao_expo: Vec<f64>,
    /// Whether the basis uses Cartesian (`6d`, `10f`, ...) functions.
    pub(crate) ao_cartesian: bool,
}

/// The fields the AO and AO-extra handlers share.
///
/// `Input_ao_extra_basis.ml` is a near copy of `Input_ao_basis.ml`: the same
/// seven arrays, the same digest, the same ordering algorithm, only the
/// group, the key names and three extra logicals differ.  This trait lets
/// both handlers run the *same* code on those arrays instead of keeping two
/// 160-line copies of the ordering logic in sync.
pub(crate) trait Shells {
    /// Number of Cartesian functions.
    fn num(&self) -> i64;
    /// Number of primitives of each function.
    fn prim_num(&self) -> &[i64];
    /// `max(prim_num)`.
    fn prim_num_max(&self) -> i64;
    /// Cartesian monomial of each function.
    fn power(&self) -> &[AngmomXyz];
    /// Nucleus index of each function.
    fn nucl(&self) -> &[i64];
    /// Contraction coefficients, column major.
    fn coef(&self) -> &[f64];
    /// Exponents, column major.
    fn expo(&self) -> &[f64];
}

impl Shells for AoBasis {
    fn num(&self) -> i64 {
        self.ao_num
    }
    fn prim_num(&self) -> &[i64] {
        &self.ao_prim_num
    }
    fn prim_num_max(&self) -> i64 {
        self.ao_prim_num_max
    }
    fn power(&self) -> &[AngmomXyz] {
        &self.ao_power
    }
    fn nucl(&self) -> &[i64] {
        &self.ao_nucl
    }
    fn coef(&self) -> &[f64] {
        &self.ao_coef
    }
    fn expo(&self) -> &[f64] {
        &self.ao_expo
    }
}

/// The per-function vectors each handler rebuilds after [`reorder_parts`].
pub(crate) struct Reordered {
    /// Permuted `prim_num`.
    pub(crate) prim_num: Vec<i64>,
    /// Permuted `nucl`.
    pub(crate) nucl: Vec<i64>,
    /// Permuted `power`.
    pub(crate) power: Vec<AngmomXyz>,
    /// Permuted coefficient matrix.
    pub(crate) coef: Vec<f64>,
    /// Permuted exponent matrix.
    pub(crate) expo: Vec<f64>,
}

/// `read_ao_basis`: the stored name, rejecting the sentinel `"None"`.
///
/// Shared with [`super::ao_extra_basis`], which differs only in its group and
/// key names.
pub(crate) fn read_name(ezfio: &EzFio, group: &str, key: &str) -> Result<String> {
    let result = ezfio.get::<String>(group, key)?;
    if result != "None" {
        qptypes::ao_basis_name_of_string(&result)
    } else {
        bail!("No basis")
    }
}

/// `read_ao_num`: the function count, validated against `AO_number.get_max`.
///
/// `AO_number.get_max` reads `ao_basis/ao_num` even when the field being
/// read is `ao_extra_num`: the wrapper's bound is global, as in the
/// reference.
pub(crate) fn read_num(ezfio: &EzFio, group: &str, key: &str) -> Result<i64> {
    qptypes::ao_number_of_int(ezfio, ezfio.get::<i64>(group, key)?)
}

/// `read_ao_prim_num`: one count per function, empty when the key is absent.
pub(crate) fn read_prim_num(ezfio: &EzFio, group: &str, key: &str) -> Result<Vec<i64>> {
    if ezfio.has(group, key)? {
        let array = ezfio.get_array::<i64>(group, key)?;
        array
            .data()
            .iter()
            .map(|&x| qptypes::ao_prim_number_of_int(x))
            .collect()
    } else {
        Ok(Vec::new())
    }
}

/// `read_ao_prim_num_max`: `max` of `ao_prim_num`, or `AO_prim_number.of_int 0`
/// when the key is absent.
///
/// Both branches fail when `ao_prim_num` is missing: the `else` branch passes
/// 0 to `AO_prim_number.of_int`, which asserts `x > 0`.  That is why a
/// directory without `ao_prim_num` reads as "no basis".
pub(crate) fn read_prim_num_max(ezfio: &EzFio, group: &str, key: &str) -> Result<i64> {
    if ezfio.has(group, key)? {
        let array = ezfio.get_array::<i64>(group, key)?;
        let maximum = array
            .data()
            .iter()
            .fold(0i64, |x, &y| if x > y { x } else { y });
        qptypes::ao_prim_number_of_int(maximum)
    } else {
        qptypes::ao_prim_number_of_int(0)
    }
}

/// `read_ao_nucl`: the nucleus index of every function, bounded by
/// `nuclei/nucl_num`.
pub(crate) fn read_nucl(ezfio: &EzFio, group: &str, key: &str) -> Result<Vec<i64>> {
    if ezfio.has(group, key)? {
        let nmax = qptypes::nucl_number_get_max(ezfio)?;
        let array = ezfio.get_array::<i64>(group, key)?;
        array
            .data()
            .iter()
            .map(|&x| qptypes::nucl_number_of_int_max(x, nmax))
            .collect()
    } else {
        Ok(Vec::new())
    }
}

/// `read_ao_power`: rebuild `"x<a>y<b>z<c>"` from the three columns.
///
/// The reference indexes all three columns unconditionally, so a `(n,3)`
/// array whose payload is short is an error even if the missing values would
/// have been zero.  The first dimension of the stored array, not `ao_num`,
/// decides how many entries are produced.
pub(crate) fn read_power(ezfio: &EzFio, group: &str, key: &str) -> Result<Vec<AngmomXyz>> {
    if !ezfio.has(group, key)? {
        return Ok(Vec::new());
    }
    let array = ezfio.get_array::<i64>(group, key)?;
    let dim = *array
        .dims()
        .first()
        .context("Invalid_argument: ao_power has no first dimension")?;
    let data = array.data();
    let mut result = Vec::with_capacity(dim);
    for i in 0..dim {
        let mut text = String::new();
        let x = *data
            .get(i)
            .context("Invalid_argument: ao_power index out of bounds (x)")?;
        if x > 0 {
            text.push('x');
            text.push_str(&x.to_string());
        }
        let y = *data
            .get(dim + i)
            .context("Invalid_argument: ao_power index out of bounds (y)")?;
        if y > 0 {
            text.push('y');
            text.push_str(&y.to_string());
        }
        let z = *data
            .get(2 * dim + i)
            .context("Invalid_argument: ao_power index out of bounds (z)")?;
        if z > 0 {
            text.push('z');
            text.push_str(&z.to_string());
        }
        result.push(AngmomXyz::of_string(&text).map_err(|error| anyhow!("{error}"))?);
    }
    Ok(result)
}

/// `read_ao_coef`: `AO_coef.of_float` is the identity.
pub(crate) fn read_coef(ezfio: &EzFio, group: &str, key: &str) -> Result<Vec<f64>> {
    if ezfio.has(group, key)? {
        Ok(ezfio.get_array::<f64>(group, key)?.data().to_vec())
    } else {
        Ok(Vec::new())
    }
}

/// `read_ao_expo`: rejects negative exponents with `AO_expo`'s message.
pub(crate) fn read_expo(ezfio: &EzFio, group: &str, key: &str) -> Result<Vec<f64>> {
    if ezfio.has(group, key)? {
        ezfio
            .get_array::<f64>(group, key)?
            .data()
            .iter()
            .map(|&x| qptypes::ao_expo_of_float(x))
            .collect()
    } else {
        Ok(Vec::new())
    }
}

/// `read_ao_cartesian`: materialise the default when absent.
///
/// Materialisation is a write, but it is the C3 contract, not a restamp:
/// the editor may create defaulted attributes that were absent.  The value is
/// `get_default "<module>" "<key>"` — a run-time lookup in
/// `$QP_ROOT/data/ezfio_defaults/` (see [`crate::defaults`]) — and only
/// `fallback`, the `EZFIO.cfg` constant, applies when that lookup finds
/// nothing (documented lenient deviation: the reference aborts).
///
/// `module` is the defaults-file section (the module id), `group` the EZFIO
/// directory; they coincide for `ao_basis` and `ao_extra_basis` but are kept
/// separate because they need not.
pub(crate) fn read_logical(
    ezfio: &EzFio,
    module: &str,
    group: &str,
    key: &str,
    fallback: bool,
) -> Result<bool> {
    if !ezfio.has(group, key)? {
        let default = match crate::defaults::lookup(module, key) {
            Some(text) => bool_of_string(&text)?,
            None => fallback,
        };
        ezfio.set::<bool>(group, key, &default)?;
    }
    Ok(ezfio.get::<bool>(group, key)?)
}

/// `bool_of_string`: exactly `"true"` and `"false"`, as in the OCaml stdlib.
///
/// The reference feeds the defaults file through this function, so a value
/// such as `T` raises `Invalid_argument` (caught by `read`'s `with _`); it
/// must not be accepted as a permissive boolean.
fn bool_of_string(text: &str) -> Result<bool> {
    match text {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => bail!("Invalid_argument: bool_of_string"),
    }
}

/// The record fields, without the `try ... with _ -> None` wrapper.
fn read_record(ezfio: &EzFio) -> Result<AoBasis> {
    Ok(AoBasis {
        ao_basis: read_name(ezfio, GROUP, "ao_basis")?,
        ao_num: read_num(ezfio, GROUP, "ao_num")?,
        ao_prim_num: read_prim_num(ezfio, GROUP, "ao_prim_num")?,
        ao_prim_num_max: read_prim_num_max(ezfio, GROUP, "ao_prim_num")?,
        ao_nucl: read_nucl(ezfio, GROUP, "ao_nucl")?,
        ao_power: read_power(ezfio, GROUP, "ao_power")?,
        ao_coef: read_coef(ezfio, GROUP, "ao_coef")?,
        ao_expo: read_expo(ezfio, GROUP, "ao_expo")?,
        ao_cartesian: read_logical(ezfio, GROUP, GROUP, "ao_cartesian", false)?,
    })
}

/// `read`: `None` on any failure of the record or of the digest.
///
/// Unlike [`crate::handlers::nuclei_by_hand::read`], the OCaml body is one
/// big `try ... with _ -> None`, so an empty name, a missing key or a
/// negative exponent all turn the section into "absent" instead of aborting
/// the editor.  `to_md5` is part of the `try` in the reference; here it is
/// computed but **not written** (decision D4).
pub fn read(ezfio: &EzFio) -> Result<Option<AoBasis>> {
    let Ok(record) = read_record(ezfio) else {
        return Ok(None);
    };
    match to_md5(&record) {
        Ok(_) => Ok(Some(record)),
        Err(_) => Ok(None),
    }
}

/// `to_long_basis`: one `(monomial, shell, nucleus)` per Cartesian function.
///
/// `ao_expo`/`ao_coef` are read as `ao_coef[ao_num*j + i]`, i.e. primitive
/// `j` of function `i`.
///
/// # Errors
///
/// The failures the reference can raise here, all swallowed by [`read`]:
/// an index out of bounds, `Gto.of_prim_coef_list` on an empty primitive
/// list (`Failure "hd"`) or on mixed symmetries, and the `assert false` when
/// `ao_power`, `ao_nucl` and the built shells do not have the same length.
pub(crate) fn to_long_basis(b: &impl Shells) -> Result<LongBasis> {
    let ao_num = usize_of(b.num())?;
    let mut shells = Vec::with_capacity(ao_num);
    for i in 0..ao_num {
        let xyz = *b
            .power()
            .get(i)
            .context("Invalid_argument: ao_power index out of bounds")?;
        let sym = xyz.to_symmetry().map_err(gto_error)?;
        let prim_num = usize_of(
            *b.prim_num()
                .get(i)
                .context("Invalid_argument: ao_prim_num index out of bounds")?,
        )?;
        let mut primitives = Vec::with_capacity(prim_num);
        for j in 0..prim_num {
            let index = row_index(ao_num, j, i)?;
            let expo = *b
                .expo()
                .get(index)
                .context("Invalid_argument: ao_expo index out of bounds")?;
            let coef = *b
                .coef()
                .get(index)
                .context("Invalid_argument: ao_coef index out of bounds")?;
            primitives.push((GaussianPrimitive::of_sym_expo(sym, expo), coef));
        }
        shells.push(Gto::of_prim_coef_list(primitives).map_err(gto_error)?);
    }
    if b.power().len() != ao_num || b.nucl().len() != ao_num {
        return Err(assertion("to_long_basis: assertion false failed"));
    }
    let mut result = Vec::with_capacity(ao_num);
    for ((xyz, gto), &nucl) in b.power().iter().copied().zip(shells).zip(b.nucl()) {
        // `Long_basis` stores `NuclNumber`, a 32-bit index; the EZFIO file
        // stores 64-bit integers, so a corrupt value can overflow here.
        let nucl =
            i32::try_from(nucl).with_context(|| format!("Nucl_number out of range: {nucl}"))?;
        result.push((xyz, gto, nucl));
    }
    Ok(LongBasis(result))
}

/// `to_basis`: collapse [`to_long_basis`] back to one entry per shell.
pub(crate) fn to_basis(b: &impl Shells) -> Result<basis::Basis> {
    Ok(to_long_basis(b)?.to_basis())
}

/// `to_md5`: `Basis.to_md5` of the collapsed basis.
pub(crate) fn to_md5(b: &impl Shells) -> Result<String> {
    Ok(basis::to_md5(&to_basis(b)?))
}

/// `ordering`: for each Cartesian function of the canonical basis, its index
/// in the stored basis.
///
/// The first pass matches the full `(monomial, Gto, nucleus)` triple.  The
/// second pass, used when the exact match fails, accepts any function with
/// the same monomial and nucleus whose primitives (symmetry and exponent)
/// agree, ignoring the contraction coefficients.  Matches are consumed, so
/// two identical shells map to two different indices; `-1` marks a function
/// of the canonical basis that is not in the stored one.
pub(crate) fn ordering(b: &impl Shells) -> Result<Vec<i64>> {
    let ordered = LongBasis::of_basis(&to_basis(b)?);
    let unordered = to_long_basis(b)?;
    let mut search: Vec<Option<usize>> = (0..unordered.0.len()).map(Some).collect();
    let mut result = Vec::with_capacity(ordered.0.len());
    for x in &ordered.0 {
        result.push(find(x, &unordered.0, &mut search));
    }
    Ok(result)
}

/// The `find`/`find2` pair of the reference, over a consumable index table.
fn find(
    x: &(AngmomXyz, Gto, i32),
    unordered: &[(AngmomXyz, Gto, i32)],
    search: &mut [Option<usize>],
) -> i64 {
    for (i, entry) in unordered.iter().enumerate() {
        if search[i].is_some() && entry == x {
            search[i] = None;
            return i as i64;
        }
    }
    for (i, entry) in unordered.iter().enumerate() {
        let Some(_) = search[i] else { continue };
        if x.0 != entry.0 || x.2 != entry.2 {
            continue;
        }
        if same_primitives(&x.1, &entry.1) {
            search[i] = None;
            return i as i64;
        }
    }
    -1
}

/// Equality of the primitive lists, ignoring contraction coefficients.
fn same_primitives(a: &Gto, b: &Gto) -> bool {
    a.lc.len() == b.lc.len()
        && a.lc
            .iter()
            .zip(&b.lc)
            .all(|((p, _), (q, _))| p.sym == q.sym && p.expo == q.expo)
}

/// `reorder`: sort the record into the canonical order of [`ordering`].
///
/// The per-function vectors are permuted; the coefficient and exponent
/// matrices are permuted **row by row** (each row is one primitive index,
/// `ao_num` entries), so the result keeps the same `ao_num × prim_num_max`
/// layout.
///
/// # Errors
///
/// `"index out of bounds"` when `ordering` contains `-1` (a function of the
/// canonical basis that the stored basis does not have), which is the OCaml
/// `Invalid_argument` from `a.(order.(i))`.
pub(crate) fn reorder(b: &AoBasis) -> Result<AoBasis> {
    let parts = reorder_parts(b)?;
    Ok(AoBasis {
        ao_prim_num: parts.prim_num,
        ao_nucl: parts.nucl,
        ao_power: parts.power,
        ao_coef: parts.coef,
        ao_expo: parts.expo,
        ..b.clone()
    })
}

/// The body of `reorder`, shared with the AO-extra handler.
pub(crate) fn reorder_parts(b: &impl Shells) -> Result<Reordered> {
    let order = ordering(b)?;
    let ao_num = usize_of(b.num())?;
    let prim_num_max = usize_of(b.prim_num_max())?;
    Ok(Reordered {
        prim_num: permute(b.prim_num(), &order)?,
        nucl: permute(b.nucl(), &order)?,
        power: permute(b.power(), &order)?,
        coef: reorder_matrix(b.coef(), ao_num, prim_num_max, &order)?,
        expo: reorder_matrix(b.expo(), ao_num, prim_num_max, &order)?,
    })
}

/// `f a = Array.init (Array.length a) (fun i -> a.(order.(i)))`.
pub(crate) fn permute<T: Clone>(values: &[T], order: &[i64]) -> Result<Vec<T>> {
    order
        .iter()
        .map(|&index| {
            let index = usize::try_from(index).map_err(|_| anyhow!("index out of bounds"))?;
            values
                .get(index)
                .cloned()
                .ok_or_else(|| anyhow!("index out of bounds"))
        })
        .collect()
}

/// The `Array.init prim_num_max (fun i -> f (row i)) |> concat` pattern.
pub(crate) fn reorder_matrix(
    values: &[f64],
    ao_num: usize,
    prim_num_max: usize,
    order: &[i64],
) -> Result<Vec<f64>> {
    let mut result = Vec::with_capacity(values.len());
    for i in 0..prim_num_max {
        let mut row = Vec::with_capacity(ao_num);
        for j in 0..ao_num {
            let index = row_index(ao_num, i, j)?;
            row.push(
                *values
                    .get(index)
                    .context("Invalid_argument: coefficient index out of bounds")?,
            );
        }
        result.extend(permute(&row, order)?);
    }
    Ok(result)
}

/// `write`: the nine fields, in record order, MD5 first.
///
/// This is a write path, so the MD5 stamp is written here (decision D4 only
/// forbids it in `read`).  The arrays are built with
/// [`super::ezfio_array_of_list`], which truncates extra data and raises
/// `"Array is not consistent"` on short data, exactly like the reference.
pub(crate) fn write(ezfio: &EzFio, b: &AoBasis) -> Result<()> {
    let ao_num = usize_of(b.ao_num)?;
    let prim_num_max = usize_of(b.ao_prim_num_max)?;

    ezfio.set::<String>(GROUP, "ao_md5", &to_md5(b)?)?;
    ezfio.set::<String>(GROUP, "ao_basis", &b.ao_basis)?;

    let prim_num = super::ezfio_array_of_list(vec![ao_num], b.ao_prim_num.clone())?;
    ezfio.set_array::<i64>(GROUP, "ao_prim_num", &prim_num)?;

    let nucl = super::ezfio_array_of_list(vec![ao_num], b.ao_nucl.clone())?;
    ezfio.set_array::<i64>(GROUP, "ao_nucl", &nucl)?;

    let mut power = Vec::with_capacity(3 * b.ao_power.len());
    power.extend(b.ao_power.iter().map(|xyz| i64::from(xyz.x)));
    power.extend(b.ao_power.iter().map(|xyz| i64::from(xyz.y)));
    power.extend(b.ao_power.iter().map(|xyz| i64::from(xyz.z)));
    let power = super::ezfio_array_of_list(vec![ao_num, 3], power)?;
    ezfio.set_array::<i64>(GROUP, "ao_power", &power)?;

    ezfio.set::<bool>(GROUP, "ao_cartesian", &b.ao_cartesian)?;

    let coef = super::ezfio_array_of_list(vec![ao_num, prim_num_max], b.ao_coef.clone())?;
    ezfio.set_array::<f64>(GROUP, "ao_coef", &coef)?;

    let expo = super::ezfio_array_of_list(vec![ao_num, prim_num_max], b.ao_expo.clone())?;
    ezfio.set_array::<f64>(GROUP, "ao_expo", &expo)?;
    Ok(())
}

/// `to_rst`: the editable block, without the section header.
///
/// `Basis.to_string` renders the collapsed basis (one entry per shell), and
/// every line is indented by two spaces.  The trailing table is one line per
/// **stored Cartesian function**, numbered from 1, with the nucleus index and
/// the shell symbol followed by the monomial (`s`, `p(x)`, `d(xy)`, …) in a
/// 9-column left-justified field.
pub fn to_rst(b: &AoBasis) -> Result<String> {
    let short_basis = to_basis(b)?;
    let basis_text = basis::to_string(Fmt::Gamess, &short_basis, None);
    let basis_text = basis_text
        .split('\n')
        .map(|line| format!("  {line}"))
        .collect::<Vec<_>>()
        .join("\n");

    let mut print_sym = String::new();
    for (i, xyz) in b.ao_power.iter().enumerate() {
        let nucl = b
            .ao_nucl
            .get(i)
            .context("Invalid_argument: ao_nucl index out of bounds")?;
        // `AngmomXyz`'s `Display` ignores the formatter width (it writes
        // through `write_str`), so the symbol string is padded, not the value.
        let xyz = xyz.to_symbol_monomial();
        print_sym.push_str(&format!(" {:>5}  {:>6}     {:<9}\n", i + 1, nucl, xyz));
    }

    Ok(format!(
        "\nName of the AO basis ::\n\n  ao_basis = {}\n\n\
         Cartesian coordinates (6d,10f,...) ::\n\n  ao_cartesian = {}\n\n\
         Basis set (read-only) ::\n\n{}\n\n\n\
         ======= ========= ===================\n Basis   Nucleus   Cartesian function\n\
         ======= ========= ===================\n{}\n\
         ======= ========= ===================\n\n",
        b.ao_basis, b.ao_cartesian, basis_text, print_sym
    ))
}

/// `get` for this section: `None` renders as `""`.
pub(crate) fn get(ezfio: &EzFio, title: &str) -> Result<String> {
    match read(ezfio)? {
        Some(record) => Ok(make_header(title) + &to_rst(&record)?),
        None => Ok(String::new()),
    }
}

/// `ao_num * j + i` in the flattened coefficient matrix.
pub(crate) fn row_index(ao_num: usize, j: usize, i: usize) -> Result<usize> {
    ao_num
        .checked_mul(j)
        .and_then(|base| base.checked_add(i))
        .ok_or_else(|| anyhow!("index out of bounds"))
}

/// The two `Gto.of_prim_coef_list` failures with the OCaml exception payloads.
pub(crate) fn gto_error(error: GtoError) -> anyhow::Error {
    match error {
        // `Failure "hd"`: the reference's `Gto.ml` matches on the empty list.
        GtoError::EmptyPrimitiveList => anyhow!("hd"),
        // `MixedSymmetry` already displays `Failure in of_prim_coef_list`.
        other => anyhow!("{other}"),
    }
}

/// A negative count is a corrupt file, not a panic.
pub(crate) fn usize_of(value: i64) -> Result<usize> {
    usize::try_from(value).map_err(|_| anyhow!("negative dimension or index: {value}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::handlers::test_support::hcn_ezfio;
    use crate::ocaml_compat::is_assert_failure;

    /// The `hcn` fixture: the digest is the one the reference binary wrote.
    #[test]
    fn get_renders_the_reference_section_byte_for_byte() {
        let (_dir, ezfio) = hcn_ezfio();
        let rendered = get(&ezfio, "AO basis").expect("get");
        assert_eq!(
            qp_types::md5::md5_hex(&rendered),
            "4403d163b97d14d1e252b860e42c5f77",
            "the AO basis section no longer matches the reference"
        );
    }

    /// D4: `read` computes the digest (and therefore fails like the
    /// reference) but never writes it.
    #[test]
    fn read_does_not_restamp_the_md5() {
        let (_dir, ezfio) = hcn_ezfio();
        let before = ezfio.get::<String>(GROUP, "ao_md5").expect("stored digest");
        assert!(read(&ezfio).expect("read").is_some());
        let after = ezfio.get::<String>(GROUP, "ao_md5").expect("stored digest");
        assert_eq!(before, after);
        assert_eq!(
            to_md5(&read(&ezfio).expect("read").expect("some")).expect("md5"),
            before
        );
    }

    #[test]
    fn read_returns_none_without_a_name_and_writes_nothing() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        assert!(read(&ezfio).expect("read").is_none());
        assert!(!ezfio.has(GROUP, "ao_md5").expect("has"));
        assert!(!ezfio.has(GROUP, "ao_basis").expect("has"));
    }

    /// The name check happens first: `"None"` is the "no basis" sentinel.
    #[test]
    fn read_returns_none_for_the_none_sentinel() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        ezfio
            .set::<String>(GROUP, "ao_basis", &"None".to_owned())
            .expect("set");
        assert!(read(&ezfio).expect("read").is_none());
    }

    /// The defaults file is parsed with OCaml's strict `bool_of_string`: only
    /// `true` and `false` are logicals, so the `T`/`F` spelling of the EZFIO
    /// files is *not* a valid default.  The error is the one `read`'s
    /// `try ... with _ -> None` swallows.
    #[test]
    fn bool_of_string_is_the_strict_ocaml_predicate() {
        assert!(bool_of_string("true").expect("true"));
        assert!(!bool_of_string("false").expect("false"));
        for text in ["T", "F", "True", "False", "0", ""] {
            let error = bool_of_string(text).expect_err(text);
            assert_eq!(error.to_string(), "Invalid_argument: bool_of_string");
        }
    }

    /// A P shell stored as `z, x, y`: `ordering` finds the canonical order
    /// and `reorder` permutes the three rows.
    #[test]
    fn ordering_and_reorder_fix_a_misordered_shell() {
        let record = AoBasis {
            ao_basis: "test".to_owned(),
            ao_num: 3,
            ao_prim_num: vec![1, 1, 1],
            ao_prim_num_max: 1,
            ao_nucl: vec![1, 1, 1],
            ao_power: vec![
                AngmomXyz::new(0, 0, 1),
                AngmomXyz::new(1, 0, 0),
                AngmomXyz::new(0, 1, 0),
            ],
            // The three components of a real P shell share the same
            // contraction: only the monomial differs.
            ao_coef: vec![1.0, 1.0, 1.0],
            ao_expo: vec![2.0, 2.0, 2.0],
            ao_cartesian: false,
        };
        assert_eq!(ordering(&record).expect("ordering"), vec![1, 2, 0]);
        let reordered = reorder(&record).expect("reorder");
        assert_eq!(
            reordered.ao_power,
            vec![
                AngmomXyz::new(1, 0, 0),
                AngmomXyz::new(0, 1, 0),
                AngmomXyz::new(0, 0, 1)
            ]
        );
        // Widths are unchanged because the components share the contraction;
        // `reorder_matrix` is exercised on distinct values below.
        assert_eq!(reordered.ao_coef, vec![1.0, 1.0, 1.0]);
        assert_eq!(reordered.ao_expo, vec![2.0, 2.0, 2.0]);
        assert_eq!(ordering(&reordered).expect("ordering"), vec![0, 1, 2]);
        // The digest depends on the shell, not on the component order.
        assert_eq!(
            to_md5(&record).expect("md5"),
            to_md5(&reordered).expect("md5")
        );
    }

    #[test]
    fn reorder_rejects_a_function_that_is_not_in_the_record() {
        // A P shell with only its `x` component: the canonical basis has
        // `x`, `y` and `z`, and only the first one is in the record.
        let record = AoBasis {
            ao_basis: "test".to_owned(),
            ao_num: 1,
            ao_prim_num: vec![1],
            ao_prim_num_max: 1,
            ao_nucl: vec![1],
            ao_power: vec![AngmomXyz::new(1, 0, 0)],
            ao_coef: vec![1.0],
            ao_expo: vec![2.0],
            ao_cartesian: false,
        };
        assert_eq!(ordering(&record).expect("ordering"), vec![0, -1, -1]);
        let error = reorder(&record).expect_err("no match");
        assert_eq!(error.to_string(), "index out of bounds");
    }

    /// The coefficient matrix is `prim_num_max` rows of `ao_num`, and each
    /// row is permuted independently: values from different rows must not
    /// migrate.
    #[test]
    fn reorder_permutes_each_primitive_row_independently() {
        let values = [10.0, 11.0, 20.0, 21.0];
        let permuted = reorder_matrix(&values, 2, 2, &[1, 0]).expect("permute");
        assert_eq!(permuted, vec![11.0, 10.0, 21.0, 20.0]);
        let error = reorder_matrix(&values, 2, 2, &[0, 2]).expect_err("out of bounds");
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

    #[test]
    fn write_reports_short_arrays_with_the_ocaml_message() {
        let (_dir, ezfio) = hcn_ezfio();
        let mut record = read(&ezfio).expect("read").expect("some");
        record.ao_coef.pop();
        let error = write(&ezfio, &record).expect_err("short");
        assert_eq!(error.to_string(), "Array is not consistent");
    }

    /// `to_long_basis` on a record whose vectors disagree hits the OCaml
    /// `assert false` (an `Assert_failure`, exit 2 if it escaped).
    #[test]
    fn to_long_basis_reports_the_length_assertion() {
        let record = AoBasis {
            ao_basis: "test".to_owned(),
            ao_num: 1,
            ao_prim_num: vec![1],
            ao_prim_num_max: 1,
            ao_nucl: vec![1],
            ao_power: vec![AngmomXyz::new(0, 0, 0)],
            ao_coef: vec![1.0],
            ao_expo: vec![1.0],
            ao_cartesian: false,
        };
        let mut longer = record.clone();
        longer.ao_power.push(AngmomXyz::new(0, 0, 0));
        let error = to_long_basis(&longer).expect_err("mismatch");
        assert!(is_assert_failure(&error), "{error}");
    }
}
