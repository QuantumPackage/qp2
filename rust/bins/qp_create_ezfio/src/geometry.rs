//! Geometry input: the XYZ reader, its Z-matrix fallback, and ghost centres.
//!
//! `ocaml/qp_create_ezfio.ml` reads the molecule with
//!
//! ```text
//! try Molecule.of_file file ~charge ~multiplicity ~units
//! with Failure msg -> failwith ("cannot parse geometry from %s: ... (parse error: %s)" ...)
//!    | Sys_error msg -> failwith ("cannot parse geometry from %s: ... (%s)" ...)
//! ```
//!
//! and `Molecule.of_file` is `try of_xyz_file with XYZError -> of_zmt_file`.
//! The two XYZ failures are **not** interchangeable and this module keeps them
//! apart:
//!
//! * a file with fewer than two lines, a non-numeric or non-positive atom
//!   count is `XYZError` → caught → the *same file* is re-parsed as a
//!   Z-matrix (that parse then fails, e.g. with `Element … unknown`);
//! * fewer atom lines than the announced count is the `Array.sub` in the XYZ
//!   reader, an `Invalid_argument`, which `of_file` does not catch. It must
//!   escape to the runtime as
//!   `Fatal error: exception Invalid_argument("Array.sub")`.
//!
//! `qp_types::molecule::of_xyz_file_with` maps the second case to its own
//! `XyzError`, so the file is read here with `input_line` semantics (a
//! trailing newline does not create an empty last line) and only the
//! extracted atom lines are handed to [`Atom::of_string`].
//!
//! # Units
//!
//! Positions are stored in **bohr**. The file is in angstrom unless `au` is
//! set. Coordinates are converted to bohr by [`Point3d::of_string`]/
//! [`Point3d::of_tuple`], which **divide** by `BOHR_TO_ANGSTROM`
//! (`0.52917721067121`) — reference commit `15e725aa` — not by multiplying
//! by the reciprocal. (The one exception is the covalent radius in
//! [`covalent_radius`], which still multiplies by `ANGSTROM_TO_BOHR`;
//! `15e725aa` did not change radii.) Z-matrix distances follow the same
//! rule as coordinates: they are angstrom unless `au` is set.
//!
//! # Ghost centres
//!
//! For every pair `i > j` closer than `(r_i + r_j) * threshold`, where `r` is
//! the covalent radius in bohr, one ghost `Element::X` of charge 0 is placed
//! at `x_i + (x_j - x_i) * r_i / (r_i + r_j)` — weighted by the radii, not
//! the midpoint. Ghosts are appended after the nuclei, in `i` then `j`
//! ascending order. A zero or negative threshold disables them.

use std::path::Path;

use qp_types::molecule::{Molecule, MultiplicityError};
use qp_types::point3d::Units;
use qp_types::point3d::units::ANGSTROM_TO_BOHR;
use qp_types::zmatrix::{Zmatrix, ZmatrixError};
use qp_types::{Atom, Charge, Element, ElementError, Point3d, element};

use crate::error::{AppError, OcamlException};
use crate::ocaml_num;

/// The fixed middle of every geometry diagnostic (`Molecule.of_file`'s
/// caller prefixes it with the file name and appends the detail).
const EXPECTED: &str =
    "expected xyz (atom count on line 1, then 'Element x y z' lines) or a Gaussian z-matrix";

/// The geometry as the reference's `run` holds it after ghost insertion.
#[derive(Debug, Clone)]
pub struct Geometry {
    /// The molecule as read: real nuclei only, with the electron counts
    /// derived from them. [`Molecule::name`] and the electron keys use this.
    pub molecule: Molecule,
    /// `molecule.nuclei @ ghost`: the list every later stage iterates.
    /// Real nuclei keep their order; ghost centres follow, `Element::X` and
    /// charge 0.
    pub nuclei: Vec<Atom>,
}

/// Read a geometry the way `Molecule.of_file` does, then add ghost centres.
///
/// * `input` — the `FILE` argument, used verbatim in diagnostics.
/// * `au` — `-a`: the coordinates and Z-matrix distances are already bohr.
/// * `charge` — `-c`, as the total molecular charge.
/// * `multiplicity` — `-m`, `2S+1`. Checked **before** the file is opened,
///   because the reference constructs the `Multiplicity` while converting
///   the options; `-m 0` therefore fails even when `input` does not exist.
/// * `ghost_threshold` — `-g`, a multiplier of the covalent radii; `0.0`
///   and negative values produce no ghosts.
///
/// # Errors
///
/// * [`AppError::Uncaught`] for the exceptions the reference does not catch:
///   `Invalid_argument("Strictly_positive_int : (x > 0) : x=…")`,
///   `Invalid_argument("Array.sub")`, `Element.ElementError`,
///   `Atom.AtomError`, the `Qptypes.ml` assertion sites and
///   `Molecule.MultiplicityError`.
/// * [`AppError::Prefixed`] for the two caught classes, wrapped as
///   `cannot parse geometry from <input>: … (parse error: <msg>)` and
///   `… (<input>: <errno text>)`.
pub fn read(
    input: &Path,
    au: bool,
    charge: i64,
    multiplicity: i64,
    ghost_threshold: f64,
) -> Result<Geometry, AppError> {
    let multiplicity = check_multiplicity(multiplicity)?;
    // `Charge.t` is a transparent float in OCaml; `Charge::of_int` takes an
    // `i32`, so go through the float directly rather than truncating.
    let charge = Charge(charge as f64);
    let units = if au { Units::Bohr } else { Units::Angstrom };

    let raw = read_text(input)?;
    // `Molecule.ml` reads with `input_line`, then removes every `\r` and
    // turns every tab into a space before parsing. Doing it on the whole
    // text is equivalent and keeps the Z-matrix path on the same bytes.
    let text = raw.replace('\r', "").replace('\t', " ");
    let molecule = parse_molecule(&text, input, charge, multiplicity, units)?;
    let ghost = ghost_centers(&molecule, ghost_threshold)?;
    let mut nuclei = molecule.nuclei.clone();
    nuclei.extend(ghost);
    Ok(Geometry { molecule, nuclei })
}

/// `Multiplicity.of_int` / `Strictly_positive_int.of_int`: `x >= 1`.
///
/// The reference stores the value in a 63-bit `int`; this returns it
/// unchanged so the later electron arithmetic cannot overflow.
fn check_multiplicity(multiplicity: i64) -> Result<i64, AppError> {
    if multiplicity >= 1 {
        Ok(multiplicity)
    } else {
        // Measured: `qp_create_ezfio -b x -m 0 mol.xyz` →
        // `Fatal error: exception Invalid_argument("Strictly_positive_int :
        // (x > 0) : x=0")`.
        Err(AppError::uncaught(OcamlException::InvalidArgument(
            format!("Strictly_positive_int : (x > 0) : x={multiplicity}"),
        )))
    }
}

/// Read the file as text, turning an I/O failure into the prefixed
/// `Sys_error` form of the geometry wrapper.
///
/// The reference calls `open_in`, so the detail is `<path>: <strerror>`.
/// Only `ENOENT` is measured; the other strings in [`strerror`] use the C
/// locale's wording.
fn read_text(input: &Path) -> Result<String, AppError> {
    std::fs::read_to_string(input).map_err(|error| {
        let detail = match error.raw_os_error() {
            Some(code) => format!("{}: {}", input.display(), strerror(code)),
            // No errno: the only such failure `read_to_string` reports is
            // invalid UTF-8, which the reference would have read as bytes.
            None => format!("{}: {error}", input.display()),
        };
        wrapper(input, &detail)
    })
}

/// The C library's `strerror` text for the errnos a geometry read can hit.
///
/// The reference prints `Sys_error`'s payload directly; Rust keeps the
/// message behind `Display`, which appends `(os error N)`. This table keeps
/// the four messages a user is likely to see exact.
///
/// `pub(crate)` because [`crate::basis_read`] needs the same wording for the
/// `Sys_error` an unreadable basis channel raises.
pub(crate) fn strerror(code: i32) -> String {
    match code {
        1 => "Operation not permitted".to_owned(),
        2 => "No such file or directory".to_owned(),
        13 => "Permission denied".to_owned(),
        20 => "Not a directory".to_owned(),
        21 => "Is a directory".to_owned(),
        40 => "Too many levels of symbolic links".to_owned(),
        other => format!("os error {other}"),
    }
}

/// The prefixed wrapper around a caught `Failure`/`Sys_error`.
fn wrapper(input: &Path, detail: &str) -> AppError {
    AppError::fail(format!(
        "cannot parse geometry from {}: {EXPECTED} ({detail})",
        input.display()
    ))
}

/// Why the XYZ parse did not produce a molecule.
enum XyzFailure {
    /// The reference's `XYZError`: fall back to the Z-matrix reader.
    NotXyz,
    /// Anything the reference would not have caught, already classified.
    Fatal(AppError),
}

/// `Molecule.of_file`: XYZ first, Z-matrix on `XYZError` only.
fn parse_molecule(
    text: &str,
    input: &Path,
    charge: Charge,
    multiplicity: i64,
    units: Units,
) -> Result<Molecule, AppError> {
    match parse_xyz(text, input, charge, multiplicity, units) {
        Ok(molecule) => Ok(molecule),
        Err(XyzFailure::Fatal(error)) => Err(error),
        Err(XyzFailure::NotXyz) => {
            // `of_zmt_file`: whole file → `Zmatrix.of_string` →
            // `to_xyz_string` → `of_xyz_string`. `Zmatrix::of_string` drops
            // blank lines and trims, so the trailing-newline question does
            // not reach it.
            let zmatrix = Zmatrix::of_string(text).map_err(|error| zmatrix_error(input, error))?;
            let xyz = zmatrix
                .to_xyz_string()
                .map_err(|error| zmatrix_error(input, error))?;
            parse_atoms(&xyz, input, charge, multiplicity, units)
        }
    }
}

/// `Molecule.of_xyz_file`: the `natoms :: title :: rest` match.
///
/// `input` is only used by [`atom_error`]; the atom lines themselves carry
/// no file name.
fn parse_xyz(
    text: &str,
    input: &Path,
    charge: Charge,
    multiplicity: i64,
    units: Units,
) -> Result<Molecule, XyzFailure> {
    let lines = input_lines(text);
    if lines.len() < 2 {
        return Err(XyzFailure::NotXyz);
    }
    let Ok(natoms) = ocaml_num::int_of_string(lines[0]) else {
        return Err(XyzFailure::NotXyz);
    };
    if natoms <= 0 {
        return Err(XyzFailure::NotXyz);
    }
    let rest = &lines[2..];
    // `Array.sub a 0 natoms` with `natoms > Array.length a` is an
    // `Invalid_argument`, not an `XYZError`; `of_file` does not catch it.
    if i64::try_from(rest.len()).unwrap_or(i64::MAX) < natoms {
        return Err(XyzFailure::Fatal(AppError::uncaught(
            OcamlException::InvalidArgument("Array.sub".to_owned()),
        )));
    }
    let natoms = usize::try_from(natoms).unwrap_or(usize::MAX);
    let body = rest[..natoms].join("\n");
    parse_atoms(&body, input, charge, multiplicity, units).map_err(XyzFailure::Fatal)
}

/// Replicate `input_line`: split on `\n` and drop the empty tail a final
/// newline leaves behind.
///
/// The distinction matters: `3\nt\nH\nH\n` has **two** atom lines for the
/// reference (whose `input_line` never returns the empty tail), so it is the
/// `Array.sub` crash. A plain `split('\n')` would see a third, empty line
/// and silently parse two atoms instead.
fn input_lines(text: &str) -> Vec<&str> {
    let mut lines: Vec<&str> = text.split('\n').collect();
    if text.ends_with('\n') {
        lines.pop();
    }
    lines
}

/// `Molecule.of_xyz_string`: one atom per non-empty line, then the electron
/// bookkeeping of `Qptypes` and `Molecule.ml:25-88`.
///
/// The four failure payloads are the measured ones:
///
/// * `n_elec <= 0` → `Qptypes.ml:57:20-26` assertion (`Elec_number (x > 0)`);
/// * `elec_beta < 0` → `Qptypes.ml:70:20-26` (`Elec_beta (x >= 0)`) — this
///   fires before the alpha assert because OCaml evaluates the tuple
///   right-to-left;
/// * `elec_alpha < 1` → `Qptypes.ml:83:20-26` (`Elec_alpha (x > 0)`);
/// * the implied multiplicity differs → `Molecule.MultiplicityError`.
fn parse_atoms(
    text: &str,
    input: &Path,
    charge: Charge,
    multiplicity: i64,
    units: Units,
) -> Result<Molecule, AppError> {
    let mut nuclei = Vec::new();
    for line in text.split('\n') {
        if line.is_empty() {
            continue;
        }
        // `Atom.of_string` matches a four- or five-token line and raises
        // `AtomError s` on anything else; check the shape first so the
        // payload is the reference's untrimmed line.
        let tokens: Vec<&str> = line.split(' ').filter(|token| !token.is_empty()).collect();
        if tokens.len() != 4 && tokens.len() != 5 {
            return Err(AppError::uncaught(OcamlException::AtomError(
                line.to_owned(),
            )));
        }
        nuclei.push(Atom::of_string(units, line).map_err(|error| atom_error(&error, input))?);
    }

    // `Molecule.of_xyz_string`: ne = trunc(Σ Z − 1) + 1 − trunc(charge),
    // summed left to right over the nuclei.
    let mut sum = 0.0_f64;
    for atom in &nuclei {
        sum += atom.charge.to_float();
    }
    let ne = Charge::of_float(sum - 1.0).to_int() + 1 - charge.to_int();
    if ne <= 0 {
        return Err(AppError::uncaught(OcamlException::Assert {
            file: "Qptypes.ml".to_owned(),
            line: 57,
            start: 20,
            end: 26,
        }));
    }
    // `Multiplicity.to_alpha_beta` in 64-bit arithmetic: the reference's
    // `int` is 63-bit, so `-m` values past `i32` must not wrap here.
    let beta = (i64::from(ne) - multiplicity + 1) / 2;
    if beta < 0 {
        return Err(AppError::uncaught(OcamlException::Assert {
            file: "Qptypes.ml".to_owned(),
            line: 70,
            start: 20,
            end: 26,
        }));
    }
    let beta = i32::try_from(beta).expect("invariant violated: 0 <= beta <= n_elec");
    let alpha = ne - beta;
    if alpha < 1 {
        return Err(AppError::uncaught(OcamlException::Assert {
            file: "Qptypes.ml".to_owned(),
            line: 83,
            start: 20,
            end: 26,
        }));
    }
    if 1 + alpha - beta != i32::try_from(multiplicity).unwrap_or(i32::MAX) {
        return Err(AppError::uncaught(OcamlException::MultiplicityError(
            MultiplicityError {
                electrons: ne,
                multiplicity: i32::try_from(multiplicity).unwrap_or(i32::MAX),
            }
            .to_string(),
        )));
    }
    Ok(Molecule {
        nuclei,
        elec_alpha: alpha,
        elec_beta: beta,
    })
}

/// Classify an `Atom::of_string` failure.
///
/// An unknown element stays an unprefixed `Element.ElementError`; a
/// coordinate or charge that does not parse is OCaml's
/// `Failure("float_of_string")` under the wrapper's `parse error:` detail.
/// The reference additionally prints `Bad string: …` / `name = …`
/// diagnostics from `Point3d.of_string`/`Atom.of_string`, which the crate's
/// parser does not; that is a known deviation (see `NOTES.md`).
fn atom_error(error: &anyhow::Error, input: &Path) -> AppError {
    if let Some(element) = error.downcast_ref::<ElementError>() {
        return AppError::uncaught(OcamlException::ElementError(element.to_string()));
    }
    if error
        .chain()
        .any(|cause| cause.downcast_ref::<std::num::ParseFloatError>().is_some())
    {
        return wrapper(input, "parse error: float_of_string");
    }
    let message = error
        .chain()
        .last()
        .map_or_else(|| error.to_string(), ToString::to_string);
    wrapper(input, &format!("parse error: {message}"))
}

/// Map a Z-matrix failure to the exception the reference raises.
///
/// Measured: `Q.zmt` with an unknown element prints
/// `Fatal error: exception Element.ElementError("Element Q unknown")`, and a
/// malformed line prints the prefixed
/// `… or a Gaussian z-matrix (parse error: Syntax error: H 0 0 0)`.
fn zmatrix_error(input: &Path, error: ZmatrixError) -> AppError {
    match error {
        ZmatrixError::Element(element) => {
            AppError::uncaught(OcamlException::ElementError(element.to_string()))
        }
        // `assert (x > 0)` in `atom_id_of_int`, `Zmatrix.ml:22:13-19`.
        ZmatrixError::NonPositiveAtomId { .. } => AppError::uncaught(OcamlException::Assert {
            file: "Zmatrix.ml".to_owned(),
            line: 22,
            start: 13,
            end: 19,
        }),
        // The reference indexes past the end of its array and aborts with
        // OCaml's own `Invalid_argument("index out of bounds")`.
        ZmatrixError::AtomOutOfBounds { .. } => AppError::uncaught(
            OcamlException::InvalidArgument("index out of bounds".to_owned()),
        ),
        // `Failure` messages, which `run`'s first handler catches. The
        // `NegativeDistance` arm is unreachable through `Zmatrix::of_string`
        // (the string parser turns a negative number into a label).
        other => wrapper(input, &format!("parse error: {other}")),
    }
}

/// `Element.covalent_radius` in **bohr**.
///
/// The crate returns angstrom; the reference's radii are bohr, so this
/// multiplies by `ANGSTROM_TO_BOHR` (it does not divide). Elements from
/// berkelium on have no radius and the reference raises
/// `Element.ElementError` — unprefixed, because the ghost pass is outside
/// the geometry `try`.
fn covalent_radius(element: Element) -> Result<f64, AppError> {
    element::covalent_radius(element)
        .map(|angstrom| angstrom * ANGSTROM_TO_BOHR)
        .map_err(|error| AppError::uncaught(OcamlException::ElementError(error.to_string())))
}

/// `ghost_centers`: bond-centred ghost atoms.
///
/// The reference walks `(i, j)` with `i` from `n-1` down to 1 and `j` from
/// `i-1` down to 0, prepending each match; the result is therefore ordered
/// by `i` then `j` ascending, which is what this function returns after the
/// explicit reverse.
///
/// # Cost
///
/// `O(n²)` time and — through [`Molecule::distance_matrix`] — `O(n²)`
/// memory, exactly like the reference. The radius lookups happen inside the
/// pair loop, so a nucleus without a covalent radius reports the highest
/// `(i, j)` pair first, as the reference does.
fn ghost_centers(molecule: &Molecule, threshold: f64) -> Result<Vec<Atom>, AppError> {
    let distances = molecule.distance_matrix();
    let nuclei = &molecule.nuclei;
    let n = nuclei.len();
    let mut pairs: Vec<(usize, f64, usize, f64)> = Vec::new();
    for i in (1..n).rev() {
        for j in (0..i).rev() {
            let radius_i = covalent_radius(nuclei[i].element)?;
            let radius_j = covalent_radius(nuclei[j].element)?;
            if distances[i][j] < (radius_i + radius_j) * threshold {
                pairs.push((i, radius_i, j, radius_j));
            }
        }
    }
    pairs.reverse();
    Ok(pairs
        .into_iter()
        .map(|(i, radius_i, j, radius_j)| {
            let weight = radius_i / (radius_i + radius_j);
            let from = nuclei[i].coord;
            let to = nuclei[j].coord;
            let coord = Point3d::of_tuple(
                Units::Bohr,
                (
                    from.x + (to.x - from.x) * weight,
                    from.y + (to.y - from.y) * weight,
                    from.z + (to.z - from.z) * weight,
                ),
            );
            Atom {
                element: Element::X,
                charge: Charge::of_int(0),
                coord,
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    // Coordinates divide by this constant (`15e725aa`); the module-level
    // `ANGSTROM_TO_BOHR` is the radius factor, used by `covalent_radius`.
    use qp_types::point3d::units::BOHR_TO_ANGSTROM;

    use tempfile::tempdir;

    /// Write `content` to a file in a fresh directory and return the path and
    /// the directory (kept alive by the caller).
    fn scratch(name: &str, content: &str) -> (tempfile::TempDir, std::path::PathBuf) {
        let dir = tempdir().expect("tempdir");
        let path = dir.path().join(name);
        fs::write(&path, content).expect("write fixture");
        (dir, path)
    }

    /// `read` with the defaults the reference uses (`-a` off, charge 0,
    /// singlet, no ghosts).
    fn read_default(path: &Path) -> Result<Geometry, AppError> {
        read(path, false, 0, 1, 0.0)
    }

    #[test]
    fn valid_xyz_is_read_in_angstrom() {
        let (_dir, path) = scratch("h2.xyz", "2\nh2\nH 0 0 0\nH 0 0 0.7\n");
        let geometry = read_default(&path).expect("valid xyz");
        assert_eq!(geometry.molecule.nuclei.len(), 2);
        assert_eq!(geometry.nuclei.len(), 2);
        // 0.7 angstrom is 0.7 / 0.52917721067121 bohr: coordinates are
        // converted by division since reference commit `15e725aa`, not by
        // multiplying by the reciprocal.
        let expected = 0.7 / BOHR_TO_ANGSTROM;
        assert!((geometry.molecule.nuclei[1].coord.z - expected).abs() < 1e-12);
        assert_eq!(geometry.molecule.elec_alpha, 1);
        assert_eq!(geometry.molecule.elec_beta, 1);
        assert_eq!(geometry.molecule.get_multiplicity().to_int(), 1);
    }

    #[test]
    fn au_flag_skips_the_conversion() {
        let (_dir, path) = scratch("h2.xyz", "2\nh2\nH 0 0 0\nH 0 0 0.7\n");
        let geometry = read(&path, true, 0, 1, 0.0).expect("valid xyz");
        assert!((geometry.molecule.nuclei[1].coord.z - 0.7).abs() < 1e-12);
    }

    #[test]
    fn short_xyz_is_the_uncaught_array_sub() {
        // Three announced atoms, two printed, trailing newline: `input_line`
        // yields exactly two atom lines, so `Array.sub` raises.
        let (_dir, path) = scratch("short.xyz", "3\nshort\nH 0 0 0\nH 0 0 0.7\n");
        let error = read_default(&path).expect_err("short xyz");
        assert_eq!(
            error.render(),
            "Fatal error: exception Invalid_argument(\"Array.sub\")\n"
        );
    }

    #[test]
    fn short_xyz_without_trailing_newline_is_the_same_crash() {
        let (_dir, path) = scratch("short.xyz", "3\nshort\nH 0 0 0\nH 0 0 0.7");
        let error = read_default(&path).expect_err("short xyz");
        assert_eq!(
            error.render(),
            "Fatal error: exception Invalid_argument(\"Array.sub\")\n"
        );
    }

    #[test]
    fn garbage_count_falls_back_to_the_zmatrix_reader() {
        // `zzz` is not an atom count, so the file is re-read as a Z-matrix;
        // its first line is then an element name. Measured against the
        // reference.
        let (_dir, path) = scratch("zzz.xyz", "zzz\nH 0 0 0\n");
        let error = read_default(&path).expect_err("garbage count");
        assert_eq!(
            error.render(),
            "Fatal error: exception Element.ElementError(\"Element Zzz unknown\")\n"
        );
    }

    #[test]
    fn non_positive_count_falls_back_to_the_zmatrix_reader() {
        let (_dir, path) = scratch("zero.xyz", "0\nH 0 0 0\n");
        let error = read_default(&path).expect_err("zero count");
        assert_eq!(
            error.render(),
            "Fatal error: exception Element.ElementError(\"Element 0 unknown\")\n"
        );
    }

    #[test]
    fn one_line_file_falls_back_and_reports_the_zmatrix_shape_error() {
        let (_dir, path) = scratch("empty.zmt", "");
        let error = read_default(&path).expect_err("empty file");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: cannot parse geometry from {}: expected xyz (atom count on line 1, then 'Element x y z' lines) or a Gaussian z-matrix (parse error: Syntax error)\n"
                .replace("{}", &path.display().to_string())
        );
    }

    #[test]
    fn malformed_zmatrix_line_is_a_prefixed_failure() {
        let (_dir, path) = scratch("bad4.zmt", "H 0 0 0\n");
        let error = read_default(&path).expect_err("malformed zmatrix");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: cannot parse geometry from {}: expected xyz (atom count on line 1, then 'Element x y z' lines) or a Gaussian z-matrix (parse error: Syntax error: H 0 0 0)\n"
                .replace("{}", &path.display().to_string())
        );
    }

    #[test]
    fn missing_file_matches_the_measured_bytes() {
        let dir = tempdir().expect("tempdir");
        let path = dir.path().join("missing.xyz");
        let error = read_default(&path).expect_err("missing file");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: cannot parse geometry from {}: expected xyz (atom count on line 1, then 'Element x y z' lines) or a Gaussian z-matrix ({}: No such file or directory)\n"
                .replace("{}", &path.display().to_string())
        );
    }

    #[test]
    fn zmatrix_is_converted_to_cartesian_bohr() {
        // The b2 fixture: two borons 3.0 angstrom apart, read as a Z-matrix.
        // The reference's to_xyz_string then of_xyz_string divides by
        // 0.52917721067121 exactly once (`15e725aa`).
        let (_dir, path) = scratch("b2.zmt", "b\nb 1 3.0\n");
        let geometry = read_default(&path).expect("zmatrix");
        assert_eq!(geometry.molecule.nuclei.len(), 2);
        let expected = 3.0 / BOHR_TO_ANGSTROM;
        assert!((geometry.molecule.nuclei[1].coord.z - expected).abs() < 1e-9);
    }

    #[test]
    fn zmatrix_unknown_element_escapes_unprefixed() {
        let (_dir, path) = scratch("Q.zmt", "Q\nH 1 0.7\n");
        let error = read_default(&path).expect_err("unknown element");
        assert_eq!(
            error.render(),
            "Fatal error: exception Element.ElementError(\"Element Q unknown\")\n"
        );
    }

    #[test]
    fn zero_multiplicity_fails_before_the_file_is_opened() {
        let dir = tempdir().expect("tempdir");
        let path = dir.path().join("does-not-exist.xyz");
        let error = read(&path, false, 0, 0, 0.0).expect_err("m=0");
        assert_eq!(
            error.render(),
            "Fatal error: exception Invalid_argument(\"Strictly_positive_int : (x > 0) : x=0\")\n"
        );
    }

    #[test]
    fn ghost_centre_is_r_weighted_and_between_the_nuclei() {
        // H2 at 0.7 angstrom: both radii are equal, so f = 1/2 and the ghost
        // sits at the midpoint. -g 5.0 is the h2g fixture's threshold.
        let (_dir, path) = scratch("h2.xyz", "2\nh2\nH 0 0 0\nH 0 0 0.7\n");
        let geometry = read(&path, false, 0, 1, 5.0).expect("valid xyz");
        assert_eq!(geometry.nuclei.len(), 3);
        assert_eq!(geometry.molecule.nuclei.len(), 2);
        let ghost = &geometry.nuclei[2];
        assert_eq!(ghost.element, Element::X);
        assert_eq!(ghost.charge.to_float(), 0.0);
        let expected = 0.7 / BOHR_TO_ANGSTROM / 2.0;
        assert!((ghost.coord.z - expected).abs() < 1e-12);
    }

    #[test]
    fn zero_threshold_disables_ghosts() {
        let (_dir, path) = scratch("h2.xyz", "2\nh2\nH 0 0 0\nH 0 0 0.7\n");
        let geometry = read(&path, false, 0, 1, 0.0).expect("valid xyz");
        assert_eq!(geometry.nuclei.len(), 2);
    }

    #[test]
    fn ghost_pair_order_is_i_then_j_ascending() {
        // Three hydrogens equally spaced on a line; all pairs are within the
        // threshold, so the order (1,0), (2,0), (2,1) is observable.
        let (_dir, path) = scratch("h3.xyz", "3\nh3\nH 0 0 0\nH 0 0 0.7\nH 0 0 1.4\n");
        let geometry = read(&path, false, 0, 2, 5.0).expect("valid xyz");
        let ghosts: Vec<&Atom> = geometry.nuclei[3..].iter().collect();
        assert_eq!(ghosts.len(), 3);
        let expected: Vec<f64> = vec![
            0.35 / BOHR_TO_ANGSTROM,
            0.70 / BOHR_TO_ANGSTROM,
            1.05 / BOHR_TO_ANGSTROM,
        ];
        for (ghost, z) in ghosts.iter().zip(expected) {
            assert!((ghost.coord.z - z).abs() < 1e-12);
        }
    }

    #[test]
    fn atom_error_is_unprefixed() {
        let (_dir, path) = scratch("bad.xyz", "1\nt\nH 0 0\n");
        let error = read_default(&path).expect_err("three tokens");
        assert_eq!(
            error.render(),
            "Fatal error: exception Atom.AtomError(\"H 0 0\")\n"
        );
    }

    #[test]
    fn bad_coordinate_is_a_prefixed_float_of_string() {
        // Measured: the reference also prints `Bad string: "a"` and the
        // `name = …`/`xyz = …` diagnostics before this line; the crate has
        // no such diagnostics (deviation), but the final failure is this
        // one and it names the input file.
        let (_dir, path) = scratch("bad.xyz", "2\nt\nH 0 0 0\nH a 0 0\n");
        let error = read_default(&path).expect_err("bad coordinate");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: cannot parse geometry from {}: expected xyz (atom count on line 1, then 'Element x y z' lines) or a Gaussian z-matrix (parse error: float_of_string)\n"
                .replace("{}", &path.display().to_string())
        );
    }
}
