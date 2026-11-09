//! Geometry input: the XYZ reader, its Z-matrix fallback, and ghost centres.
//!
//! `ocaml/qp_create_ezfio.ml` reads the molecule with
//! `Molecule.of_file`, which is `try of_xyz_file with XYZError ->
//! of_zmt_file`: a file with fewer than two lines, a non-numeric or
//! non-positive atom count is `XYZError` and the *same file* is re-parsed
//! as a Z-matrix (that parse then fails, e.g. with an unknown element),
//! while fewer atom lines than the announced count is an `Array.sub` in
//! the XYZ reader that the reference does not catch. This port keeps the
//! same split: the count/line-count failures fall back, the missing-lines
//! failure is fatal.
//!
//! The reference reader's `input_line` semantics are replicated with the
//! private [`input_lines`] helper (a trailing newline does not create an
//! empty last line) and the `Array.sub a 0 natoms` crash when the atom
//! lines run out early; the extracted atom lines are handed to
//! [`Atom::of_string`].
//!
//! Wave S2: every failure is one plain line,
//! `qp_create_ezfio: <message>` on standard error (exit 2); the two
//! reference classes (caught `Failure`/`Sys_error` vs uncaught
//! exceptions) had the same exit code and are collapsed.
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
//! For every pair `i > j` closer than `(r_i + r_j) * threshold`, where `r`
//! is the covalent radius in bohr, one ghost `Element::X` of charge 0 is
//! placed at `x_i + (x_j - x_i) * r_i / (r_i + r_j)` — weighted by the
//! radii, not the midpoint. Ghosts are appended after the nuclei, in `i`
//! then `j` ascending order. A zero or negative threshold disables them.

use std::path::Path;

use qp_types::error::ToolError;
use qp_types::molecule::{Molecule, MultiplicityError};
use qp_types::point3d::Units;
use qp_types::point3d::units::ANGSTROM_TO_BOHR;
use qp_types::zmatrix::{Zmatrix, ZmatrixError};
use qp_types::{Atom, Charge, Element, ElementError, Point3d, element};

use crate::TOOL;

/// The geometry as the reference's `run` holds it after ghost insertion.
#[derive(Debug, Clone)]
pub struct Geometry {
    /// The molecule as read: real nuclei only, with the electron counts
    /// derived from them. [`Molecule::name`] and the electron keys use this.
    pub molecule: Molecule,
    /// `molecule.nuclei @ ghost`: the list every later stage iterates.
    /// Real nuclei keep their order; ghost centres follow, `Element::X`
    /// and charge 0.
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
/// One [`ToolError`] per failure; the message names the input file for
/// every parse failure and carries the plain detail.
pub fn read(
    input: &Path,
    au: bool,
    charge: i64,
    multiplicity: i64,
    ghost_threshold: f64,
) -> Result<Geometry, ToolError> {
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
fn check_multiplicity(multiplicity: i64) -> Result<i64, ToolError> {
    if multiplicity >= 1 {
        Ok(multiplicity)
    } else {
        Err(ToolError::new(
            TOOL,
            format!("multiplicity {multiplicity} must be at least 1"),
        ))
    }
}

/// Read the file as text, turning an I/O failure into the geometry wrapper
/// with the standard `io::Error` display.
fn read_text(input: &Path) -> Result<String, ToolError> {
    std::fs::read_to_string(input).map_err(|error| {
        ToolError::new(
            TOOL,
            format!("cannot parse geometry {}: {error}", input.display()),
        )
    })
}

/// Why the XYZ parse did not produce a molecule.
enum XyzFailure {
    /// The reference's `XYZError`: fall back to the Z-matrix reader.
    NotXyz,
    /// Anything the reference would not have caught.
    Fatal(ToolError),
}

/// `Molecule.of_file`: XYZ first, Z-matrix on `XYZError` only.
fn parse_molecule(
    text: &str,
    input: &Path,
    charge: Charge,
    multiplicity: i64,
    units: Units,
) -> Result<Molecule, ToolError> {
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
    // Wave S2.3: decimal-only `i64::from_str` (the reference's
    // `int_of_string` also accepted hex and underscores; no golden input
    // exercises the difference).
    let Ok(natoms) = lines[0].parse::<i64>() else {
        return Err(XyzFailure::NotXyz);
    };
    if natoms <= 0 {
        return Err(XyzFailure::NotXyz);
    }
    let rest = &lines[2..];
    // `Array.sub a 0 natoms` with `natoms > Array.length a` is an
    // `Invalid_argument`, not an `XYZError`; `of_file` does not catch it.
    if i64::try_from(rest.len()).unwrap_or(i64::MAX) < natoms {
        return Err(XyzFailure::Fatal(ToolError::new(
            TOOL,
            format!(
                "too few atom lines in {}: {} lines for {natoms} announced atoms",
                input.display(),
                rest.len()
            ),
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
/// reference (whose `input_line` never returns the empty tail), so it is
/// the missing-lines failure. A plain `split('\n')` would see a third,
/// empty line and silently parse two atoms instead.
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
/// The failure messages, in the reference's evaluation order:
///
/// * a line with neither four nor five tokens is a malformed atom line;
/// * `n_elec <= 0` — the electron count the nuclei imply (minus the
///   charge) must be positive;
/// * `elec_beta < 0` — the multiplicity is too large for the electron
///   count (checked before the alpha one because OCaml evaluates the
///   tuple right-to-left);
/// * `elec_alpha < 1` — at least one alpha electron;
/// * the implied multiplicity differs — [`MultiplicityError`].
fn parse_atoms(
    text: &str,
    input: &Path,
    charge: Charge,
    multiplicity: i64,
    units: Units,
) -> Result<Molecule, ToolError> {
    let mut nuclei = Vec::new();
    for line in text.split('\n') {
        if line.is_empty() {
            continue;
        }
        // `Atom.of_string` matches a four- or five-token line and fails on
        // anything else; check the shape first so the message is the
        // reference's untrimmed line.
        let tokens: Vec<&str> = line.split(' ').filter(|token| !token.is_empty()).collect();
        if tokens.len() != 4 && tokens.len() != 5 {
            return Err(ToolError::new(
                TOOL,
                format!("malformed atom line {line:?}"),
            ));
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
        return Err(ToolError::new(
            TOOL,
            format!("the molecule has {ne} electrons (must be at least 1)"),
        ));
    }
    // `Multiplicity.to_alpha_beta` in 64-bit arithmetic: the reference's
    // `int` is 63-bit, so `-m` values past `i32` must not wrap here.
    let beta = (i64::from(ne) - multiplicity + 1) / 2;
    if beta < 0 {
        return Err(ToolError::new(
            TOOL,
            format!("multiplicity {multiplicity} is too large for {ne} electrons"),
        ));
    }
    let beta = i32::try_from(beta).expect("invariant violated: 0 <= beta <= n_elec");
    let alpha = ne - beta;
    if alpha < 1 {
        return Err(ToolError::new(
            TOOL,
            format!("the molecule has {alpha} alpha electrons (must be at least 1)"),
        ));
    }
    if 1 + alpha - beta != i32::try_from(multiplicity).unwrap_or(i32::MAX) {
        return Err(ToolError::new(
            TOOL,
            MultiplicityError {
                electrons: ne,
                multiplicity: i32::try_from(multiplicity).unwrap_or(i32::MAX),
            }
            .to_string(),
        ));
    }
    Ok(Molecule {
        nuclei,
        elec_alpha: alpha,
        elec_beta: beta,
    })
}

/// Classify an `Atom::of_string` failure into the geometry wrapper.
///
/// An unknown element keeps its plain `Element … unknown` text; a
/// coordinate or charge that does not parse says so; anything else
/// carries the innermost message of the anyhow chain.
fn atom_error(error: &anyhow::Error, input: &Path) -> ToolError {
    let detail = if let Some(element) = error.downcast_ref::<ElementError>() {
        element.to_string()
    } else if error
        .chain()
        .any(|cause| cause.downcast_ref::<std::num::ParseFloatError>().is_some())
    {
        "a coordinate or charge is not a number".to_owned()
    } else {
        error
            .chain()
            .last()
            .map_or_else(|| error.to_string(), ToString::to_string)
    };
    ToolError::new(
        TOOL,
        format!("cannot parse geometry {}: {detail}", input.display()),
    )
}

/// Wrap a Z-matrix failure with the input file; the [`ZmatrixError`]
/// `Display` texts are already plain (the `Element` arm renders
/// `Element … unknown`).
fn zmatrix_error(input: &Path, error: ZmatrixError) -> ToolError {
    ToolError::new(
        TOOL,
        format!("cannot parse geometry {}: {error}", input.display()),
    )
}

/// `Element.covalent_radius` in **bohr**.
///
/// The crate returns angstrom; the reference's radii are bohr, so this
/// multiplies by `ANGSTROM_TO_BOHR` (it does not divide). Elements from
/// berkelium on have no radius; the error keeps the crate's plain text
/// (the reference raised it from outside the geometry `try`).
fn covalent_radius(element: Element) -> Result<f64, ToolError> {
    element::covalent_radius(element)
        .map(|angstrom| angstrom * ANGSTROM_TO_BOHR)
        .map_err(|error| ToolError::new(TOOL, error.to_string()))
}

/// `ghost_centers`: bond-centred ghost atoms.
///
/// The reference walks `(i, j)` with `i` from `n-1` down to 1 and `j` from
/// `i-1` down to 0, prepending each match; the result is therefore ordered
/// by `i` then `j` ascending, which is what this function returns after
/// the explicit reverse.
///
/// # Cost
///
/// `O(n²)` time and — through [`Molecule::distance_matrix`] — `O(n²)`
/// memory, exactly like the reference. The radius lookups happen inside
/// the pair loop, so a nucleus without a covalent radius reports the
/// highest `(i, j)` pair first, as the reference does.
fn ghost_centers(molecule: &Molecule, threshold: f64) -> Result<Vec<Atom>, ToolError> {
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
    fn read_default(path: &Path) -> Result<Geometry, ToolError> {
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
    fn a_short_xyz_fails_on_the_missing_atom_lines() {
        // Three announced atoms, two printed, trailing newline: `input_line`
        // yields exactly two atom lines, so the missing-lines failure fires.
        let (_dir, path) = scratch("short.xyz", "3\nshort\nH 0 0 0\nH 0 0 0.7\n");
        let error = read_default(&path).expect_err("short xyz");
        assert_eq!(
            error.to_string(),
            format!(
                "qp_create_ezfio: too few atom lines in {}: 2 lines for 3 announced atoms",
                path.display()
            )
        );
    }

    #[test]
    fn a_short_xyz_without_trailing_newline_is_the_same_failure() {
        let (_dir, path) = scratch("short.xyz", "3\nshort\nH 0 0 0\nH 0 0 0.7");
        let error = read_default(&path).expect_err("short xyz");
        assert_eq!(
            error.to_string(),
            format!(
                "qp_create_ezfio: too few atom lines in {}: 2 lines for 3 announced atoms",
                path.display()
            )
        );
    }

    #[test]
    fn a_garbage_count_falls_back_to_the_zmatrix_reader() {
        // `zzz` is not an atom count, so the file is re-read as a Z-matrix;
        // its first line is then an element name.
        let (_dir, path) = scratch("zzz.xyz", "zzz\nH 0 0 0\n");
        let error = read_default(&path).expect_err("garbage count");
        assert_eq!(
            error.to_string(),
            format!(
                "qp_create_ezfio: cannot parse geometry {}: Element Zzz unknown",
                path.display()
            )
        );
    }

    #[test]
    fn a_non_positive_count_falls_back_to_the_zmatrix_reader() {
        let (_dir, path) = scratch("zero.xyz", "0\nH 0 0 0\n");
        let error = read_default(&path).expect_err("zero count");
        assert_eq!(
            error.to_string(),
            format!(
                "qp_create_ezfio: cannot parse geometry {}: Element 0 unknown",
                path.display()
            )
        );
    }

    #[test]
    fn a_hex_count_falls_back_to_the_zmatrix_reader() {
        // Wave S2.3: the count is decimal-only, so `0x2` is not a count and
        // the file is re-read as a Z-matrix (the reference's
        // `int_of_string` would have parsed it as 2).
        let (_dir, path) = scratch("hex.xyz", "0x2\nH 0 0 0\nH 0 0 0.7\n");
        let error = read_default(&path).expect_err("hex count");
        assert!(
            error
                .to_string()
                .starts_with("qp_create_ezfio: cannot parse geometry "),
            "{error}"
        );
    }

    #[test]
    fn a_one_line_file_falls_back_and_reports_the_zmatrix_shape_error() {
        let (_dir, path) = scratch("empty.zmt", "");
        let error = read_default(&path).expect_err("empty file");
        assert_eq!(
            error.to_string(),
            format!(
                "qp_create_ezfio: cannot parse geometry {}: Syntax error",
                path.display()
            )
        );
    }

    #[test]
    fn a_malformed_zmatrix_line_is_reported() {
        let (_dir, path) = scratch("bad4.zmt", "H 0 0 0\n");
        let error = read_default(&path).expect_err("malformed zmatrix");
        assert_eq!(
            error.to_string(),
            format!(
                "qp_create_ezfio: cannot parse geometry {}: Syntax error: H 0 0 0",
                path.display()
            )
        );
    }

    #[test]
    fn a_missing_file_is_reported_with_the_io_text() {
        let dir = tempdir().expect("tempdir");
        let path = dir.path().join("missing.xyz");
        let error = read_default(&path).expect_err("missing file");
        assert_eq!(
            error.to_string(),
            format!(
                "qp_create_ezfio: cannot parse geometry {}: No such file or directory (os error 2)",
                path.display()
            )
        );
    }

    #[test]
    fn a_zmatrix_is_converted_to_cartesian_bohr() {
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
    fn a_zmatrix_unknown_element_is_reported() {
        let (_dir, path) = scratch("Q.zmt", "Q\nH 1 0.7\n");
        let error = read_default(&path).expect_err("unknown element");
        assert_eq!(
            error.to_string(),
            format!(
                "qp_create_ezfio: cannot parse geometry {}: Element Q unknown",
                path.display()
            )
        );
    }

    #[test]
    fn a_zero_multiplicity_fails_before_the_file_is_opened() {
        let dir = tempdir().expect("tempdir");
        let path = dir.path().join("does-not-exist.xyz");
        let error = read(&path, false, 0, 0, 0.0).expect_err("m=0");
        assert_eq!(
            error.to_string(),
            "qp_create_ezfio: multiplicity 0 must be at least 1"
        );
    }

    #[test]
    fn a_ghost_centre_is_r_weighted_and_between_the_nuclei() {
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
    fn a_zero_threshold_disables_ghosts() {
        let (_dir, path) = scratch("h2.xyz", "2\nh2\nH 0 0 0\nH 0 0 0.7\n");
        let geometry = read_default(&path).expect("valid xyz");
        assert_eq!(geometry.nuclei.len(), 2);
    }

    #[test]
    fn a_ghost_pair_order_is_i_then_j_ascending() {
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
    fn an_atom_line_with_the_wrong_token_count_is_malformed() {
        let (_dir, path) = scratch("bad.xyz", "1\nt\nH 0 0\n");
        let error = read_default(&path).expect_err("three tokens");
        assert_eq!(
            error.to_string(),
            "qp_create_ezfio: malformed atom line \"H 0 0\""
        );
    }

    #[test]
    fn a_bad_coordinate_is_reported_as_not_a_number() {
        let (_dir, path) = scratch("bad.xyz", "2\nt\nH 0 0 0\nH a 0 0\n");
        let error = read_default(&path).expect_err("bad coordinate");
        assert_eq!(
            error.to_string(),
            format!(
                "qp_create_ezfio: cannot parse geometry {}: a coordinate or charge is not a number",
                path.display()
            )
        );
    }
}
