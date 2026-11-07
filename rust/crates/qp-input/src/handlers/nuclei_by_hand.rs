//! `Nuclei_by_hand`: the molecule, edited as an xyz block.
//!
//! Port of `ocaml/Input_nuclei_by_hand.ml`.  The EZFIO group is `nuclei` and
//! the record is four arrays of length `nucl_num`:
//!
//! | field | EZFIO key | element type |
//! |---|---|---|
//! | `nucl_num` | `nuclei/nucl_num` | integer scalar |
//! | `nucl_label` | `nuclei/nucl_label` | `character*(32)` array |
//! | `nucl_charge` | `nuclei/nucl_charge` | `double precision` array |
//! | `nucl_coord` | `nuclei/nucl_coord` | `double precision`, `(nucl_num,3)` |
//!
//! `nucl_coord` is stored **column major**: the reference assembles it with
//! `x.(i)`, `y.(nucl_num+i)`, `z.(2*nucl_num+i)`, i.e. all the `x`, then all
//! the `y`, then all the `z`.  It is written back the same way, flattened into
//! a rank-2 `[nucl_num; 3]` array of the same order.
//!
//! Units: coordinates are stored in **bohr** and printed/parsed in
//! **angstrom** (the xyz-file convention).  Charges are in units of the
//! elementary charge and are floats (ghost atoms, pseudo-potentials).
//!
//! # The `"\n  "` in the rendered block
//!
//! `to_rst` emits `"  <nucl_num>\n  "` followed by the atom lines (the
//! master `%16.8f` coordinate fields, see [`to_rst`]), so the text between
//! the `nucl_num` line and the first atom is a line holding two spaces.
//! `of_rst` skips exactly one line for the block title
//! (`_ :: nucl_num :: title :: lines`), so the two are inverse of each other
//! for records at the master-display grid (see [`to_rst`]).  Do not "clean
//! up" that line: the parse would then read it as an atom.

use anyhow::{Context, Result};
use qp_ezfio::EzFio;
use qp_types::atom::Atom;
use qp_types::charge::Charge;
use qp_types::element::{self, Element};
use qp_types::point3d::{Point3d, Units};

use super::qptypes;
use super::{ParseResult, apply_parsed, ezfio_array_of_list};
use crate::engine::make_header;

/// The EZFIO group this handler owns.
const GROUP: &str = "nuclei";

/// One molecule, as stored in `nuclei`.
///
/// The four vectors are parallel and all have length `nucl_num`; that is the
/// only invariant, and every constructor in this module establishes it.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NucleiByHand {
    /// Number of nuclei (`>= 1`; a `Nucl_number`).
    pub(crate) nucl_num: i64,
    /// Element of each nucleus, in file order.
    pub(crate) nucl_label: Vec<Element>,
    /// Nuclear charge of each nucleus, in units of the elementary charge.
    pub(crate) nucl_charge: Vec<Charge>,
    /// Position of each nucleus, in **bohr**.
    pub(crate) nucl_coord: Vec<Point3d>,
}

/// `read_nucl_num`: the stored count, validated against itself.
///
/// `Nucl_number.of_int ~max:(get_max ()) (get_max ())` is always in range, so
/// this can only fail on the `x >= 1` assert (an empty molecule) or on a
/// `> 10 000` warning.
fn read_nucl_num(ezfio: &EzFio) -> Result<i64> {
    let nmax = qptypes::nucl_number_get_max(ezfio)?;
    qptypes::nucl_number_of_int_max(nmax, nmax)
}

/// `write_nucl_num`.
fn write_nucl_num(ezfio: &EzFio, nucl_num: i64) -> Result<()> {
    Ok(ezfio.set::<i64>(GROUP, "nucl_num", &nucl_num)?)
}

/// `read_nucl_label`: one element symbol per nucleus.
fn read_nucl_label(ezfio: &EzFio) -> Result<Vec<Element>> {
    let array = ezfio.get_array::<String>(GROUP, "nucl_label")?;
    array
        .data()
        .iter()
        .map(|text| element::of_string(text).with_context(|| format!("Element.of_string {text:?}")))
        .collect()
}

/// `write_nucl_label ~nucl_num`: truncate or fail exactly like
/// `ezfio_array_of_list`.
fn write_nucl_label(ezfio: &EzFio, nucl_num: i64, labels: &[Element]) -> Result<()> {
    let data: Vec<String> = labels
        .iter()
        .map(|label| element::to_string(*label).to_owned())
        .collect();
    let array = ezfio_array_of_list(vec![usize_of(nucl_num)?], data)?;
    Ok(ezfio.set_array::<String>(GROUP, "nucl_label", &array)?)
}

/// `read_nucl_charge`.
fn read_nucl_charge(ezfio: &EzFio) -> Result<Vec<Charge>> {
    let array = ezfio.get_array::<f64>(GROUP, "nucl_charge")?;
    Ok(array.data().iter().map(|x| Charge::of_float(*x)).collect())
}

/// `write_nucl_charge ~nucl_num`.
fn write_nucl_charge(ezfio: &EzFio, nucl_num: i64, charges: &[Charge]) -> Result<()> {
    let data: Vec<f64> = charges.iter().map(|charge| charge.to_float()).collect();
    let array = ezfio_array_of_list(vec![usize_of(nucl_num)?], data)?;
    Ok(ezfio.set_array::<f64>(GROUP, "nucl_charge", &array)?)
}

/// `read_nucl_coord`: de-interleave the `[nucl_num; 3]` column-major array.
fn read_nucl_coord(ezfio: &EzFio) -> Result<Vec<Point3d>> {
    let nucl_num = read_nucl_num(ezfio)?;
    let count = usize_of(nucl_num)?;
    let array = ezfio.get_array::<f64>(GROUP, "nucl_coord")?;
    let raw = array.data();
    let mut result = Vec::with_capacity(count);
    for i in 0..count {
        let x = *raw
            .get(i)
            .context("Invalid_argument: nucl_coord index out of bounds (x)")?;
        let y = *raw
            .get(count + i)
            .context("Invalid_argument: nucl_coord index out of bounds (y)")?;
        let z = *raw
            .get(2 * count + i)
            .context("Invalid_argument: nucl_coord index out of bounds (z)")?;
        result.push(Point3d { x, y, z });
    }
    Ok(result)
}

/// `write_nucl_coord ~nucl_num`: `x`s, then `y`s, then `z`s.
fn write_nucl_coord(ezfio: &EzFio, nucl_num: i64, coord: &[Point3d]) -> Result<()> {
    let count = usize_of(nucl_num)?;
    let mut data = Vec::with_capacity(3 * count);
    for point in coord {
        data.push(point.x);
    }
    for point in coord {
        data.push(point.y);
    }
    for point in coord {
        data.push(point.z);
    }
    let array = ezfio_array_of_list(vec![count, 3], data)?;
    Ok(ezfio.set_array::<f64>(GROUP, "nucl_coord", &array)?)
}

/// `read`: `None` when `nuclei/nucl_num` is absent — the section then
/// disappears from the editor.
///
/// Unlike the determinants handler there is no `try ... with _ -> None` here:
/// a malformed label or a short `nucl_coord` escapes, exactly as in the
/// reference.
pub(crate) fn read(ezfio: &EzFio) -> Result<Option<NucleiByHand>> {
    if !ezfio.has(GROUP, "nucl_num")? {
        return Ok(None);
    }
    Ok(Some(NucleiByHand {
        nucl_num: read_nucl_num(ezfio)?,
        nucl_label: read_nucl_label(ezfio)?,
        nucl_charge: read_nucl_charge(ezfio)?,
        nucl_coord: read_nucl_coord(ezfio)?,
    }))
}

/// `write`: the four fields, in record order.
pub(crate) fn write(ezfio: &EzFio, record: &NucleiByHand) -> Result<()> {
    write_nucl_num(ezfio, record.nucl_num)?;
    write_nucl_label(ezfio, record.nucl_num, &record.nucl_label)?;
    write_nucl_charge(ezfio, record.nucl_num, &record.nucl_charge)?;
    write_nucl_coord(ezfio, record.nucl_num, &record.nucl_coord)?;
    Ok(())
}

/// `to_rst`: the xyz block, without the section header.
///
/// Layout: the master `Input_nuclei_by_hand.ml` `to_rst`
/// (`Printf.sprintf "  %-3s  %3d %s"`) is the target, restored by the w12
/// user decision (2026-11-06): the w10 fixed-width `%.17g` print and the
/// w9-era no-quantisation display (reference commit `15e725aa`) are both
/// retracted for the nuclei display.  The charge is `Charge.to_int`
/// (`%3d`, unsigned, no `+`); the coordinates are the master
/// `Point3d.to_string ~units:Units.Angstrom` form `%16.8f %16.8f %16.8f`
/// (the quantising display, DEVIATIONS C1 as rewritten in w12).
/// `Charge.to_string` / `Point3d::to_string` (`%.17g`) remain the
/// non-display form, e.g. in [`qp_types::atom::Atom`].
///
/// Because the display quantises, `to_rst`/`of_rst` are inverse of each
/// other for records at the master-display grid
/// ([`Point3d::to_master_grid`]); an off-grid record reads back as its grid
/// value (that is the check-mode write-back, contract C2's nuclei
/// exception, and the survival clause of contract C5 as rewritten in w12).
///
/// ```text
/// \nNuclear coordinates in xyz format (Angstroms) ::\n\n
///   <nucl_num>\n
///   \n
///   <sym padded to 3>  <charge %3d> <x y z in angstrom, %16.8f each>
///   ...
/// \n\n
/// ```
pub(crate) fn to_rst(record: &NucleiByHand) -> String {
    let mut lines: Vec<String> = Vec::with_capacity(record.nucl_coord.len());
    for i in 0..usize::try_from(record.nucl_num).unwrap_or(0) {
        let label = record.nucl_label.get(i);
        let charge = record.nucl_charge.get(i);
        let coord = record.nucl_coord.get(i);
        let (Some(label), Some(charge), Some(coord)) = (label, charge, coord) else {
            // `Array.get` would raise `Invalid_argument`; the reference's
            // `List.init` has the same bound as the record, so this is only
            // reachable for a hand-built record.  Keep the line shape.
            lines.push(String::new());
            continue;
        };
        // `Printf.sprintf "  %-3s  %3d %s"`: left-justify the symbol to 3,
        // right-justify the integer charge to 3, and the coordinates in the
        // master `%16.8f %16.8f %16.8f` Angstrom form (w12 user decision).
        lines.push(format!(
            "  {:<3}  {:>3} {}",
            element::to_string(*label),
            charge.to_int(),
            Point3d::to_string_master(Units::Angstrom, coord)
        ));
    }
    let text = format!("  {}\n  \n{}", record.nucl_num, lines.join("\n"));
    format!("\nNuclear coordinates in xyz format (Angstroms) ::\n\n{text}\n\n")
}

/// `of_rst`: parse the xyz block back into a record.
///
/// The block is located by the first line whose trimmed form is longer than
/// three characters and ends in `"::"`; the line after the next one is skipped
/// (`_ :: nucl_num :: title :: lines`, see the module docs).
///
/// # Errors
///
/// * `Err` when no `"::"` line exists (the reference's `Not_found`) or when
///   `nucl_num` does not parse as an integer (its `Failure`): both escape
///   `of_rst` and make the driver print `Info: Read error in <title>`, retry,
///   then fail;
/// * `Ok(None)`? — never: this handler reports all its failures as
///   exceptions, so it has no "already printed" path.  The return type is the
///   shared [`ParseResult`] for uniformity with the other handlers.
pub(crate) fn of_rst(ezfio: &EzFio, body: &str) -> ParseResult<NucleiByHand> {
    let lines: Vec<&str> = body.split('\n').collect();

    let mut start = None;
    for (index, line) in lines.iter().enumerate() {
        let trimmed = line.trim();
        if trimmed.len() > 3 && trimmed.ends_with("::") {
            start = Some(index + 1);
            break;
        }
    }
    let Some(start) = start else {
        // OCaml `raise Not_found`.
        anyhow::bail!("Not_found");
    };

    let nmax = qptypes::nucl_number_get_max(ezfio)?;
    let rest = &lines[start..];
    if rest.len() < 3 {
        anyhow::bail!("Error in xyz format");
    }
    // `_ :: nucl_num :: title :: lines`: the first line is the blank line that
    // follows the `::` marker, the second the atom count, the third is
    // discarded.
    let nucl_num = rest[1].trim().parse::<i64>().context("int_of_string")?;
    let nucl_num = qptypes::nucl_number_of_int_max(nucl_num, nmax)?;
    let atom_lines = &rest[3..];

    let mut atoms = Vec::with_capacity(usize_of(nucl_num)?);
    for i in 0..usize_of(nucl_num)? {
        let line = atom_lines
            .get(i)
            .context("Invalid_argument: index out of bounds in the xyz block")?;
        atoms.push(Atom::of_string(Units::Angstrom, line).context("Atom.of_string")?);
    }

    Ok(Some(NucleiByHand {
        // `List.length atom_list |> Nucl_number.of_int ~max:nmax`, i.e. the
        // count of atoms actually parsed, which is `nucl_num` here.
        nucl_num: qptypes::nucl_number_of_int_max(atoms.len() as i64, nmax)?,
        nucl_label: atoms.iter().map(|atom| atom.element).collect(),
        nucl_charge: atoms.iter().map(|atom| atom.charge).collect(),
        nucl_coord: atoms.iter().map(|atom| atom.coord).collect(),
    }))
}

/// `get` for this section: `None` renders as `""`.
pub(crate) fn get(ezfio: &EzFio, title: &str) -> Result<String> {
    match read(ezfio)? {
        Some(record) => Ok(make_header(title) + &to_rst(&record)),
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
        |record| write(ezfio, record),
    )
}

/// A negative `nucl_num` is a corrupt file, not a panic.
fn usize_of(value: i64) -> Result<usize> {
    usize::try_from(value).map_err(|_| anyhow::anyhow!("negative dimension or index: {value}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use qp_types::point3d::units;
    use std::path::Path;

    use crate::handlers::test_support;

    /// The `hcn` fixture: C at the origin, H and N on the z axis.
    ///
    /// The z values are 8-decimal angstrom readings converted back to bohr
    /// by division (`a / K`).  Each value is a fixed point of both round
    /// trips the tests exercise: the master display (`%16.8f` angstrom, the
    /// w12 grid: `to_rst`/`of_rst` are the identity) and the EZFIO array
    /// storage (`%24.15e`, 16 significant digits: `set_array`/`get_array`
    /// are the identity, so a read-back equals the written `a / K`).
    /// `2 * K` rounds to `1.05835442`, whose bohr double
    /// (`1.9999999974631935`) has 17 significant digits and loses 1 ulp in
    /// storage; `1.05835443` (one grid step over, `2.000000016360455`
    /// bohr) is the nearest 8-decimal value that survives both.
    /// `0.52917721` (`0.9999999987315967` bohr) already does.
    fn hcn_like() -> (tempfile::TempDir, EzFio) {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        ezfio.set::<i64>(GROUP, "nucl_num", &3).expect("nucl_num");
        let labels = ["C", "H", "N"];
        ezfio
            .set_array::<String>(
                GROUP,
                "nucl_label",
                &ezfio_array_of_list(vec![3], labels.map(str::to_owned).to_vec()).expect("labels"),
            )
            .expect("nucl_label");
        let charges = [6.0, 1.0, 7.0];
        ezfio
            .set_array::<f64>(
                GROUP,
                "nucl_charge",
                &ezfio_array_of_list(vec![3], charges.to_vec()).expect("charges"),
            )
            .expect("nucl_charge");
        // Column major: x1 x2 x3 y1 y2 y3 z1 z2 z3.  The 8-decimal angstrom
        // values divided by K (fixed points of the master display and of the
        // `%24.15e` storage; see the `hcn_like` docs).
        let k = units::BOHR_TO_ANGSTROM;
        let coord = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -1.05835443 / k,
            1.05835443 / k,
            0.52917721 / k,
        ];
        ezfio
            .set_array::<f64>(
                GROUP,
                "nucl_coord",
                &ezfio_array_of_list(vec![3, 3], coord.to_vec()).expect("coord"),
            )
            .expect("nucl_coord");
        (dir, ezfio)
    }

    #[test]
    fn read_assembles_the_columns_in_bohr() {
        let (_dir, ezfio) = hcn_like();
        let record = read(&ezfio).expect("read").expect("some");
        assert_eq!(record.nucl_num, 3);
        assert_eq!(record.nucl_label, vec![Element::C, Element::H, Element::N]);
        assert_eq!(record.nucl_charge[2], Charge::of_int(7));
        // The z values are the grid values: the 8-decimal angstrom readings
        // converted back by division (storage-stable, so the read-back
        // equals the written value).
        assert_eq!(
            record.nucl_coord[0],
            Point3d {
                x: 0.0,
                y: 0.0,
                z: -1.05835443 / units::BOHR_TO_ANGSTROM
            }
        );
        assert_eq!(
            record.nucl_coord[1],
            Point3d {
                x: 0.0,
                y: 0.0,
                z: 1.05835443 / units::BOHR_TO_ANGSTROM
            }
        );
        assert_eq!(
            record.nucl_coord[2],
            Point3d {
                x: 0.0,
                y: 0.0,
                z: 0.52917721 / units::BOHR_TO_ANGSTROM
            }
        );
    }

    #[test]
    fn read_returns_none_without_nucl_num() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        assert!(read(&ezfio).expect("read").is_none());
    }

    #[test]
    fn to_rst_uses_the_master_atom_line() {
        let (_dir, ezfio) = hcn_like();
        let record = read(&ezfio).expect("read").expect("some");
        // The master form (w12): `  %-3s  %3d` plus `%16.8f` Angstrom
        // coordinates.  The literal is the independent Python
        // reconstruction (`$HOME/OpenCode/tmp/w12-nuclei/rust_literals.py`):
        // K read from `point3d.rs`, `f"{v * K:16.8f}"` per coordinate,
        // `f"{charge:>3d}"` (unsigned integer charge), diffed byte for byte
        // against this literal.
        assert_eq!(
            to_rst(&record),
            "\nNuclear coordinates in xyz format (Angstroms) ::\n\n  3\n  \n  \
             C      6       0.00000000       0.00000000      -1.05835443\n  \
             H      1       0.00000000       0.00000000       1.05835443\n  \
             N      7       0.00000000       0.00000000       0.52917721\n\n"
        );
        // A ragged user-edited line (coordinate fields shorter than the
        // rendered width) still parses: the parser is unchanged (w12), and
        // the `interactive.rs` awk tests rely on it.
        let ragged = "\nNuclear coordinates in xyz format (Angstroms) ::\n\n  1\n  \n  \
                      C      6 0.01 0 0\n\n";
        let parsed = of_rst(&ezfio, ragged)
            .expect("ragged line parses")
            .expect("some");
        assert_eq!(parsed.nucl_num, 1);
        assert_eq!(parsed.nucl_label, vec![Element::C]);
        assert_eq!(parsed.nucl_charge, vec![Charge::of_int(6)]);
        assert_eq!(
            parsed.nucl_coord,
            vec![Point3d::of_tuple(Units::Angstrom, (0.01, 0.0, 0.0))]
        );
    }

    /// Gate B (w12): for every corpus fixture, the rendered block is
    /// byte-identical to the independent Python reconstruction of the
    /// master form (`$HOME/OpenCode/tmp/w12-nuclei/rust_literals.py`, with
    /// K read from `point3d.rs`, `f"{v * K:16.8f}"` per coordinate,
    /// `f"{charge:>3d}"`, the `  <n>\n  ` count prefix and the `::`
    /// header).  Block byte counts: h2 183, h3 245, c2h2 307, b2 183,
    /// h2g 245, hcn 245.
    #[test]
    fn to_rst_matches_the_independent_master_reconstruction_for_every_fixture() {
        for (name, expected) in [
            (
                "h2",
                "\nNuclear coordinates in xyz format (Angstroms) ::\n\n  2\n  \n  \
                    H      1       0.00000000       0.00000000      -0.74000000\n  \
                    H      1       0.00000000       0.00000000       0.74000000\n\n",
            ),
            (
                "h3",
                "\nNuclear coordinates in xyz format (Angstroms) ::\n\n  3\n  \n  \
                    H      1       0.00000000       0.00000000      -0.74000000\n  \
                    H      1       0.00000000       0.00000000       0.74000000\n  \
                    H      1       0.00000000       0.00000000       0.00000000\n\n",
            ),
            (
                "c2h2",
                "\nNuclear coordinates in xyz format (Angstroms) ::\n\n  4\n  \n  \
                    C      4       0.00000000       0.00000000       0.60130000\n  \
                    C      4       0.00000000       0.00000000      -0.60130000\n  \
                    H      1       0.00000000       0.00000000       1.66440000\n  \
                    H      1       0.00000000       0.00000000      -1.66440000\n\n",
            ),
            (
                "b2",
                "\nNuclear coordinates in xyz format (Angstroms) ::\n\n  2\n  \n  \
                    B      5       0.00000000       0.00000000       0.00000000\n  \
                    B      5       0.00000000       0.00000000       3.00000000\n\n",
            ),
            (
                "h2g",
                "\nNuclear coordinates in xyz format (Angstroms) ::\n\n  3\n  \n  \
                    H      1       0.00000000       0.00000000      -0.74000000\n  \
                    H      1       0.00000000       0.00000000       0.74000000\n  \
                    X      0       0.00000000       0.00000000       0.00000000\n\n",
            ),
            (
                "hcn",
                "\nNuclear coordinates in xyz format (Angstroms) ::\n\n  3\n  \n  \
                    C      6       0.00000000       0.00000000       0.00000000\n  \
                    H      1       0.00000000       0.00000000       1.06400000\n  \
                    N      7       0.00000000       0.00000000      -1.15600000\n\n",
            ),
        ] {
            let source = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../tests/golden/fixtures")
                .join(name);
            let dir = tempfile::tempdir().expect("temp dir");
            let root = dir.path().join("ezfio");
            test_support::copy_tree(&source, &root);
            let ezfio = EzFio::open(&root).expect("open the fixture");
            let record = read(&ezfio).expect("read").expect("some");
            assert_eq!(
                to_rst(&record),
                expected,
                "{name}: the master-form block differs from the independent \
                 reconstruction"
            );
        }
    }

    #[test]
    fn to_rst_round_trips_through_of_rst() {
        // The `hcn_like` record sits at the master-display grid **and**
        // survives the `%24.15e` storage (see `hcn_like`), so the round
        // trip is exact: this test pins the grid fixed point
        // (`of_rst(to_rst(record)) == record` for such a record).
        // Off-grid records read back as their grid value (contract C5 as
        // rewritten in w12; see `to_master_grid`).
        let (_dir, ezfio) = hcn_like();
        let record = read(&ezfio).expect("read").expect("some");
        let text = to_rst(&record);
        let parsed = of_rst(&ezfio, &text).expect("parse").expect("some");
        assert_eq!(parsed, record);
    }

    #[test]
    fn of_rst_takes_the_count_and_atoms_from_the_edited_text() {
        let (_dir, ezfio) = hcn_like();
        let body = "\nNuclear coordinates in xyz format (Angstroms) ::\n\n  \
                    2\n  \n  \
                    O 0. 0. 1.2\n  \
                    H 0. 0. 0.\n\n";
        let parsed = of_rst(&ezfio, body).expect("parse").expect("some");
        assert_eq!(parsed.nucl_num, 2);
        assert_eq!(parsed.nucl_label, vec![Element::O, Element::H]);
        assert_eq!(
            parsed.nucl_charge,
            vec![Charge::of_int(8), Charge::of_int(1)]
        );
    }

    #[test]
    fn of_rst_fails_without_a_marker_line() {
        let (_dir, ezfio) = hcn_like();
        let error = of_rst(&ezfio, "1\nC 0 0 0\n").expect_err("error");
        assert_eq!(error.to_string(), "Not_found");
    }

    #[test]
    fn of_rst_reports_the_xyz_message_when_the_block_is_too_short() {
        let (_dir, ezfio) = hcn_like();
        let error = of_rst(&ezfio, "\nCoordinates ::\n\n").expect_err("error");
        assert_eq!(error.to_string(), "Error in xyz format");
    }
}
