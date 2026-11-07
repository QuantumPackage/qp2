//! Z-matrix parsing and Cartesian reconstruction (`ocaml/Zmatrix.ml`).
//!
//! A Z-matrix is read line by line. The first atom has no reference, the
//! second gives a bond length, the third a bond length and a bond angle, and
//! every later atom a bond length, a bond angle and a dihedral angle, each
//! relative to previously defined atoms:
//!
//! ```text
//! O
//! H  1  0.96
//! H  1  0.96  2  104.5
//! ```
//!
//! Numeric fields may instead be a **label**, looked up in a map built from
//! two-token lines (`r 0.96`) placed after the atoms. A two-token line is
//! always a label definition (`Coord`), never a second atom, because the
//! reference matches the two-token case last.
//!
//! # Geometry conventions
//!
//! [`Zmatrix::to_xyz`] places the first atom at the origin and the second on
//! the `z` axis at its bond length, then builds each following atom from the
//! Euler-Rodrigues rotation matrix of [`rotation_matrix`]. The resulting
//! orientation is the reference one (a positive bond angle rotates the bond
//! vector towards negative `x`), not the "standard" one; callers must not
//! assume a particular frame.
//!
//! # Units
//!
//! Bond lengths are Bohr, angles degrees. The conversion to radians
//! (`pi/180`, with `pi = acos(-1)`) happens inside [`rotation_matrix`].
//!
//! # Preserved OCaml quirks and deviations
//!
//! * [`in_range`] wraps by **addition/subtraction only**, as the reference
//!   does. It returns the value unchanged when it is NaN, and does not
//!   terminate for an infinite input; both behaviours are inherited.
//! * Angles are wrapped to `(-180, 180]`, dihedrals to `(-360, 360]`; a value
//!   exactly equal to the lower bound is mapped to the upper bound.
//! * `distance_of_string` catches *every* exception of the reference, so a
//!   negative number is **not** an error: it becomes a label (`"-1.0"`) and
//!   is looked up in the map. Reproduced by [`distance_of_string`].
//! * [`center_of_mass`] prints `"%f %f %f\n"` to standard output before
//!   returning, exactly like the reference (which uses it as a debug trace).
//! * Two-token lines split on spaces only (`"C  1.0"` works, `"C\t1.0"` does
//!   not); a tab-separated coordinate line is a syntax error.
//! * The reference indexes its result array without a bounds check; an atom
//!   id past the end raises `Invalid_argument` and aborts. Here it is
//!   [`ZmatrixError::AtomOutOfBounds`].
//! * OCaml `int_of_string`/`float_of_string` accept hexadecimal literals,
//!   `_` separators and (for floats) hexadecimal floats; this port uses the
//!   decimal Rust parsers. Z-matrix files use plain decimals.
//! * `Zmatrix.to_string` does not exist in the reference either: lines are
//!   rendered by [`string_of_line`], which needs the label map.

use std::collections::BTreeMap;
use std::fmt;

use crate::gaussian_primitive::ocaml_fmt::format_f;
use crate::pseudo::{ElementError, ElementSymbol};

/// One-based index of an atom in the Z-matrix, as in the reference
/// `type atom_id = int`.
pub type AtomId = usize;

/// A Cartesian position in Bohr.
pub type Vector = (f64, f64, f64);

/// A 3x3 rotation matrix stored row by row; `rot[0]` is the first row.
pub type RotationMatrix = [[f64; 3]; 3];

/// An atom with its Cartesian position, in the order the reference uses:
/// `(element, x, y, z)`.
pub type Atom = (ElementSymbol, f64, f64, f64);

/// A bond length: either a value in Bohr or a label to look up.
#[derive(Debug, Clone, PartialEq)]
pub enum Distance {
    /// A literal value.
    Value(f64),
    /// A name resolved against the label map.
    Label(String),
}

/// A bond angle in degrees: either a value in `(-180, 180]` or a label.
#[derive(Debug, Clone, PartialEq)]
pub enum Angle {
    /// A literal value, already wrapped by [`angle_of_float`].
    Value(f64),
    /// A name resolved against the label map.
    Label(String),
}

/// A dihedral angle in degrees: either a value in `(-360, 360]` or a label.
#[derive(Debug, Clone, PartialEq)]
pub enum Dihedral {
    /// A literal value, already wrapped by [`dihedral_of_float`].
    Value(f64),
    /// A name resolved against the label map.
    Label(String),
}

/// One Z-matrix line.
///
/// The first four variants mirror the reference constructors `First`,
/// `Second`, `Third` and `Other`; [`Line::Coord`] is the reference `Coord`
/// (a `name value` pair that defines a label and contributes no atom).
#[derive(Debug, Clone, PartialEq)]
pub enum Line {
    /// A single atom, placed at the origin.
    First(ElementSymbol),
    /// An atom at a distance from atom 1 (the distance is not tied to an id).
    Second(ElementSymbol, Distance),
    /// An atom at a distance from `i`, with the angle `i-j-a` at `j`.
    Third(ElementSymbol, AtomId, Distance, AtomId, Angle),
    /// An atom with a distance from `i`, an angle at `j` and a dihedral
    /// around the `i-j` axis.
    Other(
        ElementSymbol,
        AtomId,
        Distance,
        AtomId,
        Angle,
        AtomId,
        Dihedral,
    ),
    /// A label definition; `Coord(name, value)`.
    Coord(String, f64),
}

/// A parsed Z-matrix: the atom lines and the label map.
///
/// This is the reference `type t = line array * float StringMap.t`. The
/// fields are public because the reference tuple is public and has no
/// invariant: an empty label map, or a line array that would fail
/// [`Zmatrix::of_string`]'s shape check, can be constructed directly.
#[derive(Debug, Clone, PartialEq)]
pub struct Zmatrix {
    /// Atom and label-definition lines, in input order.
    pub lines: Vec<Line>,
    /// Label definitions, sorted by name. Duplicate names keep the last
    /// value, like `StringMap.add`.
    pub labels: BTreeMap<String, f64>,
}

/// Failure modes of the Z-matrix parser and of [`Zmatrix::to_xyz`].
#[derive(Debug, Clone, PartialEq)]
pub enum ZmatrixError {
    /// `Element.of_string` failed while reading an atom line.
    Element(ElementError),
    /// OCaml `Failure "int_of_string"` (atom id).
    InvalidInteger {
        /// The text that failed to parse.
        text: String,
    },
    /// OCaml `assert (x > 0)` in `atom_id_of_int`, which carries no message.
    NonPositiveAtomId {
        /// The rejected value.
        value: i64,
    },
    /// OCaml `assert (x >= 0.)` in `distance_of_float`, which carries no
    /// message. Only the direct constructor reports this; the string parser
    /// turns the value into a label instead.
    NegativeDistance {
        /// The rejected value.
        value: f64,
    },
    /// OCaml `Failure "float_of_string"` (a `Coord` value, which is not
    /// protected by the label fallback).
    InvalidFloat {
        /// The text that failed to parse.
        text: String,
    },
    /// OCaml `failwith ("Syntax error: " ^ line)` from `line_of_string`.
    Syntax {
        /// The offending line, untrimmed as passed to `line_of_string`.
        text: String,
    },
    /// OCaml `failwith "Syntax error"` from `of_string`'s shape check.
    InvalidShape,
    /// Reference message `"Zmatrix error: distance %s undefined"`.
    UndefinedDistance {
        /// The label that was not found.
        label: String,
    },
    /// Reference message `"Zmatrix error: angle %s undefined"`.
    UndefinedAngle {
        /// The label that was not found.
        label: String,
    },
    /// Reference message `"Zmatrix error: dihedral %s undefined"`.
    UndefinedDihedral {
        /// The label that was not found.
        label: String,
    },
    /// Reference message `"Atom %d is defined in the future"`.
    AtomInFuture {
        /// The one-based atom id that has not been defined yet.
        atom: AtomId,
    },
    /// The reference would index out of bounds and abort (deviation).
    AtomOutOfBounds {
        /// The one-based atom id that was requested.
        atom: AtomId,
        /// Number of atom lines, the size of the reference array.
        defined: usize,
    },
    /// Reference message `"Some atoms were not defined"`.
    AtomNotDefined,
}

impl fmt::Display for ZmatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Element(error) => error.fmt(f),
            Self::InvalidInteger { .. } => f.write_str("int_of_string"),
            Self::NonPositiveAtomId { value } => {
                write!(f, "atom_id: expected a value > 0, got {value}")
            }
            Self::NegativeDistance { value } => {
                write!(f, "distance_of_float: expected a value >= 0, got {value}")
            }
            Self::InvalidFloat { .. } => f.write_str("float_of_string"),
            Self::Syntax { text } => write!(f, "Syntax error: {text}"),
            Self::InvalidShape => f.write_str("Syntax error"),
            Self::UndefinedDistance { label } => {
                write!(f, "Zmatrix error: distance {label} undefined")
            }
            Self::UndefinedAngle { label } => {
                write!(f, "Zmatrix error: angle {label} undefined")
            }
            Self::UndefinedDihedral { label } => {
                write!(f, "Zmatrix error: dihedral {label} undefined")
            }
            Self::AtomInFuture { atom } => {
                write!(f, "Atom {atom} is defined in the future")
            }
            Self::AtomOutOfBounds { atom, defined } => write!(
                f,
                "Zmatrix error: atom {atom} out of bounds ({defined} atom lines)"
            ),
            Self::AtomNotDefined => f.write_str("Some atoms were not defined"),
        }
    }
}

impl std::error::Error for ZmatrixError {}

impl From<ElementError> for ZmatrixError {
    fn from(error: ElementError) -> Self {
        Self::Element(error)
    }
}

/// Wrap `x` into `(xmin, xmax]` by adding or subtracting the period
/// `xmax - xmin`, exactly like the reference recursion.
///
/// The reference implementation is
///
/// ```text
/// if x <= xmin then in_range (x - xmin + xmax)
/// else if x > xmax then in_range (x - xmax + xmin)
/// else x
/// ```
///
/// so each step is rounded separately. A NaN input is returned unchanged
/// (both comparisons are false); an infinite input loops forever, as the
/// reference does.
#[must_use]
pub fn in_range(xmin: f64, xmax: f64, x: f64) -> f64 {
    let mut value = x;
    loop {
        if value <= xmin {
            value = value - xmin + xmax;
        } else if value > xmax {
            value = value - xmax + xmin;
        } else {
            return value;
        }
    }
}

/// `atom_id_of_int`: a one-based atom index.
///
/// # Errors
///
/// [`ZmatrixError::NonPositiveAtomId`] when `value <= 0` (the reference
/// asserts).
pub fn atom_id_of_int(value: i64) -> Result<AtomId, ZmatrixError> {
    if value <= 0 {
        return Err(ZmatrixError::NonPositiveAtomId { value });
    }
    // A positive i64 fits in a 64-bit usize; on smaller targets `try_from`
    // reports the overflow instead of wrapping.
    usize::try_from(value).map_err(|_| ZmatrixError::NonPositiveAtomId { value })
}

/// `distance_of_float`: wraps a bond length.
///
/// # Errors
///
/// [`ZmatrixError::NegativeDistance`] when `value < 0` or `value` is NaN
/// (the reference asserts `x >= 0.`).
pub fn distance_of_float(value: f64) -> Result<Distance, ZmatrixError> {
    if value >= 0.0 {
        Ok(Distance::Value(value))
    } else {
        Err(ZmatrixError::NegativeDistance { value })
    }
}

/// `angle_of_float`: wraps a bond angle into `(-180, 180]`.
#[must_use]
pub fn angle_of_float(value: f64) -> Angle {
    Angle::Value(in_range(-180.0, 180.0, value))
}

/// `dihedral_of_float`: wraps a dihedral angle into `(-360, 360]`.
#[must_use]
pub fn dihedral_of_float(value: f64) -> Dihedral {
    Dihedral::Value(in_range(-360.0, 360.0, value))
}

/// OCaml `atom_id_of_string`: parse a decimal integer (optional sign, `_`
/// separators) and require it to be positive.
///
/// # Errors
///
/// [`ZmatrixError::InvalidInteger`] when the text is not an integer, or the
/// error of [`atom_id_of_int`].
pub fn atom_id_of_string(text: &str) -> Result<AtomId, ZmatrixError> {
    parse_ocaml_int(text)
        .ok_or_else(|| ZmatrixError::InvalidInteger {
            text: text.to_owned(),
        })
        .and_then(atom_id_of_int)
}

/// OCaml `distance_of_string`: parse a float, else keep the text as a label.
///
/// The reference's `try ... with _ -> Label s` also catches the
/// non-negativity assertion, so a negative number becomes a label rather
/// than an error; that is reproduced here.
#[must_use]
pub fn distance_of_string(text: &str) -> Distance {
    match text.parse::<f64>() {
        Ok(value) if value >= 0.0 => Distance::Value(value),
        _ => Distance::Label(text.to_owned()),
    }
}

/// OCaml `angle_of_string`: parse a float, else keep the text as a label.
#[must_use]
pub fn angle_of_string(text: &str) -> Angle {
    match text.parse::<f64>() {
        Ok(value) => angle_of_float(value),
        Err(_) => Angle::Label(text.to_owned()),
    }
}

/// OCaml `dihedral_of_string`: parse a float, else keep the text as a label.
#[must_use]
pub fn dihedral_of_string(text: &str) -> Dihedral {
    match text.parse::<f64>() {
        Ok(value) => dihedral_of_float(value),
        Err(_) => Dihedral::Label(text.to_owned()),
    }
}

/// Resolve a bond length against the label map.
///
/// # Errors
///
/// [`ZmatrixError::UndefinedDistance`] for a missing label.
pub fn float_of_distance(
    labels: &BTreeMap<String, f64>,
    distance: &Distance,
) -> Result<f64, ZmatrixError> {
    match distance {
        Distance::Value(value) => Ok(*value),
        Distance::Label(label) => {
            labels
                .get(label)
                .copied()
                .ok_or_else(|| ZmatrixError::UndefinedDistance {
                    label: label.clone(),
                })
        }
    }
}

/// Resolve a bond angle against the label map.
///
/// # Errors
///
/// [`ZmatrixError::UndefinedAngle`] for a missing label.
pub fn float_of_angle(labels: &BTreeMap<String, f64>, angle: &Angle) -> Result<f64, ZmatrixError> {
    match angle {
        Angle::Value(value) => Ok(*value),
        Angle::Label(label) => {
            labels
                .get(label)
                .copied()
                .ok_or_else(|| ZmatrixError::UndefinedAngle {
                    label: label.clone(),
                })
        }
    }
}

/// Resolve a dihedral angle against the label map.
///
/// # Errors
///
/// [`ZmatrixError::UndefinedDihedral`] for a missing label.
pub fn float_of_dihedral(
    labels: &BTreeMap<String, f64>,
    dihedral: &Dihedral,
) -> Result<f64, ZmatrixError> {
    match dihedral {
        Dihedral::Value(value) => Ok(*value),
        Dihedral::Label(label) => {
            labels
                .get(label)
                .copied()
                .ok_or_else(|| ZmatrixError::UndefinedDihedral {
                    label: label.clone(),
                })
        }
    }
}

/// `int_of_atom_id`: the identity, kept for symmetry with the reference.
#[must_use]
pub fn int_of_atom_id(atom: AtomId) -> AtomId {
    atom
}

/// Render one line with the reference formats.
///
/// * `First` is `"%-3s"` (element symbol left-padded to three columns),
/// * `Second`/`Third`/`Other` use `%5d` for the ids and `%f` for the values,
/// * `Coord` is `"%s  %f"` (two spaces, no padding).
///
/// # Errors
///
/// Propagates the label-resolution errors of [`float_of_distance`],
/// [`float_of_angle`] and [`float_of_dihedral`].
pub fn string_of_line(labels: &BTreeMap<String, f64>, line: &Line) -> Result<String, ZmatrixError> {
    match line {
        Line::First(element) => Ok(format!("{:<3}", element.symbol())),
        Line::Second(element, distance) => {
            let distance = float_of_distance(labels, distance)?;
            Ok(format!(
                "{:<3} {:>5} {}",
                element.symbol(),
                1,
                format_f(distance, 0)
            ))
        }
        Line::Third(element, i, distance, j, angle) => {
            let distance = float_of_distance(labels, distance)?;
            let angle = float_of_angle(labels, angle)?;
            Ok(format!(
                "{:<3} {:>5} {} {:>5} {}",
                element.symbol(),
                int_of_atom_id(*i),
                format_f(distance, 0),
                int_of_atom_id(*j),
                format_f(angle, 0)
            ))
        }
        Line::Other(element, i, distance, j, angle, k, dihedral) => {
            let distance = float_of_distance(labels, distance)?;
            let angle = float_of_angle(labels, angle)?;
            let dihedral = float_of_dihedral(labels, dihedral)?;
            Ok(format!(
                "{:<3} {:>5} {} {:>5} {} {:>5} {}",
                element.symbol(),
                int_of_atom_id(*i),
                format_f(distance, 0),
                int_of_atom_id(*j),
                format_f(angle, 0),
                int_of_atom_id(*k),
                format_f(dihedral, 0)
            ))
        }
        Line::Coord(name, value) => Ok(format!("{name}  {}", format_f(*value, 0))),
    }
}

/// `line_of_string`: split one line on spaces and classify it.
///
/// The token count decides the variant; a two-token line is always a
/// [`Line::Coord`]. Every atom line validates its element; only `Coord` can
/// fail on its float (there is no label fallback for label values).
///
/// # Errors
///
/// [`ZmatrixError::Syntax`] for any other token count, plus the errors of
/// the element, integer and float parsers.
pub fn line_of_string(line: &str) -> Result<Line, ZmatrixError> {
    let tokens: Vec<&str> = line.split(' ').filter(|token| !token.is_empty()).collect();
    match tokens.as_slice() {
        [element] => Ok(Line::First(ElementSymbol::of_string(element)?)),
        // The middle token is an atom id in the file syntax, but the
        // reference discards it without parsing it (`e :: i :: r :: [] ->
        // Second (e, distance_of_string r)`).
        [element, _ignored_id, distance] => Ok(Line::Second(
            ElementSymbol::of_string(element)?,
            distance_of_string(distance),
        )),
        [element, i, distance, j, angle] => Ok(Line::Third(
            ElementSymbol::of_string(element)?,
            atom_id_of_string(i)?,
            distance_of_string(distance),
            atom_id_of_string(j)?,
            angle_of_string(angle),
        )),
        [element, i, distance, j, angle, k, dihedral] => Ok(Line::Other(
            ElementSymbol::of_string(element)?,
            atom_id_of_string(i)?,
            distance_of_string(distance),
            atom_id_of_string(j)?,
            angle_of_string(angle),
            atom_id_of_string(k)?,
            dihedral_of_string(dihedral),
        )),
        [name, value] => Ok(Line::Coord(
            (*name).to_owned(),
            parse_ocaml_float(value).ok_or_else(|| ZmatrixError::InvalidFloat {
                text: (*value).to_owned(),
            })?,
        )),
        _ => Err(ZmatrixError::Syntax {
            text: line.to_owned(),
        }),
    }
}

impl Zmatrix {
    /// Parse a Z-matrix.
    ///
    /// Blank lines are dropped, lines are trimmed (ASCII blanks, as OCaml
    /// `String.trim`), and the result must start with one of
    ///
    /// ```text
    /// First Second Third ...
    /// First Second Coord
    /// First Second
    /// First
    /// ```
    ///
    /// `Coord` lines are removed from [`Zmatrix::lines`] and folded into
    /// [`Zmatrix::labels`]; a later definition of the same name wins.
    ///
    /// # Errors
    ///
    /// Any error of [`line_of_string`], or
    /// [`ZmatrixError::InvalidShape`] when the shape check fails.
    pub fn of_string(text: &str) -> Result<Self, ZmatrixError> {
        let parsed: Vec<Line> = text
            .split('\n')
            .map(trim_ocaml)
            .filter(|line| !line.is_empty())
            .map(line_of_string)
            .collect::<Result<_, _>>()?;

        let shape_ok = matches!(
            parsed.as_slice(),
            [Line::First(_), Line::Second(..), Line::Third(..), ..]
                | [Line::First(_), Line::Second(..), Line::Coord(..)]
                | [Line::First(_), Line::Second(..)]
                | [Line::First(_)]
        );
        if !shape_ok {
            return Err(ZmatrixError::InvalidShape);
        }

        let mut lines = Vec::new();
        let mut labels = BTreeMap::new();
        for line in parsed {
            match line {
                Line::Coord(name, value) => {
                    labels.insert(name, value);
                }
                other => lines.push(other),
            }
        }
        Ok(Self { lines, labels })
    }

    /// Reconstruct the Cartesian coordinates.
    ///
    /// The first atom is at the origin, the second at `(0, 0, r)`, and each
    /// later atom is built from the Euler-Rodrigues rotation described in
    /// the module docs. Label-based distances, angles and dihedrals are
    /// resolved against [`Zmatrix::labels`].
    ///
    /// # Errors
    ///
    /// The resolution errors of [`string_of_line`], the reference messages
    /// [`ZmatrixError::AtomInFuture`] and [`ZmatrixError::AtomNotDefined`],
    /// and [`ZmatrixError::AtomOutOfBounds`] where the reference aborts.
    pub fn to_xyz(&self) -> Result<Vec<Atom>, ZmatrixError> {
        let mut result: Vec<Option<Atom>> = vec![None; self.lines.len()];
        for (index, line) in self.lines.iter().enumerate() {
            match line {
                Line::First(element) => {
                    result[index] = Some((*element, 0.0, 0.0, 0.0));
                }
                Line::Second(element, distance) => {
                    let distance = float_of_distance(&self.labels, distance)?;
                    result[index] = Some((*element, 0.0, 0.0, distance));
                }
                Line::Third(element, i, distance, j, angle) => {
                    let distance = float_of_distance(&self.labels, distance)?;
                    let angle = float_of_angle(&self.labels, angle)?;
                    let ui = get_cartesian_coord(&result, *i)?;
                    let uj = get_cartesian_coord(&result, *j)?;
                    let u_ij = sub(uj, ui);
                    let rotation = rotation_matrix((0.0, 1.0, 0.0), angle);
                    let new_vec =
                        apply_rotation_matrix(&rotation, scale(distance, normalized(u_ij)));
                    let (x, y, z) = add(new_vec, ui);
                    result[index] = Some((*element, x, y, z));
                }
                Line::Other(element, i, distance, j, angle, k, dihedral) => {
                    let distance = float_of_distance(&self.labels, distance)?;
                    let angle = float_of_angle(&self.labels, angle)?;
                    let dihedral = float_of_dihedral(&self.labels, dihedral)?;
                    let ui = get_cartesian_coord(&result, *i)?;
                    let uj = get_cartesian_coord(&result, *j)?;
                    let uk = get_cartesian_coord(&result, *k)?;
                    let u_ij = sub(uj, ui);
                    let u_kj = sub(uj, uk);
                    let normal = cross(u_ij, u_kj);
                    let new_vec = scale(distance, normalized(u_ij));
                    let new_vec = apply_rotation_matrix(&rotation_matrix(normal, angle), new_vec);
                    let new_vec = apply_rotation_matrix(&rotation_matrix(u_ij, dihedral), new_vec);
                    let (x, y, z) = add(new_vec, ui);
                    result[index] = Some((*element, x, y, z));
                }
                Line::Coord(..) => {}
            }
        }
        result
            .into_iter()
            .map(|atom| atom.ok_or(ZmatrixError::AtomNotDefined))
            .collect()
    }

    /// `to_xyz_string`: one `"%s %f %f %f\n"` per atom, joined with `"\n"`.
    ///
    /// The join means consecutive atoms are separated by a blank line and
    /// the string ends with a newline, exactly like the reference.
    ///
    /// # Errors
    ///
    /// Propagates [`Zmatrix::to_xyz`].
    pub fn to_xyz_string(&self) -> Result<String, ZmatrixError> {
        let atoms = self.to_xyz()?;
        let rendered: Vec<String> = atoms
            .iter()
            .map(|(element, x, y, z)| {
                format!(
                    "{} {} {} {}\n",
                    element.symbol(),
                    format_f(*x, 0),
                    format_f(*y, 0),
                    format_f(*z, 0)
                )
            })
            .collect();
        Ok(rendered.join("\n"))
    }
}

/// Look up an already reconstructed atom, reproducing the reference error
/// messages.
fn get_cartesian_coord(result: &[Option<Atom>], atom: AtomId) -> Result<Vector, ZmatrixError> {
    match result.get(atom - 1) {
        Some(Some((_, x, y, z))) => Ok((*x, *y, *z)),
        Some(None) => Err(ZmatrixError::AtomInFuture { atom }),
        None => Err(ZmatrixError::AtomOutOfBounds {
            atom,
            defined: result.len(),
        }),
    }
}

/// `u - v`, component by component.
#[must_use]
pub fn sub(u: Vector, v: Vector) -> Vector {
    (u.0 - v.0, u.1 - v.1, u.2 - v.2)
}

/// `u + v`, component by component.
#[must_use]
pub fn add(u: Vector, v: Vector) -> Vector {
    (u.0 + v.0, u.1 + v.1, u.2 + v.2)
}

/// Scalar product `s * u`, component by component.
#[must_use]
pub fn scale(s: f64, u: Vector) -> Vector {
    (s * u.0, s * u.1, s * u.2)
}

/// Dot product.
#[must_use]
pub fn dot(u: Vector, v: Vector) -> f64 {
    u.0 * v.0 + u.1 * v.1 + u.2 * v.2
}

/// Euclidean norm.
#[must_use]
pub fn norm(u: Vector) -> f64 {
    dot(u, u).sqrt()
}

/// `u / norm(u)`, or infinities/NaN for the zero vector, as the reference.
#[must_use]
pub fn normalized(u: Vector) -> Vector {
    scale(1.0 / norm(u), u)
}

/// Cross product `u x v`.
#[must_use]
pub fn cross(u: Vector, v: Vector) -> Vector {
    (
        u.1 * v.2 - u.2 * v.1,
        -(u.0 * v.2 - u.2 * v.0),
        u.0 * v.1 - u.1 * v.0,
    )
}

/// Euler-Rodrigues rotation matrix for a rotation of `angle` degrees around
/// `axis`.
///
/// Reproduces the reference transcription of
/// `github.com/jevandezande/zmatrix/converter.py`, including its sign
/// convention: the quaternion vector part is `-sin(angle/2) * axis`, so a
/// positive rotation about `+z` maps `+x` to `-y`. The axis is normalised
/// internally; a zero axis gives NaN entries rather than an error, as in the
/// reference.
#[must_use]
pub fn rotation_matrix(axis: Vector, angle: f64) -> RotationMatrix {
    let a = (angle * to_radian() * 0.5).cos();
    let (b, c, d) = scale(-(angle * to_radian() * 0.5).sin(), normalized(axis));
    [
        [
            a * a + b * b - c * c - d * d,
            2.0 * (b * c - a * d),
            2.0 * (b * d + a * c),
        ],
        [
            2.0 * (b * c + a * d),
            a * a + c * c - b * b - d * d,
            2.0 * (c * d - a * b),
        ],
        [
            2.0 * (b * d - a * c),
            2.0 * (c * d + a * b),
            a * a + d * d - b * b - c * c,
        ],
    ]
}

/// `rot * u`, i.e. the dot product of each row of `rot` with `u`.
#[must_use]
pub fn apply_rotation_matrix(rot: &RotationMatrix, u: Vector) -> Vector {
    (
        rot[0][0] * u.0 + rot[0][1] * u.1 + rot[0][2] * u.2,
        rot[1][0] * u.0 + rot[1][1] * u.1 + rot[1][2] * u.2,
        rot[2][0] * u.0 + rot[2][1] * u.1 + rot[2][2] * u.2,
    )
}

/// Centre of mass of `atoms`, in the reference units (u for the mass, Bohr
/// for the coordinates).
///
/// **Prints** `"%f %f %f\n"` for the centre to standard output, as the
/// reference does; the return value is the same triplet. A zero total mass
/// (only ghost atoms) gives infinities, as in the reference.
#[must_use]
pub fn center_of_mass(atoms: &[Atom]) -> Vector {
    let mut sum_mass = 0.0;
    let mut center = (0.0, 0.0, 0.0);
    for (element, x, y, z) in atoms {
        let mass = element.mass();
        sum_mass += mass;
        center = add(scale(mass, (*x, *y, *z)), center);
    }
    let center = scale(1.0 / sum_mass, center);
    println!(
        "{} {} {}",
        format_f(center.0, 0),
        format_f(center.1, 0),
        format_f(center.2, 0)
    );
    center
}

/// `pi / 180`, with `pi = acos(-1)` as in the reference.
fn to_radian() -> f64 {
    std::f64::consts::PI / 180.0
}

/// OCaml `String.trim`: strips spaces, form feeds, line feeds, carriage
/// returns and horizontal tabs (and nothing else).
fn trim_ocaml(text: &str) -> &str {
    text.trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\r' | '\x0c'))
}

/// A decimal subset of OCaml `int_of_string`: optional sign, decimal digits,
/// and `_` separators between digits.
///
/// OCaml also accepts `0x`/`0o`/`0b` prefixes; Z-matrix files use plain
/// decimals and the prefixes are rejected here (recorded in the port notes).
fn parse_ocaml_int(text: &str) -> Option<i64> {
    let (negative, digits) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text.strip_prefix('+').unwrap_or(text)),
    };
    if digits.is_empty() || digits.starts_with('_') || digits.ends_with('_') {
        return None;
    }
    let mut value: i64 = 0;
    for c in digits.chars() {
        if c == '_' {
            continue;
        }
        let digit = c.to_digit(10)?;
        value = value.checked_mul(10)?.checked_add(i64::from(digit))?;
    }
    Some(if negative { -value } else { value })
}

/// OCaml `float_of_string`, minus the hexadecimal-float and `_` separator
/// extensions. Rust accepts `inf`, `infinity` and `nan` case-insensitively,
/// as OCaml does.
fn parse_ocaml_float(text: &str) -> Option<f64> {
    if text.is_empty() {
        return None;
    }
    text.parse::<f64>().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1.0e-12
    }

    #[test]
    fn in_range_wraps_like_the_reference() {
        assert_eq!(in_range(-180.0, 180.0, 0.0), 0.0);
        assert_eq!(in_range(-180.0, 180.0, 180.0), 180.0);
        assert_eq!(in_range(-180.0, 180.0, -180.0), 180.0);
        assert_eq!(in_range(-180.0, 180.0, 540.0), 180.0);
        assert_eq!(in_range(-180.0, 180.0, -540.0), 180.0);
        assert!(close(in_range(-180.0, 180.0, 190.0), -170.0));
        assert!(close(in_range(-360.0, 360.0, 370.0), -350.0));
        assert!(in_range(-180.0, 180.0, f64::NAN).is_nan());
    }

    #[test]
    fn constructors_and_strings() {
        assert_eq!(distance_of_string("0.96"), Distance::Value(0.96));
        // A negative distance is caught by the reference's `with _` and
        // becomes a label.
        assert_eq!(
            distance_of_string("-1.5"),
            Distance::Label("-1.5".to_owned())
        );
        assert_eq!(distance_of_string("r1"), Distance::Label("r1".to_owned()));
        assert_eq!(angle_of_string("190.0"), Angle::Value(-170.0));
        assert_eq!(dihedral_of_string("360.0"), Dihedral::Value(360.0));
        assert_eq!(
            dihedral_of_string("1.0x"),
            Dihedral::Label("1.0x".to_owned())
        );
        assert_eq!(atom_id_of_string("3"), Ok(3));
        assert_eq!(atom_id_of_string("1_0"), Ok(10));
        assert_eq!(
            atom_id_of_string("0"),
            Err(ZmatrixError::NonPositiveAtomId { value: 0 })
        );
        assert!(matches!(
            atom_id_of_string("1.0"),
            Err(ZmatrixError::InvalidInteger { .. })
        ));
        assert!(matches!(
            distance_of_float(-0.5),
            Err(ZmatrixError::NegativeDistance { .. })
        ));
    }

    #[test]
    fn line_classification() {
        assert_eq!(
            line_of_string("O"),
            Ok(Line::First(ElementSymbol::of_string("O").unwrap()))
        );
        assert_eq!(
            line_of_string("  H   1  0.96 "),
            Ok(Line::Second(
                ElementSymbol::of_string("H").unwrap(),
                Distance::Value(0.96)
            ))
        );
        assert_eq!(
            line_of_string("H 1 0.96 2 104.5"),
            Ok(Line::Third(
                ElementSymbol::of_string("H").unwrap(),
                1,
                Distance::Value(0.96),
                2,
                Angle::Value(104.5)
            ))
        );
        // Two tokens are always a label definition, even when the first
        // token is an element symbol.
        assert_eq!(
            line_of_string("H 1.0"),
            Ok(Line::Coord("H".to_owned(), 1.0))
        );
        assert!(matches!(
            line_of_string("H 1 1.0 2"),
            Err(ZmatrixError::Syntax { .. })
        ));
        assert!(matches!(
            line_of_string("Qq"),
            Err(ZmatrixError::Element(_))
        ));
    }

    #[test]
    fn water_geometry_matches_the_analytic_rotation() {
        let zmatrix = Zmatrix::of_string("O\nH 1 0.96\nH 1 0.96 2 104.5\n").unwrap();
        assert_eq!(zmatrix.lines.len(), 3);
        assert!(zmatrix.labels.is_empty());
        let xyz = zmatrix.to_xyz().unwrap();
        assert_eq!(
            xyz[0],
            (ElementSymbol::of_string("O").unwrap(), 0.0, 0.0, 0.0)
        );
        assert_eq!(
            xyz[1],
            (ElementSymbol::of_string("H").unwrap(), 0.0, 0.0, 0.96)
        );
        let angle = 104.5_f64.to_radians();
        assert!(close(xyz[2].1, -0.96 * angle.sin()));
        assert!(close(xyz[2].2, 0.0));
        assert!(close(xyz[2].3, 0.96 * angle.cos()));
    }

    #[test]
    fn dihedral_geometry_matches_the_analytic_rotation() {
        let zmatrix =
            Zmatrix::of_string("O\nH 1 1.0\nH 1 1.0 2 90.0\nH 1 1.0 2 90.0 3 180.0\n").unwrap();
        let xyz = zmatrix.to_xyz().unwrap();
        assert!(close(xyz[2].1, -1.0));
        assert!(close(xyz[3].1, 1.0));
        assert!(close(xyz[3].2, 0.0));
        assert!(close(xyz[3].3, 0.0));
    }

    #[test]
    fn labels_and_rendering() {
        let text = "O\nH 1 r\nH 1 r 2 a\nr 0.96\na 104.5\n";
        let zmatrix = Zmatrix::of_string(text).unwrap();
        assert_eq!(zmatrix.lines.len(), 3);
        assert_eq!(zmatrix.labels.get("r"), Some(&0.96));
        assert_eq!(zmatrix.labels.get("a"), Some(&104.5));
        assert_eq!(
            zmatrix.lines[1],
            Line::Second(
                ElementSymbol::of_string("H").unwrap(),
                Distance::Label("r".to_owned())
            )
        );
        assert_eq!(
            string_of_line(&zmatrix.labels, &zmatrix.lines[1]),
            Ok("H       1 0.960000".to_owned())
        );
        assert_eq!(
            string_of_line(&zmatrix.labels, &zmatrix.lines[2]),
            Ok("H       1 0.960000     2 104.500000".to_owned())
        );
        assert_eq!(
            string_of_line(&zmatrix.labels, &Line::Coord("x".to_owned(), -1.5)),
            Ok("x  -1.500000".to_owned())
        );
        // The first atom pads its symbol to three columns.
        assert_eq!(
            string_of_line(&zmatrix.labels, &zmatrix.lines[0]),
            Ok("O  ".to_owned())
        );
        let xyz = zmatrix.to_xyz().unwrap();
        assert!(close(xyz[2].3, 0.96 * 104.5_f64.to_radians().cos()));
    }

    #[test]
    fn missing_label_reports_the_reference_message() {
        let zmatrix = Zmatrix::of_string("O\nH 1 d\n").unwrap();
        assert_eq!(
            zmatrix.to_xyz(),
            Err(ZmatrixError::UndefinedDistance {
                label: "d".to_owned()
            })
        );
        assert_eq!(
            ZmatrixError::UndefinedDistance {
                label: "d".to_owned()
            }
            .to_string(),
            "Zmatrix error: distance d undefined"
        );
    }

    #[test]
    fn shape_check_and_undefined_atoms() {
        assert_eq!(
            Zmatrix::of_string("H 1 1.0\nH 1 1.0\n"),
            Err(ZmatrixError::InvalidShape)
        );
        let future = Zmatrix::of_string("O\nH 1 1.0\nH 1 1.0 3 90.0\n").unwrap();
        assert_eq!(future.to_xyz(), Err(ZmatrixError::AtomInFuture { atom: 3 }));
        let undefined = Zmatrix {
            lines: vec![
                Line::First(ElementSymbol::of_string("O").unwrap()),
                Line::Coord("x".to_owned(), 1.0),
            ],
            labels: BTreeMap::new(),
        };
        assert_eq!(undefined.to_xyz(), Err(ZmatrixError::AtomNotDefined));
        let out_of_bounds = Zmatrix::of_string("O\nH 1 1.0\nH 1 1.0 2 90.0\n").unwrap();
        let mut lines = out_of_bounds.lines.clone();
        lines[2] = Line::Third(
            ElementSymbol::of_string("H").unwrap(),
            2,
            Distance::Value(1.0),
            9,
            Angle::Value(90.0),
        );
        let broken = Zmatrix {
            lines,
            labels: BTreeMap::new(),
        };
        assert_eq!(
            broken.to_xyz(),
            Err(ZmatrixError::AtomOutOfBounds {
                atom: 9,
                defined: 3
            })
        );
    }

    #[test]
    fn rotation_matrix_is_orthonormal_and_uses_the_reference_sign() {
        let rot = rotation_matrix((0.0, 0.0, 1.0), 90.0);
        let rotated = apply_rotation_matrix(&rot, (1.0, 0.0, 0.0));
        assert!(close(rotated.0, 0.0));
        assert!(close(rotated.1, -1.0));
        assert!(close(rotated.2, 0.0));
        // Identity at zero angle.
        let identity = rotation_matrix((1.0, 2.0, 3.0), 0.0);
        for (row, expected) in
            identity
                .iter()
                .zip([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        {
            for (got, want) in row.iter().zip(expected) {
                assert!(close(*got, want));
            }
        }
    }

    #[test]
    fn xyz_string_has_blank_line_separators() {
        let zmatrix = Zmatrix::of_string("O\nH 1 0.96\n").unwrap();
        assert_eq!(
            zmatrix.to_xyz_string(),
            Ok("O 0.000000 0.000000 0.000000\n\nH 0.000000 0.000000 0.960000\n".to_owned())
        );
    }

    #[test]
    fn center_of_mass_averages_by_mass() {
        let h = ElementSymbol::of_string("H").unwrap();
        let o = ElementSymbol::of_string("O").unwrap();
        let atoms = [(o, 0.0, 0.0, 0.0), (h, 0.0, 0.0, 1.0), (h, 0.0, 0.0, -1.0)];
        let com = center_of_mass(&atoms);
        assert!(close(com.0, 0.0));
        assert!(close(com.1, 0.0));
        assert!(close(com.2, 0.0));
    }
}
