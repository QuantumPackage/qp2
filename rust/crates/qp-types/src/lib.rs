//! Physics-adjacent data types for the QP2 Rust tool layer.
//!
//! Direct ports of the `Qptypes`-level OCaml sources in `ocaml/`, one
//! module per source file (plus a few packagings). The contract is the
//! reference's output bytes, not its internal representation. This is a
//! **pure data** crate: no EZFIO I/O, no `unsafe`.
//!
//! # Modules
//!
//! | Rust module | OCaml source | Contents |
//! |---|---|---|
//! | [`atom`] | `Atom.ml` | a nucleus (element + coordinates) |
//! | [`basis`] | `Basis.ml`, `Long_basis.ml` | the md5-canonical basis type |
//! | [`bitlist`] | `Bitlist.ml` | fixed-width bit strings, 64 bits per word |
//! | [`charge`] | `Charge.ml` | charge states |
//! | [`command_line`] | `Command_line.ml` | CLI argument structs and usage errors |
//! | [`determinant`] | `Determinant.ml` | an alpha/beta pair of [`bitlist::Bitlist`] |
//! | [`element`] | `Element.ml` | the 119-element periodic table |
//! | [`error`] | — | the shared one-line user-facing error of the tools |
//! | [`gto`] | `Gto.ml` | Gaussian shells (re-exports the canonical [`angmom`]/[`gaussian_primitive`] types) |
//! | [`md5`] | `To_md5.ml` | the digest helper (`MD5_OF_NONE`) |
//! | [`molecule`] | `Molecule.ml`, `Multiplicity.ml` | molecule + spin multiplicity |
//! | [`point3d`] | `Point3d.ml`, `Units.ml` | coordinates + the unit boundary |
//! | [`sexplib_compat`] | sexplib0 | `Sexp.to_string` + float rules (the md5 contract) |
//! | [`angmom`] | `Angmom.ml` | shell labels `S..L`, Cartesian powers |
//! | [`gaussian_primitive`] | `GaussianPrimitive.ml` | one `(sym, expo)` Gaussian + `ocaml_fmt` |
//! | [`pseudo`] | `Pseudo.ml` | GAMESS-format pseudopotentials |
//! | [`zmatrix`] | `Zmatrix.ml` | Z-matrix parsing and Cartesian reconstruction |
//!
//! # What this crate guarantees
//!
//! 1. **The md5 contract** (`rust/CONTRACTS.md` §3): [`basis::to_md5`] is
//!    `MD5(Sexp.to_string (sexp_of_t b))` reproduced byte for byte — the
//!    serialization traps (no separator between sibling lists, one space
//!    between sibling atoms, `%.15G`/`%.17G` float rule, atom quoting)
//!    live in [`sexplib_compat`] and are pinned by the 27-row golden
//!    corpus.
//! 2. **The reference's printing and its quirks** — documented per module
//!    in `NOTES-group-a.md` / `NOTES-group-b.md`; quirks are preserved,
//!    never silently fixed.
//!
//! # Conventions
//!
//! * Bit strings are little-endian **inside** a word: bit `i` of word `w`
//!   (0-based from the least significant bit) is element `64*w + i` of the
//!   string; [`bitlist::Bitlist::to_string`] prints words in storage order
//!   (OCaml and Fortran `bit_kind_size = 64` layout).
//! * MO indices are 1-based everywhere they cross the public API.
//! * Floating-point text output reproduces the OCaml `Printf` formats used
//!   by the reference (`%f`, `%e`, `%G`); some of these strings are hashed
//!   into EZFIO md5 stamps. Helpers: [`gaussian_primitive::ocaml_fmt`].
//! * Errors are typed enums with one variant per reference failure mode,
//!   carrying the reference message text where the OCaml code exposes one.
//!
//! # Units
//!
//! Lengths are **bohr** (coordinates, [`point3d`]); energies Hartree;
//! Gaussian exponents bohr⁻²; charges in units of the elementary charge.
//! [`Units::Angstrom`] converts at the boundary only, by *division* by
//! `point3d::units::BOHR_TO_ANGSTROM` (the reference's post-fix behaviour:
//! division is a fixed point of the Angstrom round trip).
//!
#![deny(missing_docs)]

pub mod angmom;
pub mod atom;
pub mod basis;
pub mod bitlist;
pub mod charge;
pub mod command_line;
pub mod determinant;
pub mod element;
pub mod error;
pub mod gaussian_primitive;
pub mod gto;
pub mod md5;
pub mod molecule;
pub mod point3d;
pub mod pseudo;
pub mod sexplib_compat;
pub mod zmatrix;

pub use atom::Atom;
pub use basis::{Basis, LongBasis, NuclNumber};
pub use charge::Charge;
pub use element::{Element, ElementError};
pub use gto::{Angmom, AngmomXyz, Fmt, GaussianPrimitive, Gto, GtoError};
pub use md5::MD5_OF_NONE;
pub use molecule::{Molecule, Multiplicity, MultiplicityError};
pub use point3d::{Point3d, Units};
pub use sexplib_compat::Sexp;
