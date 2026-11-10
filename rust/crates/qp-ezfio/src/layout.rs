//! The element model the QP2 tool layer exchanges with `ezfio-core`.
//!
//! EZFIO stores every array element as one line of text. The *readers* are
//! whitespace-tolerant (the reference OCaml client does `String.trim` on every
//! line, the historical Fortran client reads the line into `character*(32)`
//! and then uses an `(I20)` / `(E24.15)` / `(L1)` / `(A)` internal read), and
//! the two reference *writers* disagree on the exact bytes:
//!
//! | element | reference OCaml client (`ocaml/ezfio.ml`) | historical Fortran client |
//! |---|---|---|
//! | `integer`, `integer*8` | `%20d` — 20 columns | `(I20)` in a 24-char slot: 20 columns + 4 blanks |
//! | `double precision` | `%24.15e` — `0.dddddddddddddddde±ee` (lower-case `e`) | `(E24.15)` — `0.dddddddddddddddE±ee` (upper-case `E`) |
//! | `logical` | scalar file: `T\n\n` (`print_bool` emits `T\n`, the write format appends a second `\n`) | scalar file: `T\n`; `(L1)` in a 24-char slot: `T` + 23 blanks |
//! | `character*(n)` | the string, unpadded | `(A)` in an `n`-char slot: blank-padded to the declared length |
//!
//! Both are mutually readable, so reading is layout-agnostic. **Writing** is
//! not done here: [`crate::EzFio`] delegates every byte to `ezfio_core` with
//! [`ezfio_core::WriterDialect::Ocaml`], because this crate replaces the OCaml
//! tool layer and the fixture harness compares EZFIO state reference-vs-Rust.
//! This module only maps a Rust value onto the typed value `ezfio-core` takes.
//!
//! Measured, not guessed: see `tests/fixtures/README.md` for the exact
//! reference binaries and commands that produced the golden payloads.
//!
//! The reference OCaml client's extra newline is part of the on-disk
//! contract: a reference-OCaml **logical scalar** file is 3 bytes
//! (`T\n\n` / `F\n\n`) — `print_bool` writes `"T\n"` and the write format
//! appends one more `\n`. Every OCaml-client-written EZFIO in the golden
//! corpus carries this form (e.g. `work/empty`, `ao_basis/ao_cartesian`), so
//! `WriterDialect::Ocaml` reproduces it.
//!
//! Note on the reference OCaml `logical` **arrays**: `print_bool` writes
//! `"T\n"` and `write_array` appends another `\n` per element, so a
//! reference-OCaml logical array contains a blank line after every element
//! — which its own reader cannot read back (`fortran_bool_of_string ""`
//! raises). No QP2 tool writes logical arrays (they come from Fortran
//! modules), so `ezfio-core` never reproduces that layout: logicals are
//! always written one field per line.

use std::fmt;

use ezfio_core::{ArrayData, Scalar, Ty};

/// The element types EZFIO can store.
///
/// This mirrors the type names of the EZFIO config format
/// (`external/ezfio/config/*.ezfio_*config`: `integer`, `integer*8`, `real`,
/// `logical`, `character`), whose storage-side counterpart is
/// [`ezfio_core::Ty`]. Fortran `integer` and
/// `integer*8` are distinct here because the EZFIO config declares them
/// separately, but they share the same on-disk representation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ElementType {
    /// Fortran `integer` (32-bit in the core, stored in `I20` fields).
    Integer,
    /// Fortran `integer*8` (64-bit, stored in `I20` fields).
    Integer64,
    /// Fortran `double precision`.
    Real,
    /// Fortran `logical`.
    Logical,
    /// Fortran `character*(n)`.
    Character,
}

impl ElementType {
    /// Human-readable name used in diagnostics.
    pub fn as_str(self) -> &'static str {
        match self {
            ElementType::Integer => "integer",
            ElementType::Integer64 => "integer*8",
            ElementType::Real => "double precision",
            ElementType::Logical => "logical",
            ElementType::Character => "character",
        }
    }

    /// The `ezfio_core` type this element is stored as.
    ///
    /// **Trap:** QP's [`ElementType::Real`] is Fortran `double precision`, so
    /// it maps to [`Ty::Double`], not to [`Ty::Real`]. EZFIO's `real` is the
    /// single-precision declaration of the config grammar; no QP2 module uses
    /// it, and `ezfio-core` stores both at double precision, so the two would
    /// differ only in the name of the ABI entry point.
    pub(crate) fn to_core_ty(self) -> Ty {
        match self {
            ElementType::Integer => Ty::Int,
            ElementType::Integer64 => Ty::I64,
            ElementType::Real => Ty::Double,
            ElementType::Logical => Ty::Logical,
            ElementType::Character => Ty::Char,
        }
    }
}

impl fmt::Display for ElementType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// A Rust value that can be stored as one EZFIO field.
///
/// Implemented for [`f64`], [`i64`], [`bool`] and [`String`]. Each
/// implementation names the [`ezfio_core`] value the element is exchanged as;
/// the on-disk bytes themselves are produced by `ezfio_core`'s writer.
pub trait Element: Sized + Clone + PartialEq + fmt::Debug {
    /// EZFIO element type this Rust type maps to.
    const TYPE: ElementType;

    /// The `ezfio_core` scalar this value is stored as.
    ///
    /// The variant must be the one [`ElementType::to_core_ty`] of
    /// [`Self::TYPE`] selects; [`Self::from_scalar`] relies on that.
    fn to_scalar(&self) -> Scalar;

    /// Extract `Self` from a scalar read with [`Self::TYPE`]'s core type.
    ///
    /// Returns `None` when `value` is not the variant [`Self::to_scalar`]
    /// produces. Callers ask `ezfio-core` for exactly that type, so `None` is
    /// an internal invariant violation, never a property of the payload.
    fn from_scalar(value: Scalar) -> Option<Self>;

    /// Pack `values` as the `ezfio_core` payload for [`Self::TYPE`].
    ///
    /// `values` are in Fortran order (first dimension fastest). This copies,
    /// because [`ArrayData`] owns its elements; the caller's array is borrowed
    /// from [`crate::Array`], which it cannot hand over.
    fn to_array_data(values: &[Self]) -> ArrayData;

    /// Unpack the elements of an array read with [`Self::TYPE`]'s core type.
    ///
    /// `None` has the same meaning as in [`Self::from_scalar`].
    fn from_array_data(data: ArrayData) -> Option<Vec<Self>>;
}

impl Element for f64 {
    const TYPE: ElementType = ElementType::Real;

    fn to_scalar(&self) -> Scalar {
        // QP's `Real` is Fortran `double precision`: `Double`, never `Real`.
        Scalar::Double(*self)
    }

    fn from_scalar(value: Scalar) -> Option<Self> {
        match value {
            Scalar::Double(v) => Some(v),
            _ => None,
        }
    }

    fn to_array_data(values: &[Self]) -> ArrayData {
        ArrayData::Double(values.to_vec())
    }

    fn from_array_data(data: ArrayData) -> Option<Vec<Self>> {
        match data {
            ArrayData::Double(v) => Some(v),
            _ => None,
        }
    }
}

impl Element for i64 {
    const TYPE: ElementType = ElementType::Integer64;

    fn to_scalar(&self) -> Scalar {
        Scalar::I64(*self)
    }

    fn from_scalar(value: Scalar) -> Option<Self> {
        match value {
            Scalar::I64(v) => Some(v),
            _ => None,
        }
    }

    fn to_array_data(values: &[Self]) -> ArrayData {
        ArrayData::I64(values.to_vec())
    }

    fn from_array_data(data: ArrayData) -> Option<Vec<Self>> {
        match data {
            ArrayData::I64(v) => Some(v),
            _ => None,
        }
    }
}

impl Element for bool {
    const TYPE: ElementType = ElementType::Logical;

    fn to_scalar(&self) -> Scalar {
        Scalar::Logical(*self)
    }

    fn from_scalar(value: Scalar) -> Option<Self> {
        match value {
            Scalar::Logical(v) => Some(v),
            _ => None,
        }
    }

    fn to_array_data(values: &[Self]) -> ArrayData {
        ArrayData::Logical(values.to_vec())
    }

    fn from_array_data(data: ArrayData) -> Option<Vec<Self>> {
        match data {
            ArrayData::Logical(v) => Some(v),
            _ => None,
        }
    }
}

impl Element for String {
    const TYPE: ElementType = ElementType::Character;

    fn to_scalar(&self) -> Scalar {
        Scalar::Str(self.clone())
    }

    fn from_scalar(value: Scalar) -> Option<Self> {
        match value {
            Scalar::Str(v) => Some(v),
            _ => None,
        }
    }

    fn to_array_data(values: &[Self]) -> ArrayData {
        ArrayData::Str(values.to_vec())
    }

    fn from_array_data(data: ArrayData) -> Option<Vec<Self>> {
        match data {
            ArrayData::Str(v) => Some(v),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The trap the delegation must not fall into: QP's `Real` is Fortran
    /// `double precision`, and EZFIO's `real` declaration is a different type.
    #[test]
    fn qp_real_maps_to_core_double() {
        assert_eq!(ElementType::Real.to_core_ty(), Ty::Double);
        assert_ne!(ElementType::Real.to_core_ty(), Ty::Real);
        assert_eq!(ElementType::Integer.to_core_ty(), Ty::Int);
        assert_eq!(ElementType::Integer64.to_core_ty(), Ty::I64);
        assert_eq!(ElementType::Logical.to_core_ty(), Ty::Logical);
        assert_eq!(ElementType::Character.to_core_ty(), Ty::Char);
    }

    #[test]
    fn scalar_conversions_round_trip() {
        assert_eq!(i64::from_scalar((-7i64).to_scalar()), Some(-7));
        assert_eq!(f64::from_scalar(2.5f64.to_scalar()), Some(2.5));
        assert_eq!(bool::from_scalar(true.to_scalar()), Some(true));
        let text = "a b".to_string();
        assert_eq!(String::from_scalar(text.to_scalar()), Some(text));
    }

    #[test]
    fn array_conversions_round_trip_in_fortran_order() {
        let ints = vec![3i64, -1, 0];
        assert_eq!(i64::from_array_data(i64::to_array_data(&ints)), Some(ints));

        let reals = vec![1.0f64, -0.0];
        assert_eq!(
            f64::from_array_data(f64::to_array_data(&reals)),
            Some(reals)
        );

        let logicals = vec![true, false];
        assert_eq!(
            bool::from_array_data(bool::to_array_data(&logicals)),
            Some(logicals)
        );

        let chars = vec!["H".to_string(), "O".to_string()];
        assert_eq!(
            String::from_array_data(String::to_array_data(&chars)),
            Some(chars)
        );
    }

    /// A value of another core type is refused, never reinterpreted.
    #[test]
    fn mismatched_core_values_are_refused() {
        assert_eq!(i64::from_scalar(Scalar::Double(1.0)), None);
        assert_eq!(f64::from_scalar(Scalar::Real(1.0)), None);
        assert_eq!(bool::from_scalar(Scalar::I64(1)), None);
        assert_eq!(String::from_scalar(Scalar::I64(1)), None);
        assert_eq!(i64::from_array_data(ArrayData::Double(vec![1.0])), None);
    }
}
