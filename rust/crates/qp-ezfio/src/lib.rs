//! Native Rust implementation of the EZFIO on-disk format (QP2 tool layer).
//!
//! # What EZFIO is
//!
//! An EZFIO directory is a flat set of *groups*, one subdirectory each, and a
//! group holds one file per *attribute*. There is no index: the file name is
//! the attribute name, and the payload file itself is the presence flag. A
//! `.version` file (this crate writes `2.0.8\n`) marks the root directory and
//! every group directory as EZFIO; a directory without it is not a group, even
//! if it contains files.
//!
//! * **Scalars** are plain-text files holding one value on one line:
//!   `integer`/`integer*8` in `%20d`, `double precision` in `%24.15e` (OCaml)
//!   or `E24.15` (Fortran), `logical` as `T`/`F`, `character*(n)` as the raw
//!   bytes of the string. The Fortran scalar writers are `I20`, `E24.15`,
//!   `L1` and `A`: unlike an array element, a Fortran scalar integer is *not*
//!   followed by four blanks (see [`FieldKind`]).
//! * **Arrays** are gzip-compressed files named `<key>.gz`. The decompressed
//!   payload is
//!
//!   ```text
//!   <rank: I3>\n
//!   <dim 1: I20> <dim 2: I20> ...\n
//!   <element 1>\n
//!   <element 2>\n
//!   ...
//!   ```
//!
//!   with the elements flattened in **Fortran order** (first dimension
//!   fastest). The Fortran array writers put every element in a 24-column
//!   slot: `I20` followed by 4 blanks for integers, `E24.15` for reals, `L1`
//!   followed by 23 blanks for logicals, and the raw string blank-padded to
//!   the declared `character*(n)` length. The reference OCaml client also
//!   writes one blank after the last dimension, which is recorded in
//!   [`ByteLayout`].
//!
//! # Precision
//!
//! Reals are **not** stored losslessly. The OCaml client prints 16 significant
//! digits (`%24.15e`), the Fortran client 15 (`E24.15`); a value such as
//! `1234567.8901234567` (17 significant digits) comes back as
//! `1234567.890123457` or `1234567.89012346`. Integers and logicals are exact,
//! and characters are exact up to the declared width.
//!
//! The two reference implementations (the OCaml client `ocaml/ezfio.ml` and the
//! Fortran client `external/ezfio/src/libezfio_util.irp.f`) produce
//! byte-different but mutually readable payloads. The readers here accept
//! both; [`ByteLayout`] selects which one the writers reproduce, and defaults
//! to [`ByteLayout::Ocaml`] because this crate replaces the OCaml tool layer
//! and the P0 fixture harness compares EZFIO state reference-vs-Rust.
//!
//! # Units and numeric types
//!
//! EZFIO carries no unit metadata: a value means what the Fortran module that
//! wrote it means. In QP2 that is atomic units (Hartree for energies, Bohr for
//! lengths), but the format only defines the storage classes:
//!
//! | EZFIO type | Rust type | Storage |
//! |---|---|---|
//! | `integer` (32-bit) | [`i64`] | `I20` field |
//! | `integer*8` | [`i64`] | `I20` field |
//! | `double precision` | [`f64`] | IEEE binary64, `E24.15` |
//! | `logical` | [`bool`] | `T`/`F` |
//! | `character*(n)` | [`String`] | bytes, truncated/padded to `n` by Fortran |
//!
//! Dimensions in array payloads are element counts, never byte sizes.
//!
//! # Reading, writing, and the registry
//!
//! [`EzFio`] is the directory handle: [`EzFio::open`] and [`EzFio::create`]
//! manage the root, [`EzFio::get`]/[`EzFio::set`] do scalars,
//! [`EzFio::get_array`]/[`EzFio::set_array`] do gzip arrays, and
//! [`EzFio::groups`], [`EzFio::keys`], [`EzFio::has`] and [`EzFio::remove`]
//! list and delete attributes.
//!
//! [`StaticRegistry`] holds the 332 attributes declared by the QP2 modules and
//! is used by the typed accessors to reject a call whose Rust **type**
//! contradicts the `EZFIO.cfg` declaration, and — on the write path only — a
//! write whose **shape** does. Reads never consult the declared extents: the
//! reference getters (`read_float_array`, `read_int64_array`, …) do not
//! either, so a tree whose scalars moved (e.g. `determinants.n_states` grown
//! after the arrays were written) stays readable; a section reader that needs
//! a shape carries its own check or repair. The registry is **advisory**: keys
//! written by plugins are unknown, and unknown keys are read and written
//! without checks, exactly like both reference clients. On the write path, a
//! shape check is skipped whenever a dimension expression refers to an
//! attribute that is not set (see [`DimResolver`]).
//!
//! Writes are atomic (temporary file in the destination directory, then
//! `rename`), and gzip streams are deterministic: `mtime` is 0 and no file
//! name is embedded, so writing the same array twice yields identical bytes
//! for a given `flate2` version. D6 of the plan compares decompressed payloads.
//!
//! # Example
//!
//! ```
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! use qp_ezfio::{Array, EzFio};
//!
//! let dir = tempfile::tempdir()?;
//! let ezfio = EzFio::create(dir.path().join("h2o"))?;
//!
//! // Scalars: one file each; plugin-defined keys need no registry entry.
//! ezfio.set("nuclei", "nucl_num", &2i64)?;
//! ezfio.set("nuclei", "io_nuclear_repulsion", &"O2".to_string())?;
//! // Arrays: dimensions in Fortran order, first index varies fastest.
//! ezfio.set_array("nuclei", "nucl_coord", &Array::new(vec![2, 3], vec![0.0; 6])?)?;
//!
//! assert_eq!(ezfio.get::<i64>("nuclei", "nucl_num")?, 2);
//! assert_eq!(
//!     ezfio.get_array::<f64>("nuclei", "nucl_coord")?.dims(),
//!     &[2, 3]
//! );
//! # Ok(())
//! # }
//! ```

#![deny(missing_docs)]
#![forbid(unsafe_code)]

mod array;
mod error;
mod ezfio;
mod keys;
mod layout;
mod registry;

pub use array::Array;
pub use error::{Error, Result};
pub use ezfio::{EZFIO_VERSION, EzFio, MAX_RANK};
pub use layout::{ByteLayout, Element, ElementType, FieldKind};
pub use registry::{DimExpr, DimResolver, KeySpec, Registry, StaticRegistry, resolve_dims};
