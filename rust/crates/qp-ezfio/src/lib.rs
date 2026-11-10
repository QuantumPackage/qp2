//! QP2 tool-layer API over `ezfio_core` (the ezfio-rs implementation of the
//! EZFIO on-disk format).
//!
//! This crate is what the QP2 tools link: the directory handle, the 332-key
//! registry and the directory-level operations (listing, removal). **No
//! payload byte is formatted or parsed here.** Every read and write is
//! delegated to [`ezfio_core::ops`], the same implementation the Fortran,
//! Python and OCaml layers of this build use, so the build has one EZFIO
//! format implementation. The format description below is the specification
//! that implementation follows, kept here because it is what the tools'
//! behaviour depends on.
//!
//! # What EZFIO is
//!
//! An EZFIO directory is a flat set of *groups*, one subdirectory each, and a
//! group holds one file per *attribute*. There is no index: the file name is
//! the attribute name, and the payload file itself is the presence flag. A
//! `.version` file (this crate writes `2.0.8\n`, from `ezfio_core::VERSION`)
//! marks the root directory and every group directory as EZFIO; a directory
//! without it is not a group, even if it contains files.
//!
//! * **Scalars** are plain-text files holding one value on one line:
//!   `integer`/`integer*8` in `%20d`, `double precision` in `%24.15e` (OCaml)
//!   or `E24.15` (Fortran), `logical` as `T`/`F`, `character*(n)` as the raw
//!   bytes of the string. The Fortran scalar writers are `I20`, `E24.15`,
//!   `L1` and `A`: unlike an array element, a Fortran scalar integer is *not*
//!   followed by four blanks.
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
//!   fastest). A reference Fortran array writer puts every element in a
//!   24-column slot: `I20` followed by 4 blanks for integers, `E24.15` for
//!   reals, `L1` followed by 23 blanks for logicals, and the raw string
//!   blank-padded to the declared `character*(n)` length. The reference OCaml
//!   client writes one blank after the last dimension and unpadded records.
//!
//! # Writer dialect
//!
//! The two reference clients produce mutually readable but byte-different
//! payloads, and `ezfio_core` can reproduce either ([`ezfio_core::WriterDialect`]).
//! This crate always writes the **OCaml** client's layout: it replaces that
//! tool layer, and the golden harness compares EZFIO state reference-vs-Rust
//! byte for byte. The dialect is not configurable through [`EzFio`]. Reads
//! accept both layouts, because the reader is the one implementation in
//! `ezfio_core`.
//!
//! # Precision
//!
//! Reals are **not** stored losslessly. The writer prints 16 significant
//! digits (`%24.15e`, the OCaml form); a value such as `1234567.8901234567`
//! (17 significant digits) comes back as `1234567.890123457`. A Fortran
//! payload carries 15 digits (`E24.15`), so its sixteenth is gone for good.
//! Integers and logicals are exact.
//! Characters are written verbatim but **read trimmed on both ends** (the
//! reference OCaml reader's `String.trim`), so blank padding never survives a
//! read.
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
//! | `double precision` | [`f64`] | IEEE binary64, `%24.15e` |
//! | `logical` | [`bool`] | `T`/`F` |
//! | `character*(n)` | [`String`] | bytes; written verbatim, read trimmed |
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
pub use layout::{Element, ElementType};
pub use registry::{DimExpr, DimResolver, KeySpec, Registry, StaticRegistry, resolve_dims};
