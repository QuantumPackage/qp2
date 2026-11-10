//! Config-driven input layer for the Rust `qp_edit` (QP2 tool layer).
//!
//! Design (decision D3, `rust/DESIGN-qp_edit-config.md`): `qp_edit` is
//! a **fixed binary**.  It loads a declarative description of every editable
//! EZFIO section from TOML files generated at build time from the `EZFIO.cfg`
//! inputs; editing an `EZFIO.cfg` regenerates the config without recompiling
//! the binary.  Imperative special logic (the six hand-written `Input_*`
//! behaviours) lives in [`handlers`].
//!
//! The on-disk format is described in `rust/DESIGN-qp_edit-config.md`:
//!
//! * `registry.toml` lists the modules in the order `qp_edit` visits them;
//! * `<module>.toml` describes one section: title, owning handler (if any),
//!   EZFIO directory, and one `[[field]]` table per EZFIO provider.
//!
//! The generated files are **build artifacts** (`data/qp_edit_config/`,
//! gitignored).  This crate is the only consumer.
//!
//! Units: the crate does no physics and performs no I/O to EZFIO.  Field
//! `size` strings are the dimension expressions written in `EZFIO.cfg`
//! (`"1"` for a scalar, `"(ao_basis.ao_num,3)"` for an array); they are kept
//! verbatim and are not evaluated here.
//!
//! Status: P4 complete.  [`schema`] and [`loader`] implement the generated
//! format; [`engine`] renders a section from EZFIO (`get`), parses an edited
//! RST document back into values and applies them (`set`), and keeps the
//! configuration inspection helpers used by the `qp_edit` binary.  The six
//! hand-written sections are in [`handlers`].
//!
//! Byte compatibility with the reference OCaml tool is the point of the
//! crate, so nothing here invents a format: the RST layout, the OCaml value
//! printers ([`ocaml_compat`]), the sexp parse-back ([`sexp`]) and the error
//! messages all mirror `ocaml/Input_*.ml` at the post-fix commit
//! (`15e725aa`).  [`NOTES.md`](https://github.com/QuantumPackage/qp2/blob/master/ocaml)
//! in this crate records where the two trees differ and why.

#![deny(missing_docs)]

pub mod defaults;
pub mod engine;
pub mod handlers;
pub mod loader;
pub mod ocaml_compat;
pub mod schema;
pub mod types;

/// S-expression reader used by the RST parse-back path (`qp_edit` only).
mod sexp;

pub use engine::Engine;
pub use loader::Config;
pub use schema::{ModuleConfig, Registry};

/// Version of the on-disk configuration format understood by this crate.
///
/// The generator writes this value into `registry.toml` and every
/// `<module>.toml`; [`loader`] rejects anything else so that a stale build
/// artifact fails loudly instead of being misread.
pub const SCHEMA_VERSION: u32 = 1;

/// Name of the registry file inside the configuration directory.
pub const REGISTRY_FILE: &str = "registry.toml";
