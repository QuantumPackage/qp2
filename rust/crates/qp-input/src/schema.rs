//! Serde model of the generated `qp_edit` configuration files.
//!
//! The TOML format is emitted by `scripts/ezfio_interface/qp_edit_config.py`;
//! `rust/DESIGN-qp_edit-config.md` is the prose specification.  Every struct
//! here uses `deny_unknown_fields`: a key added by the generator without a
//! matching field in the binary must be a hard error, not silently ignored,
//! because the binary is never recompiled to follow an `EZFIO.cfg` edit.
//!
//! The `Serialize` implementations exist so the tests can round-trip a
//! fixture through `toml` and prove that the model and the format agree.

use serde::{Deserialize, Serialize};

/// Ordered list of the modules `qp_edit` visits (`registry.toml`).
#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Registry {
    /// Format version; must equal [`crate::SCHEMA_VERSION`].
    pub schema_version: u32,
    /// Modules in visit order.  Empty is syntactically valid (a build with no
    /// OCaml interface module), though the generator always writes 38.
    #[serde(default, rename = "module")]
    pub modules: Vec<ModuleRef>,
}

/// One entry of [`Registry::modules`].
#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ModuleRef {
    /// Module id, equal to the file stem and to [`ModuleConfig::name`].
    pub name: String,
    /// File name inside the configuration directory; a bare file name, never
    /// a path (the loader rejects separators).
    pub file: String,
    /// Whether `qp_edit` visits this module as an editable section.
    pub section: bool,
    /// Owning handler name, absent for generic (config-driven) modules.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub handler: Option<String>,
}

/// One module (`<module>.toml`).
#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ModuleConfig {
    /// Format version; must equal [`crate::SCHEMA_VERSION`].
    pub schema_version: u32,
    /// Module id, equal to the file stem.
    pub name: String,
    /// RST header printed by `qp_edit`; empty exactly when `section` is false.
    pub title: String,
    /// Whether `qp_edit` visits this module as an editable section.
    pub section: bool,
    /// Owning handler name, absent for generic modules.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub handler: Option<String>,
    /// `EZFIO.cfg` this module was generated from, relative to `$QP_ROOT`
    /// (never absolute, so the artifact does not depend on the build path).
    pub source_cfg: String,
    /// EZFIO directory (group) written by this module.  The generator asserts
    /// that a config has exactly one; the list keeps the schema extensible.
    pub ezfio_dir: Vec<String>,
    /// Fields in `EZFIO.cfg` order, which is the order of the generated OCaml
    /// record and therefore the order the editor displays them in.
    #[serde(default, rename = "field")]
    pub fields: Vec<Field>,
}

/// One `[[field]]` table: an EZFIO provider as declared in `EZFIO.cfg`.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Field {
    /// Provider name (the `[section]` name of `EZFIO.cfg`), used as the
    /// identifier inside a module.
    pub provider: String,
    /// Key inside the EZFIO group; defaults to the provider name but can be
    /// renamed with `ezfio_name:`.
    pub ezfio_name: String,
    /// Declared type, verbatim from `EZFIO.cfg` (`"integer"`,
    /// `"Threshold"`, `"character*(32)"`, ...).  [`crate::types::lookup_type`]
    /// maps it to storage and validation rules.
    #[serde(rename = "type")]
    pub type_name: String,
    /// Declared interfaces (`ezfio`, `provider`, `ocaml`).  Only `ocaml`
    /// matters to the editor: a field without it is not editable.
    pub interfaces: Vec<String>,
    /// Help text, verbatim from `EZFIO.cfg` (RST markup included).
    pub doc: String,
    /// Dimension expression, verbatim: `"1"` for a scalar,
    /// `"(ao_basis.ao_num,3)"` for an array.
    pub size: String,
    /// Default value as the generator renders it for OCaml (a TOML string).
    /// Absent when `EZFIO.cfg` declares no default.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub default: Option<String>,
    /// True when [`Field::default`] is a Fortran expression (starts with
    /// `=`), which the editor cannot evaluate.  Only ever set together with a
    /// `default` that starts with `=`.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub default_is_formula: bool,
    /// Whether the editor may change this field.  Always false for a module
    /// owned by a [`crate::Handler`]: the handler owns the whole section.
    pub editable: bool,
}

impl Field {
    /// Whether the field is an array (`size` other than `"1"`).
    ///
    /// The size expression is not evaluated here; the engine only needs to
    /// know whether the value is a scalar or a rank-≥1 array.
    pub fn is_array(&self) -> bool {
        self.size.trim() != "1"
    }
}
