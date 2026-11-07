//! Configuration discovery, parsing and semantic validation.
//!
//! [`discover`] implements the D3 lookup order (`--config` >
//! `$QP_EDIT_CONFIG` > `<exe_dir>/../data/qp_edit_config`) as a pure function
//! so it is testable without mutating the process environment.  [`load`]
//! parses `registry.toml` plus every module it lists and rejects anything that
//! does not satisfy the generator contract.
//!
//! Validation is intentionally strict: the files are build artifacts, so any
//! mismatch (unknown type, unknown handler, handler field marked editable,
//! formula default without its flag, ...) is a build bug that must fail
//! loudly rather than surface as surprising editor behaviour in a calculation.

use std::collections::HashSet;
use std::ffi::OsStr;
use std::fs;
use std::path::{Component, Path, PathBuf};

use anyhow::{Context, Result, bail};
use serde::de::DeserializeOwned;

use crate::handlers::Handler;
use crate::schema::{Field, ModuleConfig, ModuleRef, Registry};
use crate::types::lookup_type;
use crate::{REGISTRY_FILE, SCHEMA_VERSION};

/// A fully loaded and validated configuration directory.
#[derive(Debug, Clone)]
pub struct Config {
    root: PathBuf,
    registry: Registry,
    modules: Vec<ModuleConfig>,
}

impl Config {
    /// Directory the configuration was loaded from.
    pub fn root(&self) -> &Path {
        &self.root
    }

    /// Parsed `registry.toml`.
    pub fn registry(&self) -> &Registry {
        &self.registry
    }

    /// Modules in registry (visit) order.
    pub fn modules(&self) -> &[ModuleConfig] {
        &self.modules
    }

    /// Look a module up by id (exact match).
    pub fn module(&self, name: &str) -> Option<&ModuleConfig> {
        self.modules.iter().find(|module| module.name == name)
    }

    /// Modules that `qp_edit` visits as editable sections, in registry order.
    pub fn sections(&self) -> impl Iterator<Item = &ModuleConfig> {
        self.modules.iter().filter(|module| module.section)
    }

    /// Build a configuration from already validated parts (tests only).
    #[cfg(test)]
    pub(crate) fn from_parts(
        root: PathBuf,
        registry: Registry,
        modules: Vec<ModuleConfig>,
    ) -> Config {
        Config {
            root,
            registry,
            modules,
        }
    }
}

/// Resolve the configuration directory.
///
/// Order (D3):
///
/// 1. `explicit`, the value of `--config`;
/// 2. `env`, the value of `$QP_EDIT_CONFIG`;
/// 3. `<exe_dir>/../data/qp_edit_config`, the layout of an installed
///    `$QP_ROOT/bin/qp_edit`.
///
/// The first one present wins; it must exist and be a directory.  The
/// function is pure: the caller reads the environment and `current_exe`, so
/// the precedence is unit-testable without `set_var` (which is `unsafe` in
/// edition 2024).
pub fn discover(explicit: Option<&Path>, env: Option<&OsStr>, exe_dir: &Path) -> Result<PathBuf> {
    if let Some(path) = explicit {
        ensure_dir(path, "--config")?;
        return Ok(path.to_path_buf());
    }
    if let Some(value) = env {
        let path = PathBuf::from(value);
        ensure_dir(&path, "$QP_EDIT_CONFIG")?;
        return Ok(path);
    }
    let fallback = exe_dir.join("..").join("data").join("qp_edit_config");
    ensure_dir(&fallback, "<exe_dir>/../data/qp_edit_config")?;
    Ok(fallback)
}

/// Require `path` to be an existing directory, with an error naming how it
/// was obtained.
fn ensure_dir(path: &Path, origin: &str) -> Result<()> {
    match fs::metadata(path) {
        Ok(metadata) if metadata.is_dir() => Ok(()),
        Ok(_) => bail!("{origin} points to {path:?}, which is not a directory"),
        Err(error) => bail!("{origin} points to {path:?}: {error}"),
    }
}

/// Load and validate the configuration in `dir`.
///
/// `dir` must contain `registry.toml` and every file the registry lists.
pub fn load(dir: &Path) -> Result<Config> {
    ensure_dir(dir, "configuration directory")?;

    let registry_path = dir.join(REGISTRY_FILE);
    let registry: Registry = read_toml(&registry_path)?;
    if registry.schema_version != SCHEMA_VERSION {
        bail!(
            "{}: unsupported schema version {} (this binary understands {})",
            registry_path.display(),
            registry.schema_version,
            SCHEMA_VERSION
        );
    }

    let mut modules = Vec::with_capacity(registry.modules.len());
    let mut names = HashSet::with_capacity(registry.modules.len());
    let mut files = HashSet::with_capacity(registry.modules.len());

    for module_ref in &registry.modules {
        validate_ref(module_ref)?;
        if !names.insert(module_ref.name.clone()) {
            bail!(
                "{}: duplicate module `{}`",
                registry_path.display(),
                module_ref.name
            );
        }
        if !files.insert(module_ref.file.clone()) {
            bail!(
                "{}: duplicate file `{}`",
                registry_path.display(),
                module_ref.file
            );
        }

        let path = dir.join(&module_ref.file);
        let module: ModuleConfig = read_toml(&path).with_context(|| {
            format!(
                "module `{}` declared in {}",
                module_ref.name,
                registry_path.display()
            )
        })?;
        validate_module(&module, module_ref)
            .with_context(|| format!("module `{}` ({})", module_ref.name, path.display()))?;
        modules.push(module);
    }

    Ok(Config {
        root: dir.to_path_buf(),
        registry,
        modules,
    })
}

/// Read and parse one TOML file, keeping the path in the error chain.
fn read_toml<T: DeserializeOwned>(path: &Path) -> Result<T> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    toml::from_str(&text).with_context(|| format!("parsing {}", path.display()))
}

/// Validate one `[[module]]` entry of the registry.
fn validate_ref(module_ref: &ModuleRef) -> Result<()> {
    validate_identifier(&module_ref.name, "module name")?;
    let expected_file = format!("{}.toml", module_ref.name);
    if module_ref.file != expected_file {
        bail!(
            "module `{}`: file is {:?}, expected {:?}",
            module_ref.name,
            module_ref.file,
            expected_file
        );
    }
    if let Some(name) = &module_ref.handler {
        let handler = Handler::from_name(name)
            .with_context(|| format!("module `{}`: unknown handler", module_ref.name))?;
        if handler.name() != module_ref.name {
            bail!(
                "module `{}`: handler `{}` cannot own it (a handler name equals its module id)",
                module_ref.name,
                handler.name()
            );
        }
        if handler.is_section() != module_ref.section {
            bail!(
                "module `{}`: handler `{}` is {}a section",
                module_ref.name,
                handler.name(),
                if handler.is_section() { "" } else { "not " }
            );
        }
    }
    Ok(())
}

/// Validate one module file against its registry entry.
fn validate_module(module: &ModuleConfig, module_ref: &ModuleRef) -> Result<()> {
    if module.schema_version != SCHEMA_VERSION {
        bail!(
            "unsupported schema version {} (this binary understands {})",
            module.schema_version,
            SCHEMA_VERSION
        );
    }
    if module.name != module_ref.name {
        bail!(
            "name is {:?}, registry says {:?}",
            module.name,
            module_ref.name
        );
    }
    if module.section != module_ref.section {
        bail!(
            "section is {}, registry says {}",
            module.section,
            module_ref.section
        );
    }
    if module.handler != module_ref.handler {
        bail!(
            "handler is {:?}, registry says {:?}",
            module.handler,
            module_ref.handler
        );
    }
    if module.section && module.title.is_empty() {
        bail!("section has an empty title");
    }
    if !module.section && !module.title.is_empty() {
        bail!("non-section module has title {:?}", module.title);
    }
    validate_relative_path(&module.source_cfg, "source_cfg")?;

    if module.ezfio_dir.len() != 1 {
        bail!(
            "expected exactly one EZFIO directory (create_ezfio_stuff invariant), got {}",
            module.ezfio_dir.len()
        );
    }
    validate_identifier(&module.ezfio_dir[0], "EZFIO directory")?;

    let handler = match &module.handler {
        Some(name) => {
            let handler = Handler::from_name(name).with_context(|| "unknown handler")?;
            if handler.name() != module.name {
                bail!(
                    "handler `{}` does not own module `{}`",
                    handler.name(),
                    module.name
                );
            }
            if handler.is_section() != module.section {
                bail!("handler `{}` disagrees with `section`", handler.name());
            }
            Some(handler)
        }
        None => None,
    };

    let mut providers = HashSet::with_capacity(module.fields.len());
    for field in &module.fields {
        validate_field(field, handler).with_context(|| format!("field `{}`", field.provider))?;
        if !providers.insert(field.provider.as_str()) {
            bail!("duplicate field `{}`", field.provider);
        }
    }

    Ok(())
}

/// Validate one `[[field]]` table.
fn validate_field(field: &Field, handler: Option<Handler>) -> Result<()> {
    validate_identifier(&field.provider, "provider name")?;
    if field.ezfio_name.trim().is_empty() {
        bail!("empty ezfio_name");
    }
    if field.doc.trim().is_empty() {
        bail!("empty doc");
    }

    let spec = lookup_type(&field.type_name).ok_or_else(|| {
        anyhow::anyhow!(
            "unknown type {:?}; add it to crates/qp-input/src/types.rs",
            field.type_name
        )
    })?;

    if field.interfaces.is_empty() {
        bail!("empty interfaces list");
    }
    for interface in &field.interfaces {
        if !matches!(interface.as_str(), "ezfio" | "provider" | "ocaml") {
            bail!("unknown interface {interface:?}");
        }
    }

    if field.size.trim().is_empty() {
        bail!("empty size");
    }
    let size = field.size.trim();
    if size != "1" && !(size.starts_with('(') && size.ends_with(')')) {
        bail!("size {size:?} is neither \"1\" nor a parenthesized dimension list");
    }

    let declares_ocaml = field.interfaces.iter().any(|i| i == "ocaml");
    if field.editable && !declares_ocaml {
        bail!("field is editable but its interfaces do not contain `ocaml`");
    }
    if field.editable && handler.is_some() {
        bail!(
            "field is editable inside a handler-owned section (handler `{}` owns the whole section)",
            handler.map(Handler::name).unwrap_or("?")
        );
    }

    match (&field.default, field.default_is_formula) {
        (None, true) => bail!("default_is_formula is set but there is no default"),
        (Some(default), true) => {
            if !default.trim_start().starts_with('=') {
                bail!("default_is_formula is set but default {default:?} is not a formula");
            }
        }
        (Some(default), false) => {
            if default.trim_start().starts_with('=') {
                bail!("default {default:?} is a Fortran formula but default_is_formula is false");
            }
            spec.parse_scalar(default)
                .with_context(|| format!("default {default:?}"))?;
        }
        (None, false) => {}
    }

    Ok(())
}

/// Accept only lowercase ASCII identifiers (`[a-z0-9_]+`), which is what the
/// generator produces for module ids, EZFIO directories and providers.
fn validate_identifier(value: &str, what: &str) -> Result<()> {
    let valid = !value.is_empty()
        && value
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'_');
    if valid {
        Ok(())
    } else {
        bail!("invalid {what} {value:?} (expected [a-z0-9_]+)")
    }
}

/// Reject absolute paths and `..` components: build artifacts must be
/// relocatable and must not point outside the tree.
fn validate_relative_path(value: &str, what: &str) -> Result<()> {
    if value.is_empty() {
        bail!("empty {what}");
    }
    let path = Path::new(value);
    if path.is_absolute() {
        bail!("{what} {value:?} is absolute");
    }
    if path.components().any(|c| {
        matches!(
            c,
            Component::ParentDir | Component::RootDir | Component::Prefix(_)
        )
    }) {
        bail!("{what} {value:?} is not a plain relative path");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identifiers_reject_uppercase_and_paths() {
        assert!(validate_identifier("ao_basis", "name").is_ok());
        assert!(validate_identifier("utils2", "name").is_ok());
        assert!(validate_identifier("Ao_basis", "name").is_err());
        assert!(validate_identifier("../x", "name").is_err());
        assert!(validate_identifier("", "name").is_err());
    }

    #[test]
    fn relative_paths_reject_escapes() {
        assert!(validate_relative_path("src/ao_basis/EZFIO.cfg", "source_cfg").is_ok());
        assert!(validate_relative_path("/etc/passwd", "source_cfg").is_err());
        assert!(validate_relative_path("src/../etc/passwd", "source_cfg").is_err());
    }

    #[test]
    fn discover_prefers_explicit_then_env() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let explicit = root.join("tests/fixtures/minimal");
        let exe_dir = root.join("target");
        assert_eq!(
            discover(Some(&explicit), Some(OsStr::new("/nope")), &exe_dir).expect("explicit wins"),
            explicit
        );
        assert!(discover(Some(Path::new("/nope")), None, &exe_dir).is_err());
        assert!(discover(None, Some(OsStr::new("/nope")), &exe_dir).is_err());
        assert!(discover(None, None, Path::new("/nope")).is_err());
    }

    #[test]
    fn load_rejects_missing_registry() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let error = load(&root.join("tests/fixtures/does-not-exist")).expect_err("must fail");
        let message = format!("{error:#}");
        assert!(message.contains("does-not-exist"), "{message}");
        assert!(message.contains("configuration directory"), "{message}");
    }
}
