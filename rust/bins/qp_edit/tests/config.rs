//! Integration tests: the generated TOML format and the config contract.
//!
//! Every test builds its own directory under `CARGO_TARGET_TMPDIR` from the
//! checked-in fixtures and removes it on drop, so the suite is hermetic and
//! can run in parallel.  One test also drives the `qp_edit` binary end to end
//! to cover argument parsing and exit codes.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use qp_input::{Engine, ModuleConfig, Registry, loader};

/// Unique scratch directory removed on drop.
struct TempDir {
    path: PathBuf,
}

static NEXT_ID: AtomicUsize = AtomicUsize::new(0);

impl TempDir {
    fn new(tag: &str) -> TempDir {
        let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
        let path = Path::new(env!("CARGO_TARGET_TMPDIR"))
            .join(format!("qp-input-{}-{id}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).expect("create scratch directory");
        TempDir { path }
    }

    fn path(&self) -> &Path {
        &self.path
    }

    /// Copy every file of a checked-in fixture into the scratch directory.
    fn copy_fixture(&self, name: &str) {
        let source = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(name);
        for entry in fs::read_dir(&source).expect("read fixture directory") {
            let entry = entry.expect("fixture entry");
            let target = self.path.join(entry.file_name());
            fs::copy(entry.path(), target).expect("copy fixture file");
        }
    }

    fn read(&self, file: &str) -> String {
        fs::read_to_string(self.path.join(file)).expect("read scratch file")
    }

    fn write(&self, file: &str, text: &str) {
        fs::write(self.path.join(file), text).expect("write scratch file");
    }

    /// Replace the single occurrence of `from` in `file`.
    fn replace(&self, file: &str, from: &str, to: &str) {
        let text = self.read(file);
        assert_eq!(
            text.matches(from).count(),
            1,
            "expected exactly one occurrence of {from:?} in {file}"
        );
        self.write(file, &text.replacen(from, to, 1));
    }

    fn remove(&self, file: &str) {
        fs::remove_file(self.path.join(file)).expect("remove scratch file");
    }

    fn load(&self) -> anyhow::Result<qp_input::Config> {
        loader::load(self.path())
    }

    /// Load and return the error chain, asserting that loading failed.
    fn load_error(&self) -> String {
        match self.load() {
            Ok(config) => panic!(
                "configuration unexpectedly loaded with {} modules",
                config.modules().len()
            ),
            Err(error) => format!("{error:#}"),
        }
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

#[test]
fn minimal_fixture_loads_and_round_trips() {
    let dir = TempDir::new("minimal");
    dir.copy_fixture("minimal");

    let config = dir.load().expect("minimal fixture is valid");
    assert_eq!(config.modules().len(), 1);
    assert_eq!(config.sections().count(), 1);
    let module = config.module("electrons").expect("module present");
    assert_eq!(module.fields.len(), 2);
    assert!(module.fields[0].editable);
    assert!(!module.fields[1].editable);
    assert_eq!(
        module.fields[1].default.as_deref(),
        Some("= electrons.elec_alpha_num + electrons.elec_beta_num")
    );
    assert!(module.fields[1].default_is_formula);

    // The serde model and the on-disk format must agree both ways.
    let registry_text = toml::to_string(config.registry()).expect("serialize registry");
    let registry: Registry = toml::from_str(&registry_text).expect("parse registry");
    assert_eq!(&registry, config.registry());

    let module_text = toml::to_string(module).expect("serialize module");
    let reparsed: ModuleConfig = toml::from_str(&module_text).expect("parse module");
    assert_eq!(&reparsed, module);

    let engine = Engine::new(dir.load().expect("reload"));
    let text = engine.describe("electrons").expect("describe");
    assert!(text.contains("elec_alpha_num"), "{text}");
    assert!(text.contains("elec_num"), "{text}");
    assert!(text.contains("formula"), "{text}");
    assert!(text.contains("yes"), "an editable field is listed: {text}");
}

#[test]
fn warnings_name_formula_defaults() {
    let dir = TempDir::new("warnings");
    dir.copy_fixture("minimal");
    let engine = Engine::new(dir.load().expect("valid fixture"));
    let warnings = engine.warnings();
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains("elec_num"), "{warnings:?}");

    let stats = engine.stats();
    assert_eq!(stats.modules, 1);
    assert_eq!(stats.sections, 1);
    assert_eq!(stats.fields, 2);
    assert_eq!(stats.editable_fields, 1);
    assert_eq!(stats.formula_defaults, 1);
}

#[test]
fn handler_section_is_not_visited() {
    let dir = TempDir::new("handler");
    dir.copy_fixture("handler");
    let config = dir.load().expect("handler fixture is valid");
    let module = config.module("bitmasks").expect("module present");
    assert!(!module.section);
    assert_eq!(module.handler.as_deref(), Some("bitmasks"));
    let engine = Engine::new(config);
    assert!(engine.section("bitmasks").is_err());
    assert_eq!(engine.stats().sections, 0);
}

#[test]
fn unknown_key_is_rejected() {
    let dir = TempDir::new("unknown-key");
    dir.copy_fixture("minimal");
    dir.replace(
        "electrons.toml",
        "title = \"Electrons\"",
        "title = \"Electrons\"\nmystery = 1",
    );
    let error = dir.load_error();
    assert!(error.contains("unknown field"), "{error}");
    assert!(error.contains("mystery"), "{error}");
}

#[test]
fn duplicate_key_is_rejected() {
    let dir = TempDir::new("duplicate-key");
    dir.copy_fixture("minimal");
    dir.replace(
        "electrons.toml",
        "title = \"Electrons\"",
        "title = \"Electrons\"\ntitle = \"Again\"",
    );
    let error = dir.load_error();
    assert!(error.contains("duplicate key"), "{error}");
}

#[test]
fn unknown_type_is_rejected() {
    let dir = TempDir::new("bad-type");
    dir.copy_fixture("minimal");
    dir.replace(
        "electrons.toml",
        "provider = \"elec_alpha_num\"\nezfio_name = \"elec_alpha_num\"\ntype = \"Positive_int\"",
        "provider = \"elec_alpha_num\"\nezfio_name = \"elec_alpha_num\"\ntype = \"Nonsense\"",
    );
    let error = dir.load_error();
    assert!(error.contains("unknown type"), "{error}");
    assert!(error.contains("Nonsense"), "{error}");
}

#[test]
fn editable_field_must_declare_ocaml_interface() {
    let dir = TempDir::new("editable-no-ocaml");
    dir.copy_fixture("minimal");
    dir.replace(
        "electrons.toml",
        "interfaces = [\"ezfio\", \"ocaml\", \"provider\"]",
        "interfaces = [\"ezfio\", \"provider\"]",
    );
    let error = dir.load_error();
    assert!(error.contains("ocaml"), "{error}");
}

#[test]
fn editable_field_is_rejected_inside_handler_section() {
    let dir = TempDir::new("editable-handler");
    dir.copy_fixture("handler");
    dir.replace("bitmasks.toml", "editable = false", "editable = true");
    let error = dir.load_error();
    assert!(error.contains("handler-owned section"), "{error}");
    assert!(error.contains("bitmasks"), "{error}");
}

#[test]
fn formula_default_requires_the_flag() {
    let dir = TempDir::new("formula-flag");
    dir.copy_fixture("minimal");
    dir.replace("electrons.toml", "default_is_formula = true\n", "");
    let error = dir.load_error();
    assert!(error.contains("formula"), "{error}");
    assert!(error.contains("default_is_formula"), "{error}");
}

#[test]
fn formula_flag_requires_a_default() {
    let dir = TempDir::new("formula-no-default");
    dir.copy_fixture("minimal");
    dir.replace(
        "electrons.toml",
        "default = \"= electrons.elec_alpha_num + electrons.elec_beta_num\"\n",
        "",
    );
    let error = dir.load_error();
    assert!(error.contains("no default"), "{error}");
}

#[test]
fn handler_in_registry_must_match_module() {
    let dir = TempDir::new("handler-mismatch");
    dir.copy_fixture("minimal");
    dir.replace(
        "registry.toml",
        "section = true",
        "section = true\nhandler = \"electrons\"",
    );
    let error = dir.load_error();
    assert!(error.contains("registry says"), "{error}");
}

#[test]
fn unknown_handler_is_rejected() {
    let dir = TempDir::new("unknown-handler");
    dir.copy_fixture("minimal");
    dir.replace(
        "registry.toml",
        "section = true",
        "section = true\nhandler = \"nope\"",
    );
    let error = dir.load_error();
    assert!(error.contains("unknown handler"), "{error}");
}

#[test]
fn section_flag_must_match_between_registry_and_module() {
    let dir = TempDir::new("section-mismatch");
    dir.copy_fixture("minimal");
    dir.replace("registry.toml", "section = true", "section = false");
    let error = dir.load_error();
    assert!(
        error.contains("section is true, registry says false"),
        "{error}"
    );
}

#[test]
fn missing_module_file_is_rejected() {
    let dir = TempDir::new("missing-file");
    dir.copy_fixture("minimal");
    dir.remove("electrons.toml");
    let error = dir.load_error();
    assert!(error.contains("reading"), "{error}");
    assert!(error.contains("electrons.toml"), "{error}");
}

#[test]
fn schema_version_is_checked_in_registry_and_modules() {
    let registry = TempDir::new("schema-registry");
    registry.copy_fixture("minimal");
    registry.replace("registry.toml", "schema_version = 1", "schema_version = 2");
    let error = registry.load_error();
    assert!(error.contains("unsupported schema version 2"), "{error}");

    let module = TempDir::new("schema-module");
    module.copy_fixture("minimal");
    module.replace("electrons.toml", "schema_version = 1", "schema_version = 2");
    let error = module.load_error();
    assert!(error.contains("unsupported schema version 2"), "{error}");
}

#[test]
fn absolute_source_path_is_rejected() {
    let dir = TempDir::new("absolute-source");
    dir.copy_fixture("minimal");
    dir.replace(
        "electrons.toml",
        "source_cfg = \"src/electrons/EZFIO.cfg\"",
        "source_cfg = \"/etc/passwd\"",
    );
    let error = dir.load_error();
    assert!(error.contains("absolute"), "{error}");
}

#[test]
fn duplicate_module_in_registry_is_rejected() {
    let dir = TempDir::new("duplicate-module");
    dir.copy_fixture("minimal");
    dir.replace(
        "registry.toml",
        "section = true",
        "section = true\n\n[[module]]\nname = \"electrons\"\nfile = \"electrons.toml\"\nsection = true",
    );
    let error = dir.load_error();
    assert!(error.contains("duplicate module"), "{error}");
}

#[test]
fn empty_registry_is_valid() {
    let dir = TempDir::new("empty-registry");
    dir.write("registry.toml", "schema_version = 1\n");
    let config = dir
        .load()
        .expect("an empty registry is syntactically valid");
    assert!(config.modules().is_empty());
    let engine = Engine::new(config);
    assert_eq!(engine.list(), "0 modules (0 sections, 0 derived)\n");
    assert!(engine.warnings().is_empty());
}

#[test]
fn binary_lists_describes_and_checks() {
    let dir = TempDir::new("binary");
    dir.copy_fixture("minimal");
    let binary = env!("CARGO_BIN_EXE_qp_edit");

    let list = Command::new(binary)
        .arg("--config")
        .arg(dir.path())
        .arg("--list")
        .output()
        .expect("run qp_edit --list");
    assert!(list.status.success(), "{list:?}");
    let stdout = String::from_utf8(list.stdout).expect("utf-8 stdout");
    assert!(stdout.contains("electrons"), "{stdout}");

    let describe = Command::new(binary)
        .arg("--config")
        .arg(dir.path())
        .args(["--describe", "--module", "electrons"])
        .output()
        .expect("run qp_edit --describe");
    assert!(describe.status.success(), "{describe:?}");
    let stdout = String::from_utf8(describe.stdout).expect("utf-8 stdout");
    assert!(stdout.contains("elec_alpha_num"), "{stdout}");

    let check = Command::new(binary)
        .arg("--config")
        .arg(dir.path())
        .arg("--check-config")
        .output()
        .expect("run qp_edit --check-config");
    assert!(check.status.success(), "{check:?}");
    let stdout = String::from_utf8(check.stdout).expect("utf-8 stdout");
    assert!(stdout.contains("configuration OK"), "{stdout}");
    let stderr = String::from_utf8(check.stderr).expect("utf-8 stderr");
    assert!(stderr.contains("warning:"), "{stderr}");

    let missing = Command::new(binary)
        .arg("--config")
        .arg(dir.path().join("does-not-exist"))
        .arg("--list")
        .output()
        .expect("run qp_edit with a missing config");
    assert_eq!(missing.status.code(), Some(1), "{missing:?}");
}
