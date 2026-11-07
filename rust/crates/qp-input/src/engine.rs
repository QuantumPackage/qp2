//! Generic, config-driven engine.
//!
//! Two halves share one [`Config`]:
//!
//! * the P2 read-only introspection (`list`, `describe`, `stats`, `warnings`)
//!   used by `qp_edit --check-config`;
//! * the P4 `get`/`set` path: [`Engine::get`] renders a generic section to
//!   RST (materialising the defaults of absent keys exactly like the
//!   generated `Input_*.ml` `read` functions), and [`Engine::set`] parses the
//!   RST back with [`crate::sexp`] and writes it when `apply` is true.
//!
//! Handler-owned modules ([`crate::handlers::Handler`]) are never rendered or
//! parsed here; their behaviour cannot be expressed as a field table and
//! lives in the handler modules.
//!
//! Output is deterministic: modules follow the registry order, fields follow
//! the `EZFIO.cfg` order stored in the config, and every collection that is
//! aggregated for display is built in that order.  No hashing or unordered
//! iteration reaches the output.
//!
//! Units: values are rendered with the OCaml tool's own printers
//! (`string_of_float` is twelve significant digits, see
//! [`crate::ocaml_compat::string_of_float`]); the engine neither rounds nor
//! quantises anything.

use anyhow::{Context, Result, bail};
use qp_ezfio::EzFio;

use crate::loader::Config;
use crate::ocaml_compat::string_of_float;
use crate::schema::{Field, ModuleConfig};
use crate::sexp;
use crate::types::{Domain, TypeSpec, Value, ValueKind, lookup_type};

/// Counts for `qp_edit --check-config`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Stats {
    /// Modules in the registry.
    pub modules: usize,
    /// Modules visited as editable sections.
    pub sections: usize,
    /// Field tables across all modules (handler modules included).
    pub fields: usize,
    /// Fields the generic engine may edit (always 0 in handler sections).
    pub editable_fields: usize,
    /// Fields whose default is a Fortran formula.
    pub formula_defaults: usize,
}

/// Read-only view over a loaded configuration.
#[derive(Debug, Clone)]
pub struct Engine {
    config: Config,
}

impl Engine {
    /// Wrap a validated configuration.
    pub fn new(config: Config) -> Engine {
        Engine { config }
    }

    /// The underlying configuration.
    pub fn config(&self) -> &Config {
        &self.config
    }

    /// Counts over the whole configuration.
    pub fn stats(&self) -> Stats {
        let modules = self.config.modules().len();
        let sections = self.config.sections().count();
        let mut fields = 0;
        let mut editable_fields = 0;
        let mut formula_defaults = 0;
        for module in self.config.modules() {
            for field in &module.fields {
                fields += 1;
                if field.editable {
                    editable_fields += 1;
                }
                if field.default_is_formula {
                    formula_defaults += 1;
                }
            }
        }
        Stats {
            modules,
            sections,
            fields,
            editable_fields,
            formula_defaults,
        }
    }

    /// One line per module in registry order.
    pub fn list(&self) -> String {
        let stats = self.stats();
        let mut out = format!(
            "{} modules ({} sections, {} derived)\n",
            stats.modules,
            stats.sections,
            stats.modules - stats.sections
        );
        let width = self
            .config
            .modules()
            .iter()
            .map(|m| m.name.len())
            .max()
            .unwrap_or(0);
        for module in self.config.modules() {
            let kind = if module.section { "section" } else { "derived" };
            let handler = module.handler.as_deref().unwrap_or("-");
            out.push_str(&format!(
                "  {:<width$}  {kind:<8}  handler={handler}\n",
                module.name,
                width = width
            ));
        }
        out
    }

    /// Deterministic description of one module.
    ///
    /// Fields are printed in `EZFIO.cfg` order, which is the order the editor
    /// presents them in.  The value is a complete `String` (newline-terminated)
    /// so callers can print it or diff it in tests.
    pub fn describe(&self, name: &str) -> Result<String> {
        let module = self
            .config
            .module(name)
            .ok_or_else(|| anyhow::anyhow!("no module `{name}` in the configuration"))?;
        Ok(describe_module(module))
    }

    /// Warnings that do not invalidate the configuration.
    ///
    /// The only case today: a default that is a Fortran formula.  Those are
    /// kept verbatim in the config (three in the QP2 tree) but cannot be
    /// evaluated by the editor, so the owning handler must supply the value.
    pub fn warnings(&self) -> Vec<String> {
        let mut warnings = Vec::new();
        for module in self.config.modules() {
            for field in &module.fields {
                if field.default_is_formula {
                    warnings.push(format!(
                        "module `{}`: field `{}` has a Fortran-formula default that the editor \
                         cannot evaluate; the handler must supply the value",
                        module.name, field.provider
                    ));
                }
            }
        }
        warnings
    }

    /// Reject a module id that is not a section, for the interactive path
    /// (P4).  Kept here so the rule lives next to the data it needs.
    pub fn section(&self, name: &str) -> Result<&ModuleConfig> {
        match self.config.module(name) {
            Some(module) if module.section => Ok(module),
            Some(_) => bail!("module `{name}` is not visited by qp_edit"),
            None => bail!("no module `{name}` in the configuration"),
        }
    }

    /// Whether check mode (`-c`) applies this section's parsed record to the
    /// EZFIO: the per-section flag of the owning handler,
    /// [`crate::handlers::Handler::applies_in_check`].
    ///
    /// The default is `false` (contract C2: nothing is written in check
    /// mode); since w12 the only section that answers `true` is the nuclei
    /// one.  An unknown or non-section module, and a section without a
    /// handler (generic modules have no such flag), answer `false`: `set`
    /// would reject or skip it anyway.  Resolves the handler the same way
    /// [`Engine::set`] does; the loader rejects unknown handler names, so a
    /// failure is an internal inconsistency that the `false` default
    /// absorbs.
    pub fn section_applies_in_check(&self, name: &str) -> bool {
        let Ok(module) = self.section(name) else {
            return false;
        };
        let Some(handler_name) = &module.handler else {
            return false;
        };
        handler_of(module, handler_name).is_ok_and(|handler| handler.applies_in_check())
    }

    /// The sections a run visits, given the reference's canonical task list.
    ///
    /// The canonical list fixes the document order (`qp_edit.ml`'s `tasks`),
    /// which is load-bearing: the sections are located by their header, and
    /// `determinants_by_hand`'s `of_rst` consumes the rest of the document.
    /// Two adjustments make the list follow the installed tree:
    ///
    /// * a canonical name absent from this configuration is skipped — a
    ///   build without the `qmcchem` plugin has no `Qmcchem` section, and
    ///   rendering it would abort the editor;
    /// * a section of the configuration that the canonical list does not
    ///   know (a plugin added after this binary was built) is inserted in
    ///   registry order **before the last canonical section**, so that
    ///   `determinants_by_hand` stays last.  The registry order is the
    ///   generator's sorted order, and the position of an unknown section
    ///   inside the document is otherwise unconstrained.
    ///
    /// The result borrows from `self` and from `canonical`; it is
    /// deterministic, and empty when the configuration has no section.
    pub fn visit_order<'a>(&'a self, canonical: &'a [&'a str]) -> Vec<&'a str> {
        let mut order: Vec<&str> = canonical
            .iter()
            .copied()
            .filter(|name| {
                self.config
                    .module(name)
                    .is_some_and(|module| module.section)
            })
            .collect();
        let extras: Vec<&str> = self
            .config
            .sections()
            .map(|module| module.name.as_str())
            .filter(|name| !canonical.contains(name))
            .collect();
        if extras.is_empty() {
            return order;
        }
        match order.pop() {
            Some(last) => {
                order.extend(extras);
                order.push(last);
            }
            None => order = extras,
        }
        order
    }

    /// Render one section to RST, exactly like the OCaml `get`/`to_rst` pair:
    /// the section header ([`make_header`]) followed by the body.
    ///
    /// Handler-owned modules are delegated to [`crate::handlers::Handler`];
    /// the generic ones are rendered from their field table.
    ///
    /// Missing attributes are *materialised* first: the generated `read_<f>`
    /// writes `EZFIO.cfg`'s `default:` whenever the key file is absent.  This
    /// is a write, and it is intentional in check mode too (contract C3: the
    /// editor may create defaulted attributes that were absent).
    ///
    /// A failure that OCaml raises as `Sys_error` — reading an attribute that
    /// is not set, or any other I/O failure — is not fatal: `get` prints
    /// `Info: <message>` to stderr and renders the section as `""`.  That is
    /// the `try ... with Sys_error` of `qp_edit.ml`.
    ///
    /// `ezfio` must have been opened on the calculation directory; the engine
    /// never opens or closes it.
    ///
    /// # Errors
    ///
    /// Fails for a non-section or unknown module, and for any non-`Sys_error`
    /// failure (a domain assertion in a handler, a malformed EZFIO value).
    pub fn get(&self, ezfio: &EzFio, name: &str) -> Result<String> {
        let module = self.section(name)?;
        let rendered = match &module.handler {
            Some(handler_name) => handler_of(module, handler_name)?.get(ezfio, &module.title),
            None => {
                let mut out = make_header(&module.title);
                out.push_str(&generic_body(ezfio, module)?);
                Ok(out)
            }
        };
        match rendered {
            Ok(text) => Ok(text),
            Err(error) if is_sys_error(&error) => {
                eprintln!("Info: {error}");
                Ok(String::new())
            }
            Err(error) => Err(error),
        }
    }

    /// Parse one section of `text` (the whole edited RST) and, when `apply`
    /// is true, write the values back to EZFIO.
    ///
    /// This is the OCaml `set ?apply`: the section body is located by its
    /// header and bounded by the next `"=="`, exactly as
    /// `String_ext.substr_index` does.  A section whose header is absent is a
    /// silent no-op (the user may have deleted it).
    ///
    /// Handler-owned modules are delegated to [`crate::handlers::Handler`];
    /// `ao_basis`, `ao_extra_basis` and `mo_basis` parse nothing because the
    /// OCaml `set` has a `()` arm for them.
    ///
    /// A conversion failure that sexplib reports as a `Failure` is printed to
    /// stderr as `Error:  (...)` and the section is skipped: the OCaml tool
    /// returns `None` from `of_rst` in that case.
    ///
    /// Any other failure — a domain check such as a negative `Threshold`
    /// (OCaml `Invalid_argument`), or a write error — takes the reference's
    /// retry path: `Info: Read error in <title>` is printed, the section is
    /// parsed a second time without writing, and only a failure of that
    /// second parse is returned (the driver turns it into the `=-` banner).
    /// The retry's own `Error:  (...)` is printed when it is a reported
    /// sexplib failure, exactly as the generated `of_rst` does.
    pub fn set(&self, ezfio: &EzFio, name: &str, text: &str, apply: bool) -> Result<()> {
        let module = self.section(name)?;
        let header = make_header(&module.title);
        let Some(idx) = text.find(&header) else {
            return Ok(());
        };
        let index_begin = idx + header.len();
        // `Str.search_forward` is called at this position; it raises
        // `Invalid_argument "Str.search_forward"` when the position is past
        // the end of the string, and `set` does not catch it.
        let search_from = index_begin + header.len() + 1;
        if search_from > text.len() {
            bail!("Str.search_forward");
        }
        let index_end = text[search_from..]
            .find("==")
            .map(|i| search_from + i)
            .unwrap_or(text.len());
        let body = &text[index_begin..index_end];

        if let Some(handler_name) = &module.handler {
            return handler_of(module, handler_name)?.set(ezfio, &module.title, body, apply);
        }

        let outcome = match parse_body(module, body) {
            Ok(values) => {
                if apply {
                    let dir = ezfio_dir(module)?;
                    for (field, value) in values {
                        let spec = spec_of(field)?;
                        write_value(ezfio, dir, field, spec, &value)?;
                    }
                }
                Ok(())
            }
            Err(RstFailure::Reported { message, snippet }) => {
                eprint!("Error:  ({message})\n\n  {snippet}\n\n");
                return Ok(());
            }
            Err(RstFailure::Fatal(error)) => Err(error),
        };
        if outcome.is_ok() {
            return Ok(());
        }
        // The reference's `try ... with _ ->` covers both `of_rst` and the
        // write, and retries with a parse-only call whose value is ignored.
        eprintln!("Info: Read error in {}", module.title);
        match parse_body(module, body) {
            Ok(_) => Ok(()),
            Err(RstFailure::Reported { message, snippet }) => {
                eprint!("Error:  ({message})\n\n  {snippet}\n\n");
                Ok(())
            }
            Err(RstFailure::Fatal(error)) => Err(error),
        }
    }
}

/// Resolve the handler of a module.
///
/// The loader already rejects an unknown handler name, so a failure here is
/// an internal inconsistency; it is still returned rather than panicked so
/// that a hand-written config cannot abort the editor.
fn handler_of(module: &ModuleConfig, name: &str) -> Result<crate::handlers::Handler> {
    crate::handlers::Handler::from_name(name).ok_or_else(|| {
        anyhow::anyhow!(
            "module `{}` names the unknown handler `{name}`",
            module.name
        )
    })
}

/// Whether `error` belongs to the OCaml `Sys_error` class.
///
/// `qp_edit.get` catches that class, prints `Info: <message>` and renders the
/// section as `""`.  It covers every failed `Sys.*` call, in practice a
/// missing attribute (`Ezfio.get_*` opens the file) and any other I/O error;
/// [`qp_ezfio::Error::NotSet`] is the same condition on the Rust side.
fn is_sys_error(error: &anyhow::Error) -> bool {
    matches!(
        error.downcast_ref::<qp_ezfio::Error>(),
        Some(qp_ezfio::Error::Io { .. } | qp_ezfio::Error::NotSet { .. })
    )
}

/// The file header written at the top of the temporary RST file
/// (`file_header` in `ocaml/qp_edit.ml`).
///
/// `filename` is printed verbatim between backticks.
pub fn file_header(filename: &str) -> String {
    format!(
        "\n==================================================================\n                       Quantum Package\n==================================================================\n\nEditing file `{filename}`\n\n"
    )
}

/// The header that delimits one section (`make_header` in
/// `ocaml/qp_edit.ml`): two blank lines, the title, a rule of `=` as long as
/// the title, and two blank lines.
///
/// The rule length is the title's **byte** length, as OCaml's
/// `String.length`.
pub fn make_header(title: &str) -> String {
    format!("\n\n{title}\n{}\n\n", "=".repeat(title.len()))
}

/// Fields present in the generated OCaml record: exactly those declaring the
/// `ocaml` interface.  `editable` is a different flag (it excludes the keys a
/// handler owns) and must not be used to choose the record shape.
fn ocaml_fields(module: &ModuleConfig) -> impl Iterator<Item = &Field> {
    module
        .fields
        .iter()
        .filter(|field| field.interfaces.iter().any(|i| i == "ocaml"))
}

/// The single EZFIO group the generic engine reads and writes.
fn ezfio_dir(module: &ModuleConfig) -> Result<&str> {
    match module.ezfio_dir.as_slice() {
        [dir] => Ok(dir),
        other => bail!(
            "module `{}` declares {} EZFIO directories; the generic engine needs exactly one",
            module.name,
            other.len()
        ),
    }
}

/// The type table entry of a field.
///
/// The loader rejects unknown type names, so a miss here is an internal
/// inconsistency rather than user input.
fn spec_of(field: &Field) -> Result<&'static TypeSpec> {
    lookup_type(&field.type_name).ok_or_else(|| {
        anyhow::anyhow!(
            "field `{}` declares the unknown type {:?}",
            field.provider,
            field.type_name
        )
    })
}

/// Canonical spelling of a choice value (`Disk_access.of_string` followed by
/// `to_string`).  Values outside the choice set are returned unchanged; the
/// domain check rejects them before this point on the parse path.
fn canonical_choice(spec: &TypeSpec, text: &str) -> String {
    if let Domain::Choices(choices) = spec.domain()
        && let Some(choice) = choices.iter().find(|c| c.eq_ignore_ascii_case(text))
    {
        return (*choice).to_string();
    }
    text.to_string()
}

/// Render one value the way the generated `to_rst` does, i.e. through the
/// fancy type's `to_string`.
fn render_value(spec: &TypeSpec, value: &Value) -> String {
    match value {
        Value::Int(v) => v.to_string(),
        Value::Float(x) => string_of_float(*x),
        Value::Bool(b) => if *b { "true" } else { "false" }.to_string(),
        Value::Str(s) => canonical_choice(spec, s),
    }
}

/// The RST body of a generic module: one block per `ocaml` field, in
/// `EZFIO.cfg` order.
///
/// Layout of one block (the literal in the generated `to_rst`):
/// `"\n   {doc} ::\n   \n     {ezfio_name} = {value}\n   "`.  Blocks are
/// concatenated, so the body ends with `"\n   "` and has no final newline.
fn generic_body(ezfio: &EzFio, module: &ModuleConfig) -> Result<String> {
    let dir = ezfio_dir(module)?;
    let mut out = String::new();
    for field in ocaml_fields(module) {
        let spec = spec_of(field)?;
        if !ezfio.has(dir, &field.ezfio_name)? {
            materialise_default(ezfio, module, dir, field, spec)?;
        }
        let value = read_value(ezfio, dir, field, spec)?;
        out.push_str(&format!(
            "\n   {} ::\n   \n     {} = {}\n   ",
            field.doc,
            field.ezfio_name,
            render_value(spec, &value)
        ));
    }
    Ok(out)
}

/// Write the default of an absent key, the way `read_<f>` does:
/// `get_default "<name>" |> <low-level of_string> |> Ezfio.set_<dir>_<name>`.
///
/// `get_default` is `Qpackage.get_ezfio_default "<module>" "<provider>"`
/// (`ocaml/Qpackage.ml`): the value is read at run time from
/// `$QP_ROOT/data/ezfio_defaults/*.ezfio_interface_default`, so a build that
/// ships a different data tree defaults differently.  Only when no such value
/// exists does the config's `EZFIO.cfg` default apply — a documented lenient
/// deviation from the reference, which raises `Not_found` (uncaught by
/// `qp_edit`) in that case.
///
/// A formula default cannot be evaluated by the editor; the loader guarantees
/// that only handler-owned fields have one, and handlers never reach here.
fn materialise_default(
    ezfio: &EzFio,
    module: &ModuleConfig,
    dir: &str,
    field: &Field,
    spec: &TypeSpec,
) -> Result<()> {
    if field.default_is_formula {
        bail!(
            "field `{}` has a Fortran-formula default that the editor cannot evaluate",
            field.provider
        );
    }
    let runtime = crate::defaults::lookup(&module.name, &field.provider);
    let default = runtime
        .as_deref()
        .or(field.default.as_deref())
        .ok_or_else(|| {
            anyhow::anyhow!(
                "field `{}` is absent from EZFIO and has no default",
                field.provider
            )
        })?;
    let value = spec.parse_scalar(default).with_context(|| {
        format!(
            "invalid default {:?} for field `{}`",
            default, field.provider
        )
    })?;
    write_value(ezfio, dir, field, spec, &value)
}

/// Read one scalar and validate it through the type table.
fn read_value(ezfio: &EzFio, dir: &str, field: &Field, spec: &TypeSpec) -> Result<Value> {
    let key = field.ezfio_name.as_str();
    let raw = match spec.value() {
        ValueKind::Int => ezfio.get::<i64>(dir, key)?.to_string(),
        ValueKind::Float => ezfio.get::<f64>(dir, key)?.to_string(),
        ValueKind::Bool => {
            if ezfio.get::<bool>(dir, key)? {
                "true".to_string()
            } else {
                "false".to_string()
            }
        }
        ValueKind::String => ezfio.get::<String>(dir, key)?,
    };
    spec.parse_scalar(&raw)
        .with_context(|| format!("EZFIO attribute {dir}/{key}"))
}

/// Write one validated scalar, applying the fancy type's `to_*` conversion
/// on the way (`write_<f>` in the generated code).
fn write_value(
    ezfio: &EzFio,
    dir: &str,
    field: &Field,
    spec: &TypeSpec,
    value: &Value,
) -> Result<()> {
    let key = field.ezfio_name.as_str();
    let result = match value {
        Value::Int(v) => ezfio.set::<i64>(dir, key, v),
        Value::Float(x) => ezfio.set::<f64>(dir, key, x),
        Value::Bool(b) => ezfio.set::<bool>(dir, key, b),
        Value::Str(s) => {
            let canonical = canonical_choice(spec, s);
            ezfio.set::<String>(dir, key, &canonical)
        }
    };
    result.with_context(|| format!("writing EZFIO attribute {dir}/{key}"))
}

/// Why parsing one section failed.
enum RstFailure {
    /// sexplib would report this as a `Failure` inside `of_rst`: the message
    /// and snippet are printed and the section is skipped.
    Reported {
        /// Text printed after `Error:  (`.
        message: String,
        /// `##`-marked excerpt, see [`crate::sexp::snippet`].
        snippet: String,
    },
    /// A conversion exception that `fail_msg` re-raises (`Invalid_argument`
    /// from a domain check): the driver prints `Info:` and retries.
    Fatal(anyhow::Error),
}

/// Parse a section body into `(field, value)` pairs in record order.
///
/// The filter is the generated `of_rst`: keep every line containing `=`,
/// replace `=` with a space, wrap each line in parentheses, concatenate and
/// wrap the result once more.  The resulting sexp is a list of
/// `(name value)` pairs, one per `ocaml` field.
fn parse_body<'a>(
    module: &'a ModuleConfig,
    body: &str,
) -> std::result::Result<Vec<(&'a Field, Value)>, RstFailure> {
    let names: Vec<&str> = ocaml_fields(module)
        .map(|field| field.ezfio_name.as_str())
        .collect();
    let mut source = String::from("(");
    for line in body.split('\n').filter(|line| line.contains('=')) {
        source.push('(');
        source.push_str(&line.replace('=', " "));
        source.push(')');
    }
    source.push(')');

    let node = sexp::parse(&source).map_err(|error| RstFailure::Reported {
        message: error.message,
        snippet: sexp::snippet(&source, error.start, error.end),
    })?;
    let pairs = sexp::record_fields(&node, &names).map_err(|error| RstFailure::Reported {
        message: error.message,
        snippet: sexp::snippet(&source, error.start, error.end),
    })?;

    let mut values = Vec::with_capacity(names.len());
    for field in ocaml_fields(module) {
        let Some((_, pair, value_node)) =
            pairs.iter().find(|(name, _, _)| *name == field.ezfio_name)
        else {
            // `record_fields` checked completeness; a miss here is a bug.
            return Err(RstFailure::Fatal(anyhow::anyhow!(
                "internal error: field `{}` disappeared from the parsed record",
                field.ezfio_name
            )));
        };
        values.push((field, convert_field(field, pair, value_node, &source)?));
    }
    Ok(values)
}

/// Convert one value node with the field's declared type.
///
/// The primitive parse is done first so that a syntax error can be reported
/// by `fail_msg` (OCaml `Failure`) while a domain violation takes the
/// exception path (OCaml `Invalid_argument` from the fancy `of_*`).
fn convert_field(
    field: &Field,
    pair: &sexp::Sexp,
    node: &sexp::Sexp,
    source: &str,
) -> std::result::Result<Value, RstFailure> {
    let spec = spec_of(field).map_err(RstFailure::Fatal)?;
    let Some(text) = node.atom() else {
        return Err(RstFailure::Reported {
            message: format!("{}_of_sexp: expected an atom", spec.value().name()),
            snippet: {
                let (start, end) = pair.span();
                sexp::snippet(source, start, end)
            },
        });
    };
    let syntax_ok = match spec.value() {
        ValueKind::Int => text.parse::<i64>().is_ok(),
        ValueKind::Float => text.parse::<f64>().is_ok(),
        ValueKind::Bool => text.eq_ignore_ascii_case("true") || text.eq_ignore_ascii_case("false"),
        ValueKind::String => true,
    };
    if !syntax_ok {
        let (start, end) = pair.span();
        return Err(RstFailure::Reported {
            message: format!(
                "{}_of_sexp: invalid {} {:?}",
                spec.value().name(),
                spec.value().name(),
                text
            ),
            snippet: sexp::snippet(source, start, end),
        });
    }
    spec.parse_scalar(text).map_err(RstFailure::Fatal)
}

/// Render one module table.
fn describe_module(module: &ModuleConfig) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "module `{}` ({})\n",
        module.name,
        if module.section {
            "section"
        } else {
            "derived, not visited"
        }
    ));
    out.push_str(&format!("  title:     {}\n", module.title));
    out.push_str(&format!(
        "  handler:   {}\n",
        module.handler.as_deref().unwrap_or("-")
    ));
    out.push_str(&format!("  source:    {}\n", module.source_cfg));
    out.push_str(&format!("  ezfio_dir: {}\n", module.ezfio_dir.join(", ")));
    let editable = module.fields.iter().filter(|field| field.editable).count();
    out.push_str(&format!(
        "  fields:    {} ({} editable)\n",
        module.fields.len(),
        editable
    ));

    if module.fields.is_empty() {
        return out;
    }

    let provider_width = column_width(module.fields.iter().map(|f| f.provider.as_str()), 8);
    let type_width = column_width(module.fields.iter().map(|f| f.type_name.as_str()), 4);
    let storage_cells: Vec<String> = module.fields.iter().map(storage_cell).collect();
    let storage_width = column_width(storage_cells.iter().map(String::as_str), 7);
    let default_cells: Vec<String> = module.fields.iter().map(default_cell).collect();
    let default_width = column_width(default_cells.iter().map(String::as_str), 7);

    out.push('\n');
    out.push_str(&format!(
        "  {:<provider_width$}  {:<type_width$}  {:<storage_width$}  {:<8}  {:<default_width$}  {}\n",
        "provider", "type", "storage", "editable", "default", "domain",
        provider_width = provider_width,
        type_width = type_width,
        storage_width = storage_width,
        default_width = default_width,
    ));
    for (index, field) in module.fields.iter().enumerate() {
        let domain = lookup_type(&field.type_name)
            .map(|spec| spec.describe_domain())
            .unwrap_or_else(|| "unknown".to_string());
        out.push_str(&format!(
            "  {:<provider_width$}  {:<type_width$}  {:<storage_width$}  {:<8}  {:<default_width$}  {}\n",
            field.provider,
            field.type_name,
            storage_cells[index],
            if field.editable { "yes" } else { "no" },
            default_cells[index],
            domain,
            provider_width = provider_width,
            type_width = type_width,
            storage_width = storage_width,
            default_width = default_width,
        ));
    }
    out
}

/// The storage column: `scalar`, or `array <size>` for rank-≥1 providers.
fn storage_cell(field: &Field) -> String {
    if field.is_array() {
        format!("array {}", field.size)
    } else {
        "scalar".to_string()
    }
}

/// The default column: `-` (none), `formula` (unusable expression), or the
/// literal value in TOML-style quotes.
fn default_cell(field: &Field) -> String {
    match &field.default {
        None => "-".to_string(),
        Some(_) if field.default_is_formula => "formula".to_string(),
        Some(value) => format!("{value:?}"),
    }
}

/// Width of a table column: the longest cell, with a floor so headers line up.
fn column_width<'a>(cells: impl Iterator<Item = &'a str>, floor: usize) -> usize {
    cells.map(str::len).max().unwrap_or(0).max(floor)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::loader::Config;
    use crate::schema::{Field, ModuleConfig, ModuleRef, Registry};
    use std::path::PathBuf;

    fn field(provider: &str, type_name: &str, editable: bool) -> Field {
        Field {
            provider: provider.to_string(),
            ezfio_name: provider.to_string(),
            type_name: type_name.to_string(),
            interfaces: vec!["ezfio".to_string()],
            doc: "doc".to_string(),
            size: "1".to_string(),
            default: None,
            default_is_formula: false,
            editable,
        }
    }

    fn config() -> Config {
        let module = ModuleConfig {
            schema_version: crate::SCHEMA_VERSION,
            name: "alpha".to_string(),
            title: "Alpha".to_string(),
            section: true,
            handler: None,
            source_cfg: "src/alpha/EZFIO.cfg".to_string(),
            ezfio_dir: vec!["alpha".to_string()],
            fields: vec![
                field("second", "integer", true),
                field("first", "Threshold", true),
            ],
        };
        let registry = Registry {
            schema_version: crate::SCHEMA_VERSION,
            modules: vec![ModuleRef {
                name: "alpha".to_string(),
                file: "alpha.toml".to_string(),
                section: true,
                handler: None,
            }],
        };
        Config::from_parts(PathBuf::from("/config"), registry, vec![module])
    }

    #[test]
    fn describe_is_deterministic_and_keeps_field_order() {
        let engine = Engine::new(config());
        let first = engine.describe("alpha").expect("module exists");
        let second = engine.describe("alpha").expect("module exists");
        assert_eq!(first, second);
        let second_at = first.find("second").expect("second present");
        let first_at = first.find("first").expect("first present");
        assert!(second_at < first_at, "fields must keep EZFIO.cfg order");
        assert!(first.contains("Threshold"));
        assert!(first.contains("float in [0, 1]"));
    }

    #[test]
    fn stats_count_fields_and_editability() {
        let engine = Engine::new(config());
        assert_eq!(
            engine.stats(),
            Stats {
                modules: 1,
                sections: 1,
                fields: 2,
                editable_fields: 2,
                formula_defaults: 0,
            }
        );
    }

    #[test]
    fn non_section_is_rejected_for_interactive_use() {
        let engine = Engine::new(config());
        assert!(engine.section("alpha").is_ok());
        assert!(engine.section("beta").is_err());
        assert!(engine.describe("beta").is_err());
    }

    // ------------------------------------------------------------------
    // Generic get/set
    // ------------------------------------------------------------------

    /// One `ocaml` field, with the given default.
    fn ocaml_field(provider: &str, type_name: &str, doc: &str, default: Option<&str>) -> Field {
        Field {
            provider: provider.to_string(),
            ezfio_name: provider.to_string(),
            type_name: type_name.to_string(),
            interfaces: vec![
                "ezfio".to_string(),
                "ocaml".to_string(),
                "provider".to_string(),
            ],
            doc: doc.to_string(),
            size: "1".to_string(),
            default: default.map(str::to_string),
            default_is_formula: false,
            editable: true,
        }
    }

    /// A one-module configuration whose fields cover all four value kinds,
    /// plus a choice type.
    fn value_config(handler: Option<&str>) -> Config {
        let module = ModuleConfig {
            schema_version: crate::SCHEMA_VERSION,
            name: "alpha".to_string(),
            title: "Alpha".to_string(),
            section: true,
            handler: handler.map(str::to_string),
            source_cfg: "src/alpha/EZFIO.cfg".to_string(),
            ezfio_dir: vec!["alpha".to_string()],
            fields: vec![
                ocaml_field("n", "integer", "Count", Some("7")),
                ocaml_field("x", "double precision", "Tolerance", Some("1.e-15")),
                ocaml_field("flag", "logical", "Enable", Some("true")),
                ocaml_field("mode", "Disk_access", "Disk mode", Some("read")),
            ],
        };
        let registry = Registry {
            schema_version: crate::SCHEMA_VERSION,
            modules: vec![ModuleRef {
                name: "alpha".to_string(),
                file: "alpha.toml".to_string(),
                section: true,
                handler: handler.map(str::to_string),
            }],
        };
        Config::from_parts(PathBuf::from("/config"), registry, vec![module])
    }

    #[test]
    fn section_headers_match_the_ocaml_layout() {
        assert_eq!(make_header("Mrci"), "\n\nMrci\n====\n\n");
        assert_eq!(
            file_header("qp_edit_x.rst"),
            "\n==================================================================\n                       Quantum Package\n==================================================================\n\nEditing file `qp_edit_x.rst`\n\n"
        );
    }

    #[test]
    fn generic_get_materialises_defaults_and_renders_every_kind() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        let engine = Engine::new(value_config(None));

        let text = engine.get(&ezfio, "alpha").expect("render");

        // Every absent key was written from `EZFIO.cfg`'s default.
        assert_eq!(ezfio.get::<i64>("alpha", "n").expect("n"), 7);
        assert_eq!(ezfio.get::<f64>("alpha", "x").expect("x"), 1e-15);
        assert!(ezfio.get::<bool>("alpha", "flag").expect("flag"));
        // `Disk_access.of_string "read" |> to_string` is the canonical
        // spelling, which is what gets written.
        assert_eq!(ezfio.get::<String>("alpha", "mode").expect("mode"), "Read");

        assert_eq!(
            text,
            concat!(
                "\n\nAlpha\n=====\n\n",
                "\n   Count ::\n   \n     n = 7\n   ",
                "\n   Tolerance ::\n   \n     x = 1e-15\n   ",
                "\n   Enable ::\n   \n     flag = true\n   ",
                "\n   Disk mode ::\n   \n     mode = Read\n   ",
            )
        );
    }

    /// Replace an exact fragment, failing loudly when it is absent: a test
    /// that silently edits nothing would pass for the wrong reason.
    fn replace_once(text: &str, from: &str, to: &str) -> String {
        assert!(text.contains(from), "fragment {from:?} not found");
        text.replacen(from, to, 1)
    }

    #[test]
    fn generic_set_parses_without_writing_in_check_mode() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        let engine = Engine::new(value_config(None));
        let text = engine.get(&ezfio, "alpha").expect("render");
        let edited = replace_once(&text, "     n = 7", "     n = 42");

        engine
            .set(&ezfio, "alpha", &edited, false)
            .expect("check-mode parse");
        assert_eq!(
            ezfio.get::<i64>("alpha", "n").expect("n"),
            7,
            "apply=false must not write"
        );
    }

    #[test]
    fn generic_set_applies_the_edited_record() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        let engine = Engine::new(value_config(None));
        let text = engine.get(&ezfio, "alpha").expect("render");
        let edited = replace_once(&text, "     n = 7", "     n = 42");

        engine.set(&ezfio, "alpha", &edited, true).expect("apply");
        assert_eq!(ezfio.get::<i64>("alpha", "n").expect("n"), 42);
        // Untouched fields are rewritten with their current values.
        assert_eq!(ezfio.get::<f64>("alpha", "x").expect("x"), 1e-15);
        assert!(ezfio.get::<bool>("alpha", "flag").expect("flag"));

        // Re-rendering shows the new value.
        let again = engine.get(&ezfio, "alpha").expect("render");
        assert!(again.contains("     n = 42"), "{again}");
        assert!(!again.contains("     n = 7\n"));
    }

    #[test]
    fn set_without_the_section_header_is_a_silent_no_op() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        let engine = Engine::new(value_config(None));
        engine
            .set(&ezfio, "alpha", "nothing to see\n", true)
            .expect("absent header");
        assert!(!ezfio.has("alpha", "n").expect("has"));
    }

    #[test]
    fn set_past_the_end_reproduces_str_search_forward() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        let engine = Engine::new(value_config(None));
        let error = engine
            .set(&ezfio, "alpha", &make_header("Alpha"), true)
            .expect_err("Str.search_forward");
        assert_eq!(error.to_string(), "Str.search_forward");
    }

    #[test]
    fn set_reports_a_malformed_record_and_writes_nothing() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        let engine = Engine::new(value_config(None));
        let text = engine.get(&ezfio, "alpha").expect("render");
        // Drop one field: sexplib reports a missing field and `of_rst`
        // returns None, so no value is applied.
        let edited = replace_once(&text, "\n   Tolerance ::\n   \n     x = 1e-15\n   ", "");
        engine
            .set(&ezfio, "alpha", &edited, true)
            .expect("reported failure is not fatal");
        assert_eq!(ezfio.get::<i64>("alpha", "n").expect("n"), 7);
    }

    #[test]
    fn set_raises_on_a_domain_violation() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        let mut config = value_config(None);
        // Rebuild the module with a range-checked integer type.
        let module = ModuleConfig {
            fields: vec![ocaml_field(
                "n",
                "Strictly_positive_int",
                "Count",
                Some("7"),
            )],
            ..config.modules()[0].clone()
        };
        config = Config::from_parts(
            PathBuf::from("/config"),
            Registry {
                schema_version: crate::SCHEMA_VERSION,
                modules: vec![ModuleRef {
                    name: "alpha".to_string(),
                    file: "alpha.toml".to_string(),
                    section: true,
                    handler: None,
                }],
            },
            vec![module],
        );
        let engine = Engine::new(config);
        let text = engine.get(&ezfio, "alpha").expect("render");
        let edited = replace_once(&text, "     n = 7", "     n = 0");

        let error = engine
            .set(&ezfio, "alpha", &edited, true)
            .expect_err("invalid_argument path");
        assert!(
            error.to_string().contains("Strictly_positive_int"),
            "{error}"
        );
        assert_eq!(ezfio.get::<i64>("alpha", "n").expect("n"), 7);
    }

    #[test]
    fn handler_modules_are_delegated_instead_of_read_as_tables() {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        // The module is named `alpha` but its handler is `electrons`: the
        // handler reads and writes the `electrons` group, and the generic
        // table (`alpha/n`, `alpha/x`, ...) is never touched.
        let engine = Engine::new(value_config(Some("electrons")));

        assert_eq!(engine.get(&ezfio, "alpha").expect("render"), "");
        assert!(!ezfio.has("alpha", "n").expect("has"));

        let text = concat!(
            "\n\nAlpha\n=====\n\n",
            "\nSpin multiplicity is Singlet.\n\nNumber of alpha and beta electrons ::\n\n",
            "  elec_alpha_num = 7\n  elec_beta_num  = 7\n\n",
        );
        engine.set(&ezfio, "alpha", text, true).expect("apply");
        assert_eq!(
            ezfio
                .get::<i64>("electrons", "elec_alpha_num")
                .expect("alpha"),
            7
        );
        assert_eq!(
            ezfio
                .get::<i64>("electrons", "elec_beta_num")
                .expect("beta"),
            7
        );
        assert!(!ezfio.has("alpha", "n").expect("has"));
    }

    // ------------------------------------------------------------------
    // Visit order
    // ------------------------------------------------------------------

    /// A configuration of named modules; `section = false` marks a derived
    /// module that `qp_edit` never visits.
    fn named_config(names: &[(&str, bool)]) -> Config {
        let modules: Vec<ModuleConfig> = names
            .iter()
            .map(|(name, section)| ModuleConfig {
                schema_version: crate::SCHEMA_VERSION,
                name: (*name).to_string(),
                title: (*name).to_string(),
                section: *section,
                handler: None,
                source_cfg: format!("src/{name}/EZFIO.cfg"),
                ezfio_dir: vec![(*name).to_string()],
                fields: Vec::new(),
            })
            .collect();
        let registry = Registry {
            schema_version: crate::SCHEMA_VERSION,
            modules: modules
                .iter()
                .map(|module| ModuleRef {
                    name: module.name.clone(),
                    file: format!("{}.toml", module.name),
                    section: module.section,
                    handler: None,
                })
                .collect(),
        };
        Config::from_parts(PathBuf::from("/config"), registry, modules)
    }

    /// A canonical name that the configuration lacks is skipped (a build
    /// without the `qmcchem` plugin has no `qmcchem` section), and a derived
    /// module is never visited even when the canonical list names it.
    #[test]
    fn visit_order_skips_absent_and_derived_modules() {
        let engine = Engine::new(named_config(&[
            ("aaa", true),
            ("ccc", true),
            ("bitmasks", false),
        ]));
        let canonical = ["aaa", "bbb", "ccc", "bitmasks"];
        assert_eq!(engine.visit_order(&canonical), vec!["aaa", "ccc"]);
    }

    /// A section the canonical list does not know is inserted in registry
    /// order before the last canonical section, which must stay last.
    #[test]
    fn visit_order_inserts_extra_sections_before_the_last_one() {
        let engine = Engine::new(named_config(&[
            ("aaa", true),
            ("xxx", true),
            ("yyy", true),
            ("zzz", true),
        ]));
        let canonical = ["aaa", "zzz"];
        assert_eq!(
            engine.visit_order(&canonical),
            vec!["aaa", "xxx", "yyy", "zzz"]
        );
    }

    /// The last *present* canonical section stays last even when names after
    /// it are absent: in the real list that section is
    /// `determinants_by_hand`, whose `of_rst` consumes the rest of the
    /// document.
    #[test]
    fn visit_order_keeps_the_last_present_canonical_section_last() {
        let engine = Engine::new(named_config(&[("aaa", true), ("xxx", true)]));
        let canonical = ["aaa", "bbb"];
        assert_eq!(engine.visit_order(&canonical), vec!["xxx", "aaa"]);
    }

    /// An unknown section that is *derived* in the configuration is ignored.
    #[test]
    fn visit_order_ignores_derived_extras() {
        let engine = Engine::new(named_config(&[("aaa", true), ("derived", false)]));
        let canonical = ["aaa"];
        assert_eq!(engine.visit_order(&canonical), vec!["aaa"]);
    }
}
