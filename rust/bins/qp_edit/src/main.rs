//! `qp_edit` — Rust port of `ocaml/qp_edit.ml` (P4, wave 3).
//!
//! The binary is config-driven per decision D3: the per-module TOML files
//! under `data/qp_edit_config/` are build artifacts (regenerated from the
//! `EZFIO.cfg` inputs, no recompilation).  The imperative behaviour of the
//! seven hand-written sections lives in `qp-input`'s handler modules; this
//! file is the driver, a line-by-line port of the OCaml `run`:
//!
//! 1. open the EZFIO directory and make sure `work/` exists;
//! 2. clear `work/qpstop` and `work/qpkill`;
//! 3. repair an unordered AO / AO-extra basis (the only EZFIO write that
//!    check mode allows, contract C2);
//! 4. apply `-n` (truncate the wavefunction) and `-s` (extract states);
//! 5. render every section to one RST document, in `tasks` order;
//! 6. run `$EDITOR` on it (interactive runs only — there is no line REPL;
//!    "interactive" means the external editor);
//! 7. read the document back and `set` every section, with `apply` false in
//!    check mode, so malformed input is still reported but nothing is
//!    written (contract C7).
//!
//! Non-check runs wrap steps 4–7 in a `work/backup.tar` backup, restored
//! when the run fails.  Failures print the reference's 80-column `=-` banner:
//! `Failure`/`Invalid_argument` classes exit 1, an OCaml `Assert_failure`
//! exits 2, exactly as the reference's top-level handler.
//!
//! Units: none.  Values are rendered by the OCaml-compatible printers in
//! `qp_input::ocaml_compat` (`%.12g` floats, `Range` strings, `MO_label`).

mod inspect;

use std::ffi::OsString;
use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

use anyhow::{Context, Result, bail};
use clap::Parser;
use qp_ezfio::EzFio;
use qp_input::engine::file_header;
use qp_input::handlers::{extract_states, reorder_ao_basis, reorder_ao_extra_basis, update_ndet};
use qp_input::{Engine, loader, ocaml_compat};
use qp_types::command_line::{EXIT_FAILURE, EXIT_USAGE, EditArgs, EditConfig};

/// The sections the editor visits, in the reference's `tasks` order
/// (`ocaml/qp_edit.ml`, post-fix tree).
///
/// The order is load-bearing.  The document is a concatenation, and
/// `determinants_by_hand` must come **last**: its `of_rst` locates its own
/// section and then consumes the rest of the document, which is only the
/// determinant body when nothing follows it.  The registry order is
/// alphabetical and must never be used here.
///
/// `qmcchem` is a *plugin* section: it is listed here because the reference
/// lists it, but only a configuration built with the plugin has it.  The run
/// filters this list against the loaded configuration and inserts any section
/// the list does not know, see [`Engine::visit_order`].
const TASKS: [&str; 38] = [
    "nuclei_by_hand",
    "ao_basis",
    "ao_extra_basis",
    "electrons",
    "mo_two_e_ints",
    "two_body_rdm",
    "ao_two_e_ints",
    "davidson_keywords",
    "mrci",
    "scf_utils",
    "mo_one_e_ints",
    "pseudo",
    "ezfio_files",
    "dressing",
    "utils",
    "extra_nuclei",
    "ao_one_e_ints",
    "becke_numerical_grid",
    "symmetry",
    "perturbation",
    "bitmask",
    "mu_of_r",
    "utils_trust_region",
    "mo_optimization_utils",
    "qmcchem",
    "hamiltonian",
    "density_for_dft",
    "nuclei",
    "casscf_cipsi",
    "dft_keywords",
    "cipsi",
    "basis",
    "trexio",
    "determinants",
    "utils_cc",
    "mol_properties",
    "mo_basis",
    "determinants_by_hand",
];

/// The `=-` rule `qp_edit.ml` prints around every fatal message (40 copies of
/// `=-`, i.e. 80 columns).
const BANNER: &str =
    "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-";

/// Exit code of an OCaml `Assert_failure` (`qp_edit.ml` exits 2; every other
/// caught class exits 1).
const EXIT_ASSERT: u8 = 2;

/// The reference help text, with `{program}` standing for `argv[0]`.
///
/// Captured from the reference binary (`qp_edit -h`, 717 bytes for
/// `argv[0] = "qp_edit"`).  The OCaml `Format` printer re-wraps the usage
/// line for a longer `argv[0]`; this port keeps the reference layout and
/// prints `argv[0]` verbatim, a documented deviation (no test or golden
/// transcript exercises `-h`).
const HELP: &str = r#"{program} - quantum_package command 

Usage:
  
  {program} [-c] [-h] [-n <int>] [-r <string>] [-s <range>] [-w <string>] 
      [--] EZFIO_DIR
  
Arguments:
  
  EZFIO_DIR             EZFIO directory 
  
Options:
  
  -c  --check           Checks the input data 
  -h  --help            Prints the help message. 
  -n  --ndet=<int>      Truncates the wavefunction to the target number 
                        of determinants 
  -r  --read=<string>   Reads the file and applies it to the EZFIO" 
  -s  --state=<range>   Extracts selected states, for example "[1,3-5]" 
  -w  --write=<string>  Writes the qp_edit file to a file" 
  
Description:
  
  Opens a text editor to edit the parameters of an EZFIO directory. 
  
"#;

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().collect();
    if is_inspection(&argv) {
        return inspect::run(argv);
    }

    // `EditArgs` reproduces the reference option grammar.  clap prints its
    // own diagnostics for an unknown option and exits 2, which is the
    // reference's usage exit code; the text differs and is documented.
    let cli = EditArgs::parse_from(argv.clone());
    let program = argv
        .first()
        .map(|arg| arg.to_string_lossy().into_owned())
        .unwrap_or_else(|| String::from("qp_edit"));

    // The reference `-h` prints the help in `set_specs` and then *continues*
    // with the normal argument checks; it does not exit.
    if cli.help {
        print!("{}", help_text(&program));
    }

    let config = match cli.validate() {
        Ok(config) => config,
        Err(error) => {
            // The reference calls `Command_line.help ()` again before
            // failing on a missing anonymous argument (measured: `-h`
            // prints the 717-byte help twice), and fails the `-n`/`-s`
            // conversions before that check, without help.
            if error.message() == "EZFIO_DIR is missing" {
                print!("{}", help_text(&program));
            }
            eprintln!("Fatal error: exception Failure({:?})", error.message());
            return ExitCode::from(EXIT_USAGE as u8);
        }
    };

    match execute(&config) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            print_banner(&format!("{error:#}"));
            if ocaml_compat::is_assert_failure(&error) {
                ExitCode::from(EXIT_ASSERT)
            } else {
                ExitCode::from(EXIT_FAILURE as u8)
            }
        }
    }
}

/// Whether `argv` asks for one of the maintenance inspection modes.
///
/// The reference has no long options of its own beyond the six of `qp_edit`,
/// so these spellings cannot collide with a valid editor command line.
fn is_inspection(argv: &[OsString]) -> bool {
    argv.iter().skip(1).any(|arg| {
        matches!(
            arg.to_str(),
            Some("--list" | "--describe" | "--check-config")
        )
    })
}

/// The help text for `program` (see [`HELP`] for the layout caveat).
fn help_text(program: &str) -> String {
    HELP.replace("{program}", program)
}

/// Print the reference's banner-wrapped fatal message on stderr.
fn print_banner(message: &str) {
    eprint!("{BANNER}\n\n{message}\n\n{BANNER}\n\n");
}

/// Load the configuration, create the backup when needed, then edit.
///
/// The backup is created *before* `edit` and only outside check mode, in the
/// reference's order.  A successful run removes `work/backup.tar`; a failed
/// one restores it with `tar`.  Both operations swallow their errors like the
/// reference.  Check mode neither creates nor removes a backup: the
/// reference's `at_exit` would delete a pre-existing `work/backup.tar`, and
/// dropping that is the C2-clean choice.
fn execute(config: &EditConfig) -> Result<()> {
    let engine = load_engine()?;
    let directory = Path::new(&config.ezfio_dir);
    if config.check {
        return edit(&engine, config);
    }
    create_backup(directory);
    let result = edit(&engine, config);
    if result.is_ok() {
        remove_backup(directory);
    } else {
        restore_backup(directory);
    }
    result
}

/// Discover and load the generated configuration.
///
/// Discovery is `$QP_EDIT_CONFIG` > `<exe_dir>/../data/qp_edit_config`; the
/// artifact is gitignored, so a missing one fails here with the loader's
/// message rather than silently editing nothing.
fn load_engine() -> Result<Engine> {
    let executable = std::env::current_exe().context("locating the qp_edit executable")?;
    let executable_dir = executable
        .parent()
        .context("the qp_edit executable has no parent directory")?;
    let env_config = std::env::var_os("QP_EDIT_CONFIG");
    let directory = loader::discover(None, env_config.as_deref(), executable_dir)?;
    let config = loader::load(&directory)
        .with_context(|| format!("loading configuration from {}", directory.display()))?;
    Ok(Engine::new(config))
}

/// `run` in `qp_edit.ml`: everything between the backup and the exit code.
///
/// All EZFIO access goes through one [`EzFio`] handle opened here; the
/// engine and the handlers never open the directory themselves.
fn edit(engine: &Engine, config: &EditConfig) -> Result<()> {
    // The editor runs only for a plain interactive invocation.  Any of
    // `-n`, `-s`, `-r`, `-w` or `-c` disables it (the reference's
    // `open_editor`).
    let open_editor = config.ndet.is_none()
        && config.state.is_none()
        && config.read.is_none()
        && config.write.is_none()
        && !config.check;

    if !Path::new(&config.ezfio_dir).exists() {
        // Reference message, typo included.
        bail!("{} does not exists", config.ezfio_dir);
    }
    // The reference only calls `Ezfio.set_file`; requiring a `.version`
    // here is stricter (a non-EZFIO directory now fails with the banner
    // instead of being mostly a no-op) and is documented.
    let ezfio = EzFio::open(&config.ezfio_dir)?;

    ensure_work_dir(&ezfio)?;
    for name in ["qpstop", "qpkill"] {
        let stopfile = ezfio.root().join("work").join(name);
        if stopfile.exists() {
            fs::remove_file(&stopfile)
                .with_context(|| format!("removing {}", stopfile.display()))?;
        }
    }

    // Contract C2: the only write check mode allows is repairing an
    // unordered basis.  Both helpers return without writing when the
    // ordering is already sorted.
    reorder_ao_basis(&ezfio)?;
    reorder_ao_extra_basis(&ezfio)?;

    if let Some(n_det) = config.ndet {
        update_ndet(&ezfio, n_det)?;
    }
    if let Some(states) = config.state.as_deref() {
        extract_states(&ezfio, states)?;
    }

    let document = prepare_document(engine, &ezfio, config)?;
    if open_editor {
        launch_editor(document.path());
    }
    if config.write.is_none() {
        read_back(engine, &ezfio, config, document.path())?;
    }
    Ok(())
}

/// Where the RST document of this run lives.
enum Document {
    /// `Filename.temp_file "qp_edit_" ".rst"`: removed when this value is
    /// dropped, like the reference's `at_exit (Sys.remove ...)`.
    Temporary(tempfile::TempPath),
    /// A user file (`-r` or `-w`), kept after the run.
    Permanent(PathBuf),
}

impl Document {
    /// The path the editor opens and the driver reads back.
    fn path(&self) -> &Path {
        match self {
            Document::Temporary(path) => path,
            Document::Permanent(path) => path,
        }
    }
}

/// Build the document exactly like `create_temp_file`:
/// `file_header :: List.map get tasks |> String.concat "\n"`, written with
/// `%s` (no extra final newline).
///
/// `-r` supplies the document, so nothing is rendered (and no getter, hence
/// no default materialisation, runs); `-r` with `-w` is rejected before any
/// rendering, as in the reference.
fn prepare_document(engine: &Engine, ezfio: &EzFio, config: &EditConfig) -> Result<Document> {
    match (&config.read, &config.write) {
        (Some(_), Some(_)) => bail!("read and write options are incompatible"),
        (Some(path), None) => Ok(Document::Permanent(PathBuf::from(path))),
        (None, Some(path)) => {
            let text = render_document(engine, ezfio, &config.ezfio_dir)?;
            fs::write(path, text).with_context(|| format!("writing {path}"))?;
            Ok(Document::Permanent(PathBuf::from(path)))
        }
        (None, None) => {
            let text = render_document(engine, ezfio, &config.ezfio_dir)?;
            let mut file = tempfile::Builder::new()
                .prefix("qp_edit_")
                .suffix(".rst")
                .tempfile()
                .context("creating the temporary RST file")?;
            file.write_all(text.as_bytes())
                .context("writing the temporary RST file")?;
            file.flush().context("flushing the temporary RST file")?;
            Ok(Document::Temporary(file.into_temp_path()))
        }
    }
}

/// Render every visited section, one getter call each, in the order
/// [`Engine::visit_order`] derives from [`TASKS`] and the loaded
/// configuration.
fn render_document(engine: &Engine, ezfio: &EzFio, filename: &str) -> Result<String> {
    let mut text = file_header(filename);
    for name in engine.visit_order(&TASKS) {
        text.push('\n');
        text.push_str(&engine.get(ezfio, name)?);
    }
    Ok(text)
}

/// Run the external editor.  `$EDITOR` defaults to `vi`; the command is
/// interpreted by `sh`, like OCaml's `Sys.command`, and its status is
/// ignored (the reference does not check it either).
fn launch_editor(path: &Path) {
    let editor = std::env::var("EDITOR").unwrap_or_else(|_| String::from("vi"));
    let command = format!("{editor} {}", path.display());
    let _ = Command::new("sh").arg("-c").arg(command).status();
}

/// Read the document back and apply every section.
///
/// `input_lines` + `String.concat "\n"` drops exactly one final newline, so
/// the body handed to `set` never ends with one.  `apply` is false in check
/// mode: the sections are still parsed (contract C7), nothing is written.
fn read_back(engine: &Engine, ezfio: &EzFio, config: &EditConfig, path: &Path) -> Result<()> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let text = text.strip_suffix('\n').unwrap_or(&text);
    // The same order as the rendering: every getter of this run has run
    // before, so the configuration has not changed under us.
    for name in engine.visit_order(&TASKS) {
        engine.set(ezfio, name, text, !config.check)?;
    }
    Ok(())
}

/// `Qpackage.ezfio_work`: create `work/` when it does not exist.
///
/// The reference runs `Ezfio.set_work_empty false` through the OCaml writer,
/// which produces exactly `F\n\n` (the boolean printer appends the newline of
/// `fortran_string_of_bool` *and* the scalar-file newline).  That file is
/// what marks the work directory as not empty; writing it directly avoids
/// the byte difference of the typed setter (`F\n`).
fn ensure_work_dir(ezfio: &EzFio) -> Result<()> {
    let work = ezfio.root().join("work");
    if work.exists() {
        return Ok(());
    }
    ezfio.create_group("work")?;
    let empty = work.join("empty");
    fs::write(&empty, "F\n\n").with_context(|| format!("writing {}", empty.display()))?;
    Ok(())
}

/// `remove_backup`: delete `<dir>/work/backup.tar`, ignoring a missing file.
fn remove_backup(directory: &Path) {
    let _ = fs::remove_file(directory.join("work").join("backup.tar"));
}

/// `create_backup`: `tar` the directory into `<dir>/work/backup.tar`.
///
/// The reference builds `.backup.tar` in the current directory and moves it
/// into place, with `--exclude="work/*"` so the archive never contains
/// itself.  The port builds the archive in a temporary file outside the tree
/// and copies it, which avoids leaving a `.backup.tar` behind in the user's
/// working directory when `mv` fails.  Every error is swallowed: a backup
/// that cannot be made must not prevent the edit (the reference's
/// `try ... with _ -> ()`).
fn create_backup(directory: &Path) {
    remove_backup(directory);
    let work = directory.join("work");
    let _ = fs::create_dir_all(&work);
    let Ok(temporary) = tempfile::Builder::new()
        .prefix("qp_edit_backup_")
        .suffix(".tar")
        .tempfile()
    else {
        return;
    };
    let status = Command::new("tar")
        .arg("-cf")
        .arg(temporary.path())
        .arg("--exclude=work/*")
        .arg(directory)
        .status();
    if matches!(status, Ok(status) if status.success()) {
        let _ = fs::copy(temporary.path(), work.join("backup.tar"));
    }
}

/// `restore_backup`: extract `<dir>/work/backup.tar` from the current
/// directory and remove it.  Errors are swallowed, as in the reference.
fn restore_backup(directory: &Path) {
    let backup = directory.join("work").join("backup.tar");
    if backup.exists() {
        let _ = Command::new("tar").arg("-xf").arg(&backup).status();
        remove_backup(directory);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tasks_are_the_reference_order() {
        assert_eq!(TASKS.len(), 38);
        assert_eq!(TASKS[0], "nuclei_by_hand");
        assert_eq!(TASKS[37], "determinants_by_hand");
        // `Qmcchem` is the plugin section of the reference's `tasks` list; it
        // sits between `Mo_optimization_utils` and `Hamiltonian`.
        assert_eq!(TASKS.iter().position(|name| *name == "qmcchem"), Some(24));
        assert_eq!(TASKS[23], "mo_optimization_utils");
        assert_eq!(TASKS[25], "hamiltonian");
        let mut sorted = TASKS.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), TASKS.len(), "duplicate task");
        assert!(!TASKS.contains(&"bitmasks"), "bitmasks is not a section");
        // The six handler-owned sections, `bitmasks` excluded.
        for handler in [
            "ao_basis",
            "ao_extra_basis",
            "mo_basis",
            "electrons",
            "nuclei_by_hand",
            "determinants_by_hand",
        ] {
            assert!(TASKS.contains(&handler), "missing handler {handler}");
        }
        // `determinants_by_hand` must be last: its `of_rst` consumes the
        // rest of the document.
        assert_eq!(
            TASKS
                .iter()
                .filter(|name| **name == "determinants_by_hand")
                .count(),
            1
        );
    }

    #[test]
    fn banner_is_the_reference_rule() {
        assert_eq!(BANNER.len(), 80);
        assert!(BANNER.starts_with("=-"));
        assert!(BANNER.ends_with("=-"));
        assert!(BANNER.bytes().all(|b| b == b'=' || b == b'-'));
    }

    #[test]
    fn help_text_has_the_reference_byte_count() {
        // 717 bytes is the measured `qp_edit -h` output with
        // `argv[0] = "qp_edit"` (the reference prints it twice when
        // EZFIO_DIR is missing, 1434 bytes).
        assert_eq!(help_text("qp_edit").len(), 717);
        assert!(help_text("qp_edit").contains("qp_edit [-c] [-h]"));
    }

    #[test]
    fn only_the_three_maintenance_modes_are_dispatched() {
        let argv = |args: &[&str]| -> Vec<OsString> {
            std::iter::once("qp_edit")
                .chain(args.iter().copied())
                .map(OsString::from)
                .collect()
        };
        assert!(is_inspection(&argv(&["--list"])));
        assert!(is_inspection(&argv(&["--config", "d", "--check-config"])));
        assert!(is_inspection(&argv(&["--describe", "--module", "x"])));
        assert!(!is_inspection(&argv(&["-c", "ezfio"])));
        assert!(!is_inspection(&argv(&["--config", "d"])));
        assert!(!is_inspection(&argv(&["--version"])));
    }
}
