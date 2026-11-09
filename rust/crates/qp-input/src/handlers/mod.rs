//! The six sections whose edit logic is hand-written.
//!
//! Decision D3 keeps the generic engine declarative and moves the behaviours
//! that cannot be expressed as a field table into named handlers ported from
//! the tracked OCaml files:
//!
//! | handler | OCaml source | section |
//! |---|---|---|
//! | `ao_basis` | `ocaml/Input_ao_basis.ml` | yes |
//! | `ao_extra_basis` | `ocaml/Input_ao_extra_basis.ml` | yes |
//! | `mo_basis` | `ocaml/Input_mo_basis.ml` | yes |
//! | `electrons` | `ocaml/Input_electrons.ml` | yes |
//! | `nuclei_by_hand` | `ocaml/Input_nuclei_by_hand.ml` | yes |
//! | `determinants_by_hand` | `ocaml/Input_determinants_by_hand.ml` | yes |
//!
//! The names here are the contract between
//! `scripts/ezfio_interface/qp_edit_config.py` (which writes `handler = ...`)
//! and the binary.  The loader rejects an unknown name and a name that does
//! not match the module id it is attached to.
//!
//! Status: P4 ports the six handlers that `qp_edit` dispatches to; see the
//! module table above for the file that holds each one.

use anyhow::{Result, anyhow, bail};
use qp_ezfio::{Array, EzFio};

pub mod ao_basis;
mod ao_extra_basis;
mod determinants_by_hand;
mod electrons;
pub mod mo_basis;
mod nuclei_by_hand;
mod qptypes;

/// One handler's `of_rst` result (wave S2.4: the OCaml `string -> t option`
/// plus its exception class, without the `Error:  (...)` print side effect).
///
/// `Parsed` is the OCaml `Some data`.  `Reported` is a sexplib `Failure`
/// (a value the fancy type rejects): the driver reports it as `qp_edit:
/// parse error in <title>: <message>` and skips the section, so the check
/// still exits 0.  `Escaped` is an exception that escaped `of_rst` (a fancy
/// `of_int` rejecting a value, or an EZFIO read failing inside the
/// converter): the driver retries once, silently, and only a second escape
/// propagates, wrapped with the section title.
#[derive(Debug)]
pub(crate) enum ParseResult<T> {
    /// The OCaml `Some data`: parsed, applied when writing back.
    Parsed(T),
    /// A reported sexplib `Failure`; `message` is the conversion error text.
    Reported(String),
    /// An exception that escaped `of_rst`.
    Escaped(anyhow::Error),
}

/// The OCaml `Ezfio.ezfio_array_of_list`: build an array from a flat list of
/// elements in Fortran (column-major) order.
///
/// `dims[0]` varies fastest, as in the reference
/// (`flattened_ezfio` reads `dim.(0)` innermost).
///
/// # Errors
///
/// The two ways the reference can fail, with its exact messages:
///
/// * `dims` empty — the reference's `read_nd` hits its `m < 1` guard and
///   raises `Failure "dimension should not be <1"`;
/// * fewer elements than `dims.iter().product()` — `Failure "Array is not
///   consistent"`.
///
/// Extra elements beyond the product are **ignored**, because the reference's
/// `read_1d`/`read_nd` loops stop after the last extent and drop the rest of
/// the list.  Both are plain `Failure`s in OCaml, so they are ordinary
/// `anyhow` errors here (exit code 1, not the assert code 2).
pub(crate) fn ezfio_array_of_list<T>(dims: Vec<usize>, mut data: Vec<T>) -> Result<Array<T>> {
    if dims.is_empty() {
        bail!("dimension should not be <1");
    }
    let expected = dims
        .iter()
        .try_fold(1usize, |acc, &dim| acc.checked_mul(dim))
        .unwrap_or(usize::MAX);
    if data.len() < expected {
        bail!("Array is not consistent");
    }
    data.truncate(expected);
    Ok(Array::new(dims, data)?)
}

/// A negative count is a corrupt file, not a panic.
///
/// The four basis/section handlers used to carry an identical private copy;
/// they now share this one.
pub(crate) fn usize_of(value: i64) -> Result<usize> {
    usize::try_from(value).map_err(|_| anyhow!("negative dimension or index: {value}"))
}

/// Reproduce the write wrapper of `qp_edit.set` for one section.
///
/// OCaml:
///
/// ```text
/// try
///   match of_rst str with
///   | Some data -> if apply then w data
///   | None -> ()
/// with _ ->
///   ignore (of_rst str)
/// ```
///
/// Wave S2.4 wording: a reported sexplib failure prints `qp_edit: parse
/// error in <title>: <message>` (the reference's `Error:  (...)` + snippet
/// is gone) and the section is skipped; the reference's
/// `Info: Read error in <title>` line is gone too — the retry is silent,
/// and only a second escape propagates, wrapped with the section title so
/// the driver's `qp_edit: {error:#}` reads `qp_edit: <title>: <message>`.
///
/// The `try` covers the write too, so a write failure is swallowed when the
/// retry parses.  The retry's value is never written; only a second escape
/// escapes.
pub(crate) fn apply_parsed<T>(
    title: &str,
    body: &str,
    write_back: bool,
    parse: impl Fn(&str) -> ParseResult<T>,
    write: impl Fn(&T) -> Result<()>,
) -> Result<()> {
    let outcome = match parse(body) {
        ParseResult::Parsed(data) if write_back => write(&data),
        ParseResult::Parsed(_) => Ok(()),
        ParseResult::Reported(message) => {
            eprintln!("qp_edit: parse error in {title}: {message}");
            Ok(())
        }
        ParseResult::Escaped(error) => Err(error),
    };
    if outcome.is_ok() {
        return Ok(());
    }
    match parse(body) {
        ParseResult::Parsed(_) => Ok(()),
        ParseResult::Reported(message) => {
            eprintln!("qp_edit: parse error in {title}: {message}");
            Ok(())
        }
        ParseResult::Escaped(error) => Err(error.context(title.to_string())),
    }
}

/// Fixtures shared by the handler tests.
///
/// The handlers are tested against copies of the `hcn` golden fixture, never
/// against the fixture itself: the read paths materialise missing defaults
/// (contract C3), and a test must not edit an input fixture.
#[cfg(test)]
pub(crate) mod test_support {
    use qp_ezfio::EzFio;
    use std::path::{Path, PathBuf};

    /// Recursively copy the directory `from` into `to`.
    pub(crate) fn copy_tree(from: &Path, to: &Path) {
        std::fs::create_dir_all(to).expect("create the destination directory");
        for entry in std::fs::read_dir(from).expect("read the fixture directory") {
            let entry = entry.expect("read one fixture entry");
            let target = to.join(entry.file_name());
            if entry.file_type().expect("fixture entry type").is_dir() {
                copy_tree(&entry.path(), &target);
            } else {
                std::fs::copy(entry.path(), &target).expect("copy a fixture file");
            }
        }
    }

    /// A writable copy of the `hcn` fixture, opened as an EZFIO.
    ///
    /// The temporary directory is returned so that it outlives the handle.
    pub(crate) fn hcn_ezfio() -> (tempfile::TempDir, EzFio) {
        let source =
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/fixtures/hcn");
        let dir = tempfile::tempdir().expect("create a temporary directory");
        let root = dir.path().join("ezfio");
        copy_tree(&source, &root);
        let ezfio = EzFio::open(&root).expect("open the hcn fixture");
        (dir, ezfio)
    }
}

/// A hand-written section handler.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Handler {
    /// AO basis ordering and MD5 consistency.
    AoBasis,
    /// Extra (auxiliary) AO basis ordering.
    AoExtraBasis,
    /// MO basis ordering and MD5 comparator.
    MoBasis,
    /// Electron counts, spin, and derived occupation.
    Electrons,
    /// Molecule definition parsed from RST.
    NucleiByHand,
    /// Determinant list parsed from RST and sexp.
    DeterminantsByHand,
}

impl Handler {
    /// Every handler, in generator-table order (the order has no meaning).
    pub const ALL: &'static [Handler] = &[
        Handler::AoBasis,
        Handler::AoExtraBasis,
        Handler::MoBasis,
        Handler::Electrons,
        Handler::NucleiByHand,
        Handler::DeterminantsByHand,
    ];

    /// The name written in `handler = "..."` and in the registry.
    pub const fn name(self) -> &'static str {
        match self {
            Handler::AoBasis => "ao_basis",
            Handler::AoExtraBasis => "ao_extra_basis",
            Handler::MoBasis => "mo_basis",
            Handler::Electrons => "electrons",
            Handler::NucleiByHand => "nuclei_by_hand",
            Handler::DeterminantsByHand => "determinants_by_hand",
        }
    }

    /// Whether `qp_edit` visits the module as an editable section.
    pub const fn is_section(self) -> bool {
        true
    }

    /// Whether check mode (`-c`) applies this section's parsed record to the
    /// EZFIO.
    ///
    /// Contract C2: the default is that nothing is written in check mode, so
    /// the flag is `false` everywhere except one override — user decision
    /// w12: the nuclei section is written back in check mode so the file is a
    /// fixed point of the master-format print, "for reproducibility"; master
    /// achieved this with the pre-15e725aa check mode.
    pub const fn applies_in_check(self) -> bool {
        matches!(self, Handler::NucleiByHand)
    }

    /// Parse a handler name from the TOML config.
    pub fn from_name(name: &str) -> Option<Handler> {
        Handler::ALL.iter().copied().find(|h| h.name() == name)
    }

    /// Render the section (`get` in `qp_edit.ml`): the section header plus
    /// the handler's `to_rst`, or `""` when the handler's `read` finds the
    /// section absent.
    ///
    /// `title` is [`crate::schema::ModuleConfig::title`], which is the
    /// `keyword_to_string` spelling used for the header.
    pub(crate) fn get(self, ezfio: &EzFio, title: &str) -> Result<String> {
        match self {
            Handler::AoBasis => ao_basis::get(ezfio, title),
            Handler::AoExtraBasis => ao_extra_basis::get(ezfio, title),
            Handler::MoBasis => mo_basis::get(ezfio, title),
            Handler::Electrons => electrons::get(ezfio, title),
            Handler::NucleiByHand => nuclei_by_hand::get(ezfio, title),
            Handler::DeterminantsByHand => determinants_by_hand::get(ezfio, title),
        }
    }

    /// Parse the section body and, when `write_back` is true, apply it
    /// (`set ?apply` in `qp_edit.ml`).
    ///
    /// `ao_basis`, `ao_extra_basis` and `mo_basis` are no-ops: the OCaml
    /// `set` has `()` for those arms, so their RST is never read back.
    pub(crate) fn set(
        self,
        ezfio: &EzFio,
        title: &str,
        body: &str,
        write_back: bool,
    ) -> Result<()> {
        match self {
            Handler::Electrons => electrons::set(ezfio, title, body, write_back),
            Handler::NucleiByHand => nuclei_by_hand::set(ezfio, title, body, write_back),
            Handler::DeterminantsByHand => {
                determinants_by_hand::set(ezfio, title, body, write_back)
            }
            Handler::AoBasis | Handler::AoExtraBasis | Handler::MoBasis => Ok(()),
        }
    }
}

/// `qp_edit -n`: keep the first `n_det` determinants of every state.
///
/// Port of `Input_determinants_by_hand.update_ndet`: the reference prints
/// `Reducing n_det to <n_det>` on stdout before doing anything, then rewrites
/// the determinants group with `force:true`.
///
/// # Errors
///
/// `"No determinants in file"` when the group is absent (the OCaml `failwith`),
/// and any EZFIO error of the read/write pair.
pub fn update_ndet(ezfio: &EzFio, n_det: i64) -> Result<()> {
    // `qp_edit.ml` passes `Det_number.of_int n`: an assertion (exit 2) below
    // one, an `Invalid_argument "Det_number: n"` (exit 1) above the stored
    // count.
    let n_det = qptypes::det_number_of_int(ezfio, n_det)?;
    determinants_by_hand::update_ndet(ezfio, n_det)
}

/// `qp_edit -s`: keep only the requested states and copy them to the front.
///
/// Port of `Input_determinants_by_hand.extract_states`.  `states` must already
/// be sorted and deduplicated (OCaml `Range.to_int_list` does that in the
/// command-line layer, which is [`qp_types::command_line::EditConfig`] here).
/// The reference prints `Extracting states <range>` and one `OK` per state
/// that actually needed shifting.
///
/// # Errors
///
/// `"No determinants in file"` when the group is absent, and
/// `"State to extract should not be greater than n_states"` when a requested
/// state exceeds the stored count — both OCaml `failwith`s.
pub fn extract_states(ezfio: &EzFio, states: &[i64]) -> Result<()> {
    determinants_by_hand::extract_states(ezfio, states)
}

/// Reorder the AO basis into canonical order, then the MOs with it.
///
/// Port of the first block of `qp_edit.ml`'s `run`.  When
/// [`ao_basis::ordering`] is already sorted (the common case) nothing is
/// read back and nothing is written, so a well-formed calculation is left
/// byte for byte as it was — which is what check mode (contract C2) needs.
/// Otherwise the reference prints
/// `Warning: Basis set is not properly ordered. Redordering.` on **stderr**
/// (typo included) and rewrites both `ao_basis` and, when it is present,
/// `mo_basis` (the coefficients have to follow their AO order).
///
/// # Errors
///
/// Any failure of `ao_basis::read` other than "absent" is swallowed by the
/// reference (its `read` returns `Option`), so a record that cannot be read
/// is simply left alone.  Once a record is read, a failure of `ordering`,
/// `reorder`, `write`, or of the `mo_basis` read/reorder/write propagates and
/// aborts the editor.
pub fn reorder_ao_basis(ezfio: &EzFio) -> Result<()> {
    let Some(aos) = ao_basis::read(ezfio)? else {
        return Ok(());
    };
    let ordering = ao_basis::ordering(&aos)?;
    let mut sorted = ordering.clone();
    sorted.sort_unstable();
    if sorted == ordering {
        return Ok(());
    }
    eprintln!("Warning: Basis set is not properly ordered. Redordering.");
    let new_aos = ao_basis::reorder(&aos)?;
    ao_basis::write(ezfio, &new_aos)?;
    if let Some(mos) = mo_basis::read(ezfio)? {
        let new_mos = mo_basis::reorder(&mos, &ordering)?;
        mo_basis::write(ezfio, &new_mos)?;
    }
    Ok(())
}

/// Reorder the auxiliary AO basis into canonical order.
///
/// Port of the second block of `qp_edit.ml`'s `run`: same algorithm and same
/// warning as [`reorder_ao_basis`], without the MO coefficients (the
/// auxiliary basis has none).
///
/// # Errors
///
/// Same contract as [`reorder_ao_basis`]: an unreadable record is skipped,
/// anything that fails after a successful read propagates.
pub fn reorder_ao_extra_basis(ezfio: &EzFio) -> Result<()> {
    let Some(aos) = ao_extra_basis::read(ezfio)? else {
        return Ok(());
    };
    let ordering = ao_extra_basis::ordering(&aos)?;
    let mut sorted = ordering.clone();
    sorted.sort_unstable();
    if sorted == ordering {
        return Ok(());
    }
    eprintln!("Warning: Basis set is not properly ordered. Redordering.");
    let new_aos = ao_extra_basis::reorder(&aos)?;
    ao_extra_basis::write(ezfio, &new_aos)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_are_unique_and_round_trip() {
        for handler in Handler::ALL {
            assert_eq!(Handler::from_name(handler.name()), Some(*handler));
        }
        assert_eq!(Handler::ALL.len(), 6);
        assert_eq!(Handler::from_name("no_such_handler"), None);
    }

    #[test]
    fn only_nuclei_applies_in_check_mode() {
        let applying: Vec<&str> = Handler::ALL
            .iter()
            .filter(|h| h.applies_in_check())
            .map(|h| h.name())
            .collect();
        assert_eq!(applying, ["nuclei_by_hand"]);
    }
}
