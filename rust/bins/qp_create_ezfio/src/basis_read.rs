//! Basis-set and pseudopotential input: the `-b`/`-p` spec strings, the
//! library fetch, and the per-nucleus basis read.
//!
//! `ocaml/qp_create_ezfio.ml` builds two hash tables before it touches the
//! output directory:
//!
//! ```text
//! "h.basis|1,H:aux.basis|He:he.basis"     -b spec
//!    split_on_char '|' |> rev_map String.trim |> build_basis
//! ```
//!
//! Each entry is processed **right to left**. An entry without `:` is a
//! *principal* basis: its file is fetched, then every nucleus (ghosts
//! included) gets the key `Element.to_string`, but only if the key is not
//! already there — so the rightmost principal entry wins for the elements it
//! covers. An entry `key:basis` is an *auxiliary* basis: `key` is an element
//! symbol, or `i,Element` with `i` a 1-based nucleus index; a plain
//! `Hashtbl.add` shadows any earlier binding, so the **leftmost** auxiliary
//! entry wins. `fetch_channel` is called even when the key is already bound,
//! so a later name that does not exist fails the run even if nothing would
//! have read it.
//!
//! Fetching tries the name as given (relative to the current directory) and
//! then `$QP_ROOT/data/{basis,pseudo}/<name>`; `Sys.file_exists` is true for
//! a directory, and `open_in`/`File::open` accepts one, so a directory is
//! opened successfully and only fails on the first read. Basis reads run the
//! reference's three-attempt chain (see [`read_basis`]), whose error classes
//! are *not* interchangeable: attempt 2 swallows every exception, attempt 1
//! and 3 let `End_of_file` out to a handler that turns it into
//! `Element <E> not found in basis set.`, and everything else escapes with
//! its OCaml class (`Gto.GTO_Read_Failure` unprefixed, `Failure` prefixed).
//! The classification was measured against the reference build; the
//! `qb_types` errors do not carry the OCaml payloads, so
//! [`classify_read_error`] rebuilds them (`String.sub / Bytes.sub`, `hd`,
//! `int_of_string`, `Failed in of_prim_coef_list`, …).
//!
//! Pseudopotentials are stored as the **text of the library file** bound to
//! each element, because `Pseudo.read_element` takes a channel in OCaml and
//! the file may hold several elements. The text is read once per fetch and
//! kept in memory; a read error (a directory, say) yields the empty library,
//! which is what the reference's `try … with _ -> raise End_of_file` inside
//! `Pseudo.find` produces. `Pseudo::read_element` then runs at write time
//! (successor #6), with `Pseudo::empty` when no spec named the element.
//!
//! # Units
//!
//! The module converts nothing. A basis file's exponents are in bohr^-2 and
//! its contraction coefficients in Hartree; a pseudopotential's exponents
//! likewise, its coefficients in Hartree, and `r_power` is an integer power
//! of `r`. The nuclear indices attached to the shells are 1-based positions
//! in the nucleus list, ghosts included.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use qp_types::basis::{self, Basis, NuclNumber};
use qp_types::gto::GtoError;
use qp_types::{Atom, Element, element};

use crate::error::{AppError, OcamlException};
use crate::geometry::strerror;
use crate::ocaml_num;

/// `Nucl_number.get_max ()` above this count: `Nucl_number.of_int` prints
/// `More than 10_000 nuclei` on stdout (`Qptypes.ml:2` defines `warning` as
/// `print_string`, so there is no newline).
const NUCL_WARNING_THRESHOLD: i64 = 10_000;

/// The `max` `Nucl_number.get_max ()` returns while no EZFIO file is open,
/// i.e. at the `-b` spec-parsing call site.
const NUCL_DEFAULT_MAX: i64 = 10_000;

/// One nucleus's contracted shells, numbered and read from the basis file.
#[derive(Debug, Clone)]
pub struct AtomBasis {
    /// 1-based nucleus index: `atom`'s position in the list the action
    /// passes to [`read_basis`], ghosts included.
    pub nucleus: NuclNumber,
    /// The shells, in file order, each already tagged with [`Self::nucleus`].
    pub basis: Basis,
}

/// The basis library channels and the reference's `Hashtbl` keys.
///
/// Built by [`build_basis_table`], consumed by [`read_basis`]; the split is
/// the reference's own: the table (opening every file, failing on a missing
/// name) is complete before the pseudopotential table is built, and the
/// per-nucleus reads happen after both.
#[derive(Debug)]
pub struct BasisTable {
    /// One reader per fetched file, in fetch order. `Basis.find` seeks to
    /// zero before every lookup, so several keys may share one file.
    files: Vec<BufReader<File>>,
    /// `Element.to_string` for a principal entry, `"<index>,<symbol>"` for
    /// an auxiliary one. The value indexes [`Self::files`].
    by_key: HashMap<String, usize>,
}

/// Why one `Basis.read_element` attempt did not return shells.
enum ReadFailure {
    /// `Not_found`: no channel is bound to the key.
    NotFound,
    /// `End_of_file`: the element is absent from the file, or the file ends
    /// inside a shell.
    EndOfFile,
    /// Anything else, already classified for the terminal.
    Fatal(AppError),
}

impl BasisTable {
    /// `basis_channel key`: look the key up, then rewind and scan the file.
    ///
    /// `nucleus` is the tag `Basis.read_element` attaches to every shell;
    /// `wanted` is the element whose header line ends the scan.
    fn read(
        &mut self,
        key: &str,
        nucleus: NuclNumber,
        wanted: Element,
    ) -> Result<Basis, ReadFailure> {
        let Some(&index) = self.by_key.get(key) else {
            return Err(ReadFailure::NotFound);
        };
        basis::read_element(&mut self.files[index], nucleus, wanted).map_err(classify_read_error)
    }

    /// Store a freshly fetched file and return its index.
    fn push(&mut self, file: File) -> usize {
        self.files.push(BufReader::new(file));
        self.files.len() - 1
    }
}

/// `build_basis`: parse the `-b` string and open every named file.
///
/// * `spec` — the raw `-b` argument. Entries are `split '|'` then
///   `String.trim` (the OCaml five-byte trim) and processed right to left.
/// * `nuclei` — the nucleus list, ghosts included; principal entries key
///   themselves by every element in it.
/// * `root` — `$QP_ROOT`, used for the `data/basis` fallback only.
///
/// # Errors
///
/// * [`AppError::Prefixed`] for `Failure`: `Basis <name> not found`,
///   `Expected format is int,Element:basis`, `int_of_string`, a glyph that
///   is not a nucl number.
/// * [`AppError::Uncaught`] for `Element.ElementError` on an auxiliary key,
///   the `Qptypes.ml:537` assertion for a non-positive nucleus index,
///   `Invalid_argument("Nucl_number: <n>")` above [`NUCL_DEFAULT_MAX`], and
///   `Sys_error` when a name exists but cannot be opened.
///
/// The file is opened at fetch time, unconditionally: a spec whose file is
/// missing fails even when every key it would add is already bound.
pub fn build_basis_table(spec: &str, nuclei: &[Atom], root: &Path) -> Result<BasisTable, AppError> {
    let mut table = BasisTable {
        files: Vec::new(),
        by_key: HashMap::new(),
    };
    for entry in spec_entries(spec) {
        match entry.split_once(':') {
            // Principal basis: fetch, then bind every element still unbound.
            None => {
                let index = table.push(fetch(Library::Basis, entry, root)?);
                for atom in nuclei {
                    let key = element::to_string(atom.element).to_owned();
                    // `if not (Hashtbl.mem …)`: the first entry processed
                    // (the rightmost in the spec) keeps the key.
                    table.by_key.entry(key).or_insert(index);
                }
            }
            // Auxiliary basis: parse the key first (the reference does), so
            // an invalid key fails before the fetch, then shadow.
            Some((key, name)) => {
                let key = auxiliary_key(key)?;
                let index = table.push(fetch(Library::Basis, name, root)?);
                table.by_key.insert(key, index);
            }
        }
    }
    Ok(table)
}

/// Read every nucleus's shells through the three-attempt fallback chain.
///
/// `nuclei` must be the same list [`build_basis_table`] saw, in the same
/// order; the index attached to each shell is the nucleus's 1-based
/// position, ghosts included.
///
/// The chain, in the reference's nesting (the `End_of_file` handler wraps
/// the whole expression, so it also sees attempt 3):
///
/// 1. key `"<i>,<Element>"`; a missing key falls to 2, `End_of_file` fails
///    with `Element <Element> not found in basis set.`, anything else
///    escapes;
/// 2. key `"<Element>"`, with `with _` swallowing **every** error;
/// 3. only for [`Element::X`]: read [`Element::H`] from the *ghost's own*
///    channel (key `"X"`), so a ghost borrows the hydrogen section of the
///    file that a principal spec bound to `X`. `End_of_file` here is
///    `Element X not found in basis set.`, a missing key is the fallback
///    failure.
///
/// After all three, `Failure "Basis not found for atom <i> (<Element>)"`.
///
/// # Errors
///
/// See above; errors from attempts 1 and 3 keep their measured class
/// (`GTO_Read_Failure` and `Sys_error` unprefixed, `Failure` prefixed).
pub fn read_basis(table: &mut BasisTable, nuclei: &[Atom]) -> Result<Vec<AtomBasis>, AppError> {
    // `Nucl_number.get_max ()` now reads the EZFIO file, whose `nucl_num`
    // the writer sets to the nucleus count before this loop.
    let max = i64::try_from(nuclei.len()).unwrap_or(i64::MAX);
    let mut result = Vec::with_capacity(nuclei.len());
    for (offset, atom) in nuclei.iter().enumerate() {
        // `do_work` numbers from 1 and never skips; `+ 1` cannot overflow a
        // `Vec` length that fits in `i64` on any supported target.
        let index = i64::try_from(offset.saturating_add(1)).unwrap_or(i64::MAX);
        let nucleus = nucl_number_of_int(index, max)?;
        let basis = read_one_nucleus(table, nucleus, atom.element)?;
        result.push(AtomBasis { nucleus, basis });
    }
    Ok(result)
}

/// The fallback chain of [`read_basis`] for one nucleus.
fn read_one_nucleus(
    table: &mut BasisTable,
    nucleus: NuclNumber,
    wanted: Element,
) -> Result<Basis, AppError> {
    let symbol = element::to_string(wanted);
    // Attempt 1: the per-nucleus key.
    match table.read(&format!("{nucleus},{symbol}"), nucleus, wanted) {
        Ok(basis) => return Ok(basis),
        Err(ReadFailure::Fatal(error)) => return Err(error),
        Err(ReadFailure::EndOfFile) => return Err(element_not_found(wanted)),
        Err(ReadFailure::NotFound) => {}
    }
    // Attempt 2: the element key. `with _ ->` swallows everything, so all
    // three outcomes fall through to attempt 3.
    match table.read(symbol, nucleus, wanted) {
        Ok(basis) => return Ok(basis),
        Err(ReadFailure::NotFound | ReadFailure::EndOfFile | ReadFailure::Fatal(_)) => {}
    }
    // Attempt 3: a ghost reads hydrogen out of the channel bound to X.
    if wanted == Element::X {
        match table.read(symbol, nucleus, Element::H) {
            Ok(basis) => return Ok(basis),
            Err(ReadFailure::Fatal(error)) => return Err(error),
            Err(ReadFailure::EndOfFile) => return Err(element_not_found(wanted)),
            Err(ReadFailure::NotFound) => {}
        }
    }
    Err(AppError::fail(format!(
        "Basis not found for atom {nucleus} ({symbol})"
    )))
}

/// `Failure ("Element " ^ Element.to_string e ^ " not found in basis set.")`.
fn element_not_found(element: Element) -> AppError {
    AppError::fail(format!(
        "Element {} not found in basis set.",
        element::to_string(element)
    ))
}

/// Parse and canonicalise the `key` of an auxiliary spec entry.
///
/// Mirrors
///
/// ```text
/// try Element.of_string key
/// with Element.ElementError _ ->
///   match String_ext.split ~on:',' key with
///   | i :: k :: [] -> (Nucl_number.of_int @@ int_of_string i, Element.of_string k)
///   | _ -> failwith "Expected format is int,Element:basis"
/// ```
///
/// so the element form is tried first and only an `ElementError` falls back
/// to the `i,Element` form. The returned key is `Element.to_string` or
/// `"<i>,<Element>"`, exactly what `basis_channel` looks up.
fn auxiliary_key(key: &str) -> Result<String, AppError> {
    if let Ok(element) = element::of_string(key) {
        return Ok(element::to_string(element).to_owned());
    }
    let mut parts = key.split(',');
    let (Some(index), Some(symbol), None) = (parts.next(), parts.next(), parts.next()) else {
        return Err(AppError::fail("Expected format is int,Element:basis"));
    };
    let index = ocaml_num::int_of_string(index).map_err(|_| AppError::fail("int_of_string"))?;
    let nucleus = nucl_number_of_int(index, NUCL_DEFAULT_MAX)?;
    let element = element::of_string(symbol)
        .map_err(|error| AppError::uncaught(OcamlException::ElementError(error.to_string())))?;
    Ok(format!("{nucleus},{}", element::to_string(element)))
}

/// `Nucl_number.of_int ?max x` (`Qptypes.ml:535-548`).
///
/// * `x` — the requested index, an OCaml `int` (63-bit), hence `i64`.
/// * `max` — `Nucl_number.get_max ()`: [`NUCL_DEFAULT_MAX`] while no EZFIO
///   file is open, the nucleus count once `nuclei/nucl_num` is written.
///   The reference skips the upper check when `max = 1`.
///
/// The assertion `x >= 1` is `Qptypes.ml:537, characters 6-12`; above
/// [`NUCL_WARNING_THRESHOLD`] the warning goes to **stdout** with no
/// newline (`warning = print_string`), and an index above `max` is an
/// unprefixed `Invalid_argument("Nucl_number: <x>")`.
fn nucl_number_of_int(x: i64, max: i64) -> Result<NuclNumber, AppError> {
    if x < 1 {
        return Err(AppError::uncaught(OcamlException::Assert {
            file: "Qptypes.ml".to_owned(),
            line: 537,
            start: 6,
            end: 12,
        }));
    }
    if x > NUCL_WARNING_THRESHOLD {
        // A closed stdout is the only failure; the reference's `print_string`
        // would raise `Sys_error` there, which no reachable call site hits.
        let _ = io::stdout().write_all(b"More than 10_000 nuclei");
    }
    if max != 1 && x > max {
        return Err(nucl_out_of_range(x));
    }
    // `NuclNumber` is an `i32` while the reference's `int` is 63-bit. Every
    // reachable value is at most `max` or at most 10 000, so the conversion
    // succeeds; keep the failure instead of panicking if that ever changes.
    NuclNumber::try_from(x).map_err(|_| nucl_out_of_range(x))
}

/// `Invalid_argument ("Nucl_number: " ^ string_of_int x)`.
fn nucl_out_of_range(x: i64) -> AppError {
    AppError::uncaught(OcamlException::InvalidArgument(format!("Nucl_number: {x}")))
}

/// Map a `basis::read_element` error to the reference's exception class.
///
/// Measured against `~/qp2/ocaml/qp_create_ezfio.native` with an auxiliary
/// entry so that attempt 1 reports the error:
///
/// | library error | reference |
/// | --- | --- |
/// | element absent / EOF inside a shell | `End_of_file` (handled by the chain) |
/// | malformed primitive line | `Gto.GTO_Read_Failure("<line>")`, unprefixed |
/// | shell line shorter than two bytes | `Invalid_argument("String.sub / Bytes.sub")` |
/// | unknown angular momentum | `Failure("Angmom should be …")`, prefixed |
/// | non-numeric shell size | `Failure("int_of_string")`, prefixed |
/// | empty primitive list | `Failure("hd")`, prefixed |
/// | mixed symmetries | `Failure("Failed in of_prim_coef_list")`, prefixed |
/// | I/O error with an errno | `Sys_error("<strerror>")`, unprefixed |
///
/// The mixed-symmetry arm is unreachable through `Gto::read_one` (one shell
/// line has one angular momentum) and is kept for exhaustiveness. A
/// non-UTF-8 line has no OCaml counterpart — the reference carries bytes —
/// and becomes the prefixed `Display` text of the `io::Error`.
fn classify_read_error(error: anyhow::Error) -> ReadFailure {
    if let Some(error) = error.downcast_ref::<io::Error>() {
        return classify_io_error(error);
    }
    if let Some(error) = error.downcast_ref::<GtoError>() {
        return match error {
            // `Basis.read` never returns it; `find` reports its own EOF.
            GtoError::EndOfBasis => ReadFailure::EndOfFile,
            GtoError::Io(io) => classify_io_error(io),
            GtoError::ReadFailure(line) => ReadFailure::Fatal(AppError::uncaught(
                OcamlException::GtoReadFailure(line.clone()),
            )),
            // `String.sub buffer 0 2` on a one-byte line (`Gto.ml:48`).
            GtoError::ShortLine(_) => ReadFailure::Fatal(AppError::uncaught(
                OcamlException::InvalidArgument("String.sub / Bytes.sub".to_owned()),
            )),
            // `Display` is the OCaml `Failure` payload for this variant.
            GtoError::BadAngmom(_) => ReadFailure::Fatal(AppError::fail(error.to_string())),
            // `int_of_string (String_ext.strip n_str)` (`Gto.ml:51`).
            GtoError::BadShellSize(_) => ReadFailure::Fatal(AppError::fail("int_of_string")),
            // `List.hd []` in `Gto.of_prim_coef_list` (`Gto.ml:18`).
            GtoError::EmptyPrimitiveList => ReadFailure::Fatal(AppError::fail("hd")),
            GtoError::MixedSymmetry => {
                ReadFailure::Fatal(AppError::fail("Failed in of_prim_coef_list"))
            }
            // Only `Angmom.Xyz::of_string` produces it; not reachable here.
            GtoError::BadXyz(_) => ReadFailure::Fatal(AppError::fail(error.to_string())),
        };
    }
    ReadFailure::Fatal(AppError::fail(error.to_string()))
}

/// Classify an `io::Error` raised while reading a basis channel.
///
/// `UnexpectedEof` is `End_of_file` — both `Basis.find` (element absent)
/// and `Gto::read_one` (file ends inside a shell) report it that way. Any
/// other error with an errno is the reference's `Sys_error`; `input_line`
/// carries no path, so the payload is the bare `strerror` text (measured:
/// `Sys_error("Is a directory")`).
fn classify_io_error(error: &io::Error) -> ReadFailure {
    if error.kind() == io::ErrorKind::UnexpectedEof {
        return ReadFailure::EndOfFile;
    }
    match error.raw_os_error() {
        Some(code) => {
            ReadFailure::Fatal(AppError::uncaught(OcamlException::SysError(strerror(code))))
        }
        None => ReadFailure::Fatal(AppError::fail(error.to_string())),
    }
}

/// Which `$QP_ROOT/data` subdirectory a name is resolved against.
#[derive(Clone, Copy)]
enum Library {
    /// `-b`: `data/basis`, and a missing name says `Basis … not found`.
    Basis,
    /// `-p`: `data/pseudo`, and a missing name says `Pseudo file ….`.
    Pseudo,
}

impl Library {
    /// The directory name under `$QP_ROOT/data`.
    fn directory(self) -> &'static str {
        match self {
            Self::Basis => "basis",
            Self::Pseudo => "pseudo",
        }
    }
}

/// `fetch_channel`: the name as given, then `$QP_ROOT/data/<library>/<name>`.
///
/// `Sys.file_exists` is true for a directory, and both `open_in` and
/// `File::open` succeed on one; the failure then surfaces on the first read
/// (a `Sys_error` for a basis, an empty library for a pseudo). Opening is
/// what fails here, and `Sys_error("<path>: <strerror>")` is unprefixed.
fn fetch(library: Library, name: &str, root: &Path) -> Result<File, AppError> {
    let local = Path::new(name);
    let remote = long_path(library, name, root);
    let path = if local.exists() {
        local
    } else if remote.exists() {
        remote.as_path()
    } else {
        return Err(AppError::fail(match library {
            Library::Basis => format!("Basis {name} not found"),
            Library::Pseudo => format!("Pseudo file {name} not found."),
        }));
    };
    File::open(path).map_err(|error| {
        let detail = match error.raw_os_error() {
            Some(code) => format!("{}: {}", path.display(), strerror(code)),
            None => format!("{}: {error}", path.display()),
        };
        AppError::uncaught(OcamlException::SysError(detail))
    })
}

/// `Qpackage.root ^ "/data/" ^ dir ^ "/" ^ name`, byte for byte.
///
/// Built from `OsStr` pieces because `QP_ROOT` may not be UTF-8; the
/// reference concatenates strings, so nothing is canonicalised.
fn long_path(library: Library, name: &str, root: &Path) -> PathBuf {
    let mut path = root.as_os_str().to_owned();
    path.push(format!("/data/{}/{name}", library.directory()));
    PathBuf::from(path)
}

/// The pseudopotential library texts bound to the elements.
///
/// `-p` is optional: [`PseudoTable::text`] returns `None` when no entry
/// named the element, and the writer then uses `Pseudo::empty`.
#[derive(Debug, Default)]
pub struct PseudoTable {
    /// One text per fetched file, in fetch order.
    texts: Vec<String>,
    /// Element → index in [`Self::texts`]; the reference's `Hashtbl` key.
    by_element: HashMap<Element, usize>,
}

impl PseudoTable {
    /// `pseudo_channel element`: the library text bound to `element`.
    ///
    /// The text is the whole file, not just the element's block:
    /// `Pseudo::read_element` scans it. `None` means no `-p` entry covered
    /// the element.
    pub fn text(&self, element: Element) -> Option<&str> {
        self.by_element
            .get(&element)
            .map(|&index| self.texts[index].as_str())
    }

    /// Store a freshly fetched text and return its index.
    fn push(&mut self, text: String) -> usize {
        self.texts.push(text);
        self.texts.len() - 1
    }
}

/// `build_pseudo`: parse the optional `-p` string and read every named file.
///
/// * `spec` — the raw `-p` argument, or `None` when the option was absent.
/// * `nuclei` — the nucleus list, ghosts included; principal entries key
///   themselves by every element in it.
/// * `root` — `$QP_ROOT`, used for the `data/pseudo` fallback only.
///
/// The processing order is [`build_basis_table`]'s (right to left, principal
/// entries use `Hashtbl.add` unconditionally so the leftmost wins for
/// pseudo), with two differences: a principal name is used **raw**, while an
/// auxiliary name is ASCII-lowercased; and the auxiliary key must be an
/// element — the `int,Element` form is not accepted here, so `1,H:…` raises
/// `Element.ElementError("Element 1,h unknown")` before the fetch.
///
/// # Errors
///
/// * [`AppError::Prefixed`] for `Failure("Pseudo file <name> not found.")`
///   (note the trailing dot and the lowercased auxiliary name).
/// * [`AppError::Uncaught`] for the auxiliary key's `Element.ElementError`
///   and for `Sys_error` when a name exists but cannot be opened.
pub fn read_pseudo(
    spec: Option<&str>,
    nuclei: &[Atom],
    root: &Path,
) -> Result<PseudoTable, AppError> {
    let mut table = PseudoTable::default();
    let Some(spec) = spec else {
        return Ok(table);
    };
    for entry in spec_entries(spec) {
        match entry.split_once(':') {
            None => {
                let index = table.push(read_library_text(Library::Pseudo, entry, root)?);
                for atom in nuclei {
                    table.by_element.insert(atom.element, index);
                }
            }
            Some((key, name)) => {
                // `Element.of_string key` is not guarded: `1,H` lands here.
                let element = element::of_string(key).map_err(|error| {
                    AppError::uncaught(OcamlException::ElementError(error.to_string()))
                })?;
                // `String.lowercase_ascii pseudo`, then the fetch.
                let index = table.push(read_library_text(
                    Library::Pseudo,
                    &name.to_ascii_lowercase(),
                    root,
                )?);
                table.by_element.insert(element, index);
            }
        }
    }
    Ok(table)
}

/// Fetch a pseudo library and read its whole text.
///
/// `Pseudo.find` wraps its `input_line` in `try … with _ -> raise
/// End_of_file`, and the block reader swallows read errors too, so an
/// unreadable file (a directory, for instance) is an **empty library**, not
/// a failure: `--pseudo=` with `$QP_ROOT/data/pseudo/` as the name exits 0.
///
/// `Read::read_to_string` leaves the buffer untouched when it fails on a
/// first chunk — a directory read error, or invalid UTF-8 — so those yield
/// `""`. Two deviations: a non-UTF-8 file loses its bytes where OCaml would
/// parse them, and a mid-file I/O error keeps the bytes already read where
/// the reference's `input_line` would raise `End_of_file` at that point.
fn read_library_text(library: Library, name: &str, root: &Path) -> Result<String, AppError> {
    let mut file = fetch(library, name, root)?;
    let mut text = String::new();
    let _ = file.read_to_string(&mut text);
    Ok(text)
}

/// `String_ext.split ~on:'|' |> List.rev_map String.trim`.
///
/// `split_on_char` keeps empty fields (`""` is one entry, `"a|"` is `a` and
/// `""`), and `rev_map` processes them from the right.
fn spec_entries(spec: &str) -> impl DoubleEndedIterator<Item = &str> {
    spec.split('|').map(trim_ocaml).rev()
}

/// `String.trim`: exactly the bytes OCaml strips — space, form feed, `\n`,
/// `\r`, tab.
///
/// `str::trim` is not used: it also strips `\u{b}` and Unicode spaces.
fn trim_ocaml(text: &str) -> &str {
    text.trim_matches(|c| matches!(c, ' ' | '\u{0c}' | '\n' | '\r' | '\t'))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    use qp_types::point3d::Units;
    use tempfile::{TempDir, tempdir};

    /// A one-atom list, read in bohr so no conversion is involved.
    fn atoms(symbols: &[&str]) -> Vec<Atom> {
        symbols
            .iter()
            .map(|symbol| {
                Atom::of_string(Units::Bohr, &format!("{symbol} 0 0 0")).expect("fixture atom")
            })
            .collect()
    }

    /// Write `content` under `name` in a fresh directory.
    fn scratch(name: &str, content: &str) -> (TempDir, PathBuf) {
        let dir = tempdir().expect("tempdir");
        let path = dir.path().join(name);
        fs::write(&path, content).expect("write fixture");
        (dir, path)
    }

    /// A fake `$QP_ROOT` with `data/<library>/<name>` = `content`.
    fn fake_root(library: &str, name: &str, content: &str) -> (TempDir, PathBuf) {
        let dir = tempdir().expect("tempdir");
        let data = dir.path().join("data").join(library);
        fs::create_dir_all(&data).expect("create data dir");
        fs::write(data.join(name), content).expect("write fixture");
        let root = dir.path().to_path_buf();
        (dir, root)
    }

    /// The action's two phases, in the reference's order.
    fn read_all(spec: &str, nuclei: &[Atom], root: &Path) -> Result<Vec<AtomBasis>, AppError> {
        let mut table = build_basis_table(spec, nuclei, root)?;
        read_basis(&mut table, nuclei)
    }

    /// One hydrogen and one helium, no ghosts.
    fn h_he() -> Vec<Atom> {
        atoms(&["H", "He"])
    }

    // A file with a hydrogen S shell and a helium S shell.
    const HHE: &str = "H\nS 1\n1 1.0 1.0\n\nHe\nS 1\n1 2.0 1.0\n\n";
    // Hydrogen only.
    const H_ONLY: &str = "H\nS 1\n1 1.0 1.0\n\n";
    // Helium only.
    const HE_ONLY: &str = "He\nS 1\n1 2.0 1.0\n\n";

    #[test]
    fn principal_basis_is_read_for_every_nucleus() {
        let (_dir, path) = scratch("hhe.basis", HHE);
        let spec = path.display().to_string();
        let result = read_all(&spec, &h_he(), Path::new("")).expect("principal basis");
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].nucleus, 1);
        assert_eq!(result[0].basis.len(), 1);
        assert_eq!(result[0].basis[0].1, 1);
        assert_eq!(result[1].nucleus, 2);
        assert_eq!(result[1].basis[0].1, 2);
        // The exponents differ, so the two reads came from the right block.
        assert_eq!(result[0].basis[0].0.lc[0].0.expo, 1.0);
        assert_eq!(result[1].basis[0].0.lc[0].0.expo, 2.0);
    }

    #[test]
    fn principal_basis_is_found_in_the_database() {
        let (_dir, root) = fake_root("basis", "mylib.basis", H_ONLY);
        let result = read_all("mylib.basis", &atoms(&["H"]), &root).expect("database");
        assert_eq!(result[0].basis.len(), 1);
    }

    #[test]
    fn rightmost_principal_spec_wins() {
        // `a` covers H, `b` does not. Processed right to left, `b` binds H
        // first and `a` skips it, so the run fails: the rightmost wins.
        let (_a, a) = scratch("a.basis", H_ONLY);
        let (_b, b) = scratch("b.basis", HE_ONLY);
        let spec = format!("{}|{}", a.display(), b.display());
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("b wins");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Basis not found for atom 1 (H)\n"
        );
        // The other order succeeds, which proves the file really was `a`.
        let spec = format!("{}|{}", b.display(), a.display());
        read_all(&spec, &atoms(&["H"]), Path::new("")).expect("a wins");
    }

    #[test]
    fn auxiliary_spec_overrides_one_nucleus() {
        // Nucleus 1 takes the auxiliary file's hydrogen (exponent 3), the
        // principal file's (exponent 1) covers the rest.
        let (_dir, aux) = scratch("aux.basis", "H\nS 1\n1 3.0 1.0\n\n");
        let (_root, principal) = scratch("hhe.basis", HHE);
        let spec = format!("{}|1,H:{}", principal.display(), aux.display());
        let result = read_all(&spec, &h_he(), Path::new("")).expect("aux");
        assert_eq!(result[0].basis[0].0.lc[0].0.expo, 3.0);
        assert_eq!(result[1].basis[0].0.lc[0].0.expo, 2.0);
    }

    #[test]
    fn leftmost_auxiliary_spec_wins() {
        let (_a, a) = scratch("a.basis", "H\nS 1\n1 3.0 1.0\n\n");
        let (_b, b) = scratch("b.basis", "H\nS 1\n1 4.0 1.0\n\n");
        let spec = format!("1,H:{}|1,H:{}", a.display(), b.display());
        let result = read_all(&spec, &atoms(&["H"]), Path::new("")).expect("aux");
        assert_eq!(result[0].basis[0].0.lc[0].0.expo, 3.0);
    }

    #[test]
    fn a_malformed_section_in_attempt_two_is_swallowed() {
        // Principal entry: attempt 2 reads it and `with _` hides the
        // GTO_Read_Failure, so the reported failure is the fallback one.
        let (_dir, path) = scratch("bad.basis", "H\nS 1\n1 1.0\n\n");
        let spec = path.display().to_string();
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("swallowed");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Basis not found for atom 1 (H)\n"
        );
    }

    #[test]
    fn a_malformed_section_in_attempt_one_is_unprefixed() {
        // Auxiliary entry: attempt 1 reads it, and the error escapes.
        let (_dir, path) = scratch("bad.basis", "H\nS 1\n1 1.0\n\n");
        let spec = format!("1,H:{}", path.display());
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("gto");
        assert_eq!(
            error.render(),
            "Fatal error: exception Gto.GTO_Read_Failure(\"1 1.0\")\n"
        );
    }

    #[test]
    fn eof_in_attempt_one_is_element_not_found() {
        let (_dir, path) = scratch("eof.basis", "H\nS 2\n1 1.0 1.0\n");
        let spec = format!("1,H:{}", path.display());
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("eof");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Element H not found in basis set.\n"
        );
    }

    #[test]
    fn eof_in_attempt_two_is_basis_not_found() {
        let (_dir, path) = scratch("eof.basis", "H\nS 2\n1 1.0 1.0\n");
        let spec = path.display().to_string();
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("eof");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Basis not found for atom 1 (H)\n"
        );
    }

    #[test]
    fn a_ghost_reads_hydrogen_through_its_own_channel() {
        // A principal basis binds the X key to the same file; attempt 3 then
        // reads the H section with the ghost's nucleus index.
        let (_dir, path) = scratch("h.basis", H_ONLY);
        let spec = path.display().to_string();
        let result = read_all(&spec, &atoms(&["H", "X"]), Path::new("")).expect("ghost");
        assert_eq!(result.len(), 2);
        assert_eq!(result[1].nucleus, 2);
        assert_eq!(result[1].basis.len(), 1);
        assert_eq!(result[1].basis[0].1, 2);
    }

    #[test]
    fn a_ghost_without_hydrogen_is_element_x_not_found() {
        let (_dir, path) = scratch("he.basis", HE_ONLY);
        let spec = path.display().to_string();
        let error = read_all(&spec, &atoms(&["X"]), Path::new("")).expect_err("no H");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Element X not found in basis set.\n"
        );
    }

    #[test]
    fn a_missing_basis_file_is_prefixed() {
        let error = read_all("nosuchbasis", &atoms(&["H"]), Path::new("/nonexistent"))
            .expect_err("missing");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Basis nosuchbasis not found\n"
        );
    }

    #[test]
    fn a_bad_auxiliary_key_shape_is_the_expected_format_failure() {
        let error = build_basis_table("1,H,extra:x.basis", &atoms(&["H"]), Path::new(""))
            .expect_err("shape");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Expected format is int,Element:basis\n"
        );
    }

    #[test]
    fn a_zero_nucleus_index_is_the_qptypes_assertion() {
        let error = build_basis_table("0,H:x.basis", &atoms(&["H"]), Path::new(""))
            .expect_err("zero index");
        assert_eq!(
            error.render(),
            "Fatal error: exception File \"Qptypes.ml\", line 537, characters 6-12: Assertion failed\n"
        );
    }

    #[test]
    fn a_nucleus_index_above_ten_thousand_is_invalid_argument() {
        let error = build_basis_table("10001,H:x.basis", &atoms(&["H"]), Path::new(""))
            .expect_err("big index");
        assert_eq!(
            error.render(),
            "Fatal error: exception Invalid_argument(\"Nucl_number: 10001\")\n"
        );
    }

    #[test]
    fn a_non_numeric_nucleus_index_is_int_of_string() {
        let error = build_basis_table("zz,H:x.basis", &atoms(&["H"]), Path::new(""))
            .expect_err("not an int");
        assert_eq!(error.render(), "qp_create_ezfio: int_of_string\n");
    }

    #[test]
    fn a_bare_unknown_auxiliary_key_is_the_expected_format_failure() {
        // `Element.of_string "Q"` fails and the `int,Element` shape does not
        // match `"Q"`, so the reference's fallback `failwith` fires.
        let error = build_basis_table("Q:x.basis", &atoms(&["H"]), Path::new(""))
            .expect_err("unknown element");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Expected format is int,Element:basis\n"
        );
    }

    #[test]
    fn an_unknown_auxiliary_element_is_unprefixed() {
        // `1,Q`: the shape matches, `int_of_string` succeeds, and the
        // element lookup is *outside* the `try`, so its `ElementError`
        // escapes.
        let error = build_basis_table("1,Q:x.basis", &atoms(&["H"]), Path::new(""))
            .expect_err("unknown element");
        assert_eq!(
            error.render(),
            "Fatal error: exception Element.ElementError(\"Element Q unknown\")\n"
        );
    }

    #[test]
    fn an_empty_basis_name_opens_the_database_directory() {
        // `--basis=` → the name is `""`; the local path does not exist but
        // `$QP_ROOT/data/basis/` does, so the directory is opened and its
        // first read fails with EISDIR, swallowed by attempt 2.
        let dir = tempdir().expect("tempdir");
        let data = dir.path().join("data").join("basis");
        fs::create_dir_all(&data).expect("create data dir");
        let error = read_all("", &atoms(&["H"]), dir.path()).expect_err("directory");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Basis not found for atom 1 (H)\n"
        );
    }

    #[test]
    fn a_directory_in_attempt_one_is_an_unprefixed_sys_error() {
        let dir = tempdir().expect("tempdir");
        let spec = format!("1,H:{}", dir.path().display());
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("directory");
        assert_eq!(
            error.render(),
            "Fatal error: exception Sys_error(\"Is a directory\")\n"
        );
    }

    #[test]
    fn a_short_shell_line_is_the_string_sub_failure() {
        let (_dir, path) = scratch("short.basis", "H\nZ\n\n");
        let spec = format!("1,H:{}", path.display());
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("short line");
        assert_eq!(
            error.render(),
            "Fatal error: exception Invalid_argument(\"String.sub / Bytes.sub\")\n"
        );
    }

    #[test]
    fn a_bad_angular_momentum_is_prefixed() {
        let (_dir, path) = scratch("ang.basis", "H\nX 1\n1 1.0 1.0\n\n");
        let spec = format!("1,H:{}", path.display());
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("angmom");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Angmom should be S|P|D|F|G|H|I|J|K|L, not X.\n"
        );
    }

    #[test]
    fn a_non_numeric_shell_size_is_int_of_string() {
        let (_dir, path) = scratch("size.basis", "H\nS x\n1 1.0 1.0\n\n");
        let spec = format!("1,H:{}", path.display());
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("shell size");
        assert_eq!(error.render(), "qp_create_ezfio: int_of_string\n");
    }

    #[test]
    fn an_empty_primitive_list_is_hd() {
        // `S 0` reads no primitive, so `Gto.of_prim_coef_list` hits
        // `List.hd []`. Measured.
        let (_dir, path) = scratch("empty.basis", "H\nS 0\n\n");
        let spec = format!("1,H:{}", path.display());
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("no primitive");
        assert_eq!(error.render(), "qp_create_ezfio: hd\n");
    }

    #[test]
    fn no_pseudo_spec_gives_an_empty_table() {
        let table = read_pseudo(None, &atoms(&["H"]), Path::new("")).expect("no spec");
        assert!(table.text(Element::H).is_none());
    }

    #[test]
    fn principal_pseudo_covers_every_nucleus() {
        let (_dir, path) = scratch("lib.pseudo", "H GEN 1 0\n0\n\n");
        let spec = path.display().to_string();
        let table = read_pseudo(Some(&spec), &h_he(), Path::new("")).expect("principal");
        assert_eq!(table.text(Element::H), Some("H GEN 1 0\n0\n\n"));
        assert_eq!(table.text(Element::He), Some("H GEN 1 0\n0\n\n"));
        assert!(table.text(Element::Li).is_none());
    }

    #[test]
    fn auxiliary_pseudo_names_are_lowercased() {
        let (_dir, path) = scratch("bfd.pseudo", "H GEN 1 0\n0\n\n");
        // The file on disk is lowercase; the spec is not.
        let upper = path.parent().expect("parent").join("BFD.PSEUDO");
        fs::rename(&path, &upper).expect("rename fixture");
        let spec = format!("H:{}", upper.display());
        let table = read_pseudo(Some(&spec), &atoms(&["H"]), Path::new("")).expect("lowercase");
        assert!(table.text(Element::H).is_some());
    }

    #[test]
    fn an_auxiliary_pseudo_key_must_be_an_element() {
        let error = read_pseudo(Some("1,H:bfd"), &atoms(&["H"]), Path::new(""))
            .expect_err("int,Element rejected");
        assert_eq!(
            error.render(),
            "Fatal error: exception Element.ElementError(\"Element 1,h unknown\")\n"
        );
    }

    #[test]
    fn a_missing_pseudo_file_keeps_the_trailing_dot() {
        let error =
            read_pseudo(Some("nosuchpseudo"), &atoms(&["H"]), Path::new("")).expect_err("missing");
        assert_eq!(
            error.render(),
            "qp_create_ezfio: Pseudo file nosuchpseudo not found.\n"
        );
    }

    #[test]
    fn leftmost_pseudo_spec_wins() {
        let (_a, a) = scratch("a.pseudo", "H GEN 1 0\n0\n\n");
        let (_b, b) = scratch("b.pseudo", "H GEN 2 0\n0\n\n");
        let spec = format!("H:{}|H:{}", a.display(), b.display());
        let table = read_pseudo(Some(&spec), &atoms(&["H"]), Path::new("")).expect("first");
        assert_eq!(table.text(Element::H), Some("H GEN 1 0\n0\n\n"));
    }

    #[test]
    fn a_pseudo_is_found_in_the_database() {
        let (_dir, root) = fake_root("pseudo", "mylib", "H GEN 1 0\n0\n\n");
        let table = read_pseudo(Some("mylib"), &atoms(&["H"]), &root).expect("database");
        assert!(table.text(Element::H).is_some());
    }

    #[test]
    fn a_directory_pseudo_is_an_empty_library() {
        // `--pseudo=` resolves to `$QP_ROOT/data/pseudo/`, a directory: the
        // reference swallows EISDIR and exits 0 with empty potentials.
        let dir = tempdir().expect("tempdir");
        let data = dir.path().join("data").join("pseudo");
        fs::create_dir_all(&data).expect("create data dir");
        let table = read_pseudo(Some(""), &atoms(&["H"]), dir.path()).expect("directory");
        assert_eq!(table.text(Element::H), Some(""));
    }

    #[test]
    fn a_late_missing_basis_name_fails_even_when_unused() {
        // The rightmost entry binds H; the leftmost name still has to exist.
        let (_dir, path) = scratch("h.basis", H_ONLY);
        let spec = format!("nosuch|{}", path.display());
        let error = read_all(&spec, &atoms(&["H"]), Path::new("")).expect_err("unused");
        assert_eq!(error.render(), "qp_create_ezfio: Basis nosuch not found\n");
    }

    #[test]
    fn spec_entries_are_trimmed_and_reversed() {
        let entries: Vec<&str> = spec_entries(" a |\tb\t| ").collect();
        assert_eq!(entries, vec!["", "b", "a"]);
    }

    #[test]
    fn the_ocaml_trim_keeps_a_vertical_tab() {
        // OCaml's `String.trim` does not strip U+000B; `str::trim` would.
        assert_eq!(trim_ocaml("\u{b}x\u{b}"), "\u{b}x\u{b}");
        assert_eq!(trim_ocaml(" \t\r\n\u{c}x "), "x");
    }
}
