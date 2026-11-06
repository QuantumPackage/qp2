//! GAMESS input-deck writers (`ocaml/Gamess.ml`).
//!
//! **This module is dead code in QP2**: `ocaml/Gamess.ml` refers to a
//! `Sym.t` module that does not exist anywhere in the tree, so it cannot
//! compile and has no callers. It is ported for completeness, as a set of
//! pure string formatters.
//!
//! The reference file-reading helpers (`read_mos`, `read_until_found`,
//! `read_natural_mos`, `read_canonical_mos`, `string_of_vec`) are **not**
//! ported: they use `Str` regular expressions and file I/O, and this crate is
//! I/O-free. [`Guess::Canonical`] and [`Guess::Natural`] therefore carry the
//! already-extracted `$VEC` block text instead of a file name; [`MoVec`] is
//! just that text.
//!
//! # Units
//!
//! Coordinates are Angstrom in a GAMESS deck (`UNITS=ANGS`), even though the
//! rest of QP2 works in Bohr; the reference writes `Element.to_charge` (an
//! integer number of electrons) as the nuclear charge and does **not**
//! convert the coordinates it receives. Callers are responsible for the
//! unit of the `f64` coordinate triplets.
//!
//! # Preserved OCaml quirks and deviations
//!
//! * `string_of_basis` is defined in the reference as a *partial
//!   application* (`Printf.sprintf " $BASIS\n  GBASIS=%s\n $END"` with no
//!   argument); here it is a plain function of the basis name.
//! * `string_of_drt` takes a symmetry argument and ignores it. Reproduced.
//! * `make_data`'s constructors set `Data::sym` to `"D4h\n"`, `"C4v\n"` or
//!   `"C1\n"`. The reference calls the missing `Sym.to_data`, so the exact
//!   symmetry text is unknowable from this tree; the trailing newline is
//!   required for `string_of_data` to separate the symmetry line from the
//!   coordinates. Pass a different `sym` for anything else.
//! * `create_cas_input` catches `Not_found` around `string_of_vec` and falls
//!   back to a canonical vector; with the file I/O removed nothing can raise,
//!   so only `Natural` is used.
//! * Integer division in `make_drt` and `n_elec_alpha_beta` truncates toward
//!   zero, as OCaml's `int` division does.

use std::fmt;

use crate::gaussian_primitive::ocaml_fmt::format_f;
use crate::pseudo::ElementSymbol;

/// SCF type, the reference `scftyp_t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scftyp {
    /// Restricted Hartree-Fock.
    Rhf,
    /// Restricted open-shell Hartree-Fock.
    Rohf,
    /// Multi-configuration SCF.
    Mcscf,
    /// No SCF (the reference `NONE`).
    None,
}

/// `string_of_scftyp`.
#[must_use]
pub fn string_of_scftyp(scftyp: Scftyp) -> &'static str {
    match scftyp {
        Scftyp::Rhf => "RHF",
        Scftyp::Rohf => "ROHF",
        Scftyp::Mcscf => "MCSCF",
        Scftyp::None => "NONE",
    }
}

/// `$CONTRL` block, the reference `contrl`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Contrl {
    /// SCF type.
    pub scftyp: Scftyp,
    /// Maximum number of SCF iterations.
    pub maxit: i64,
    /// `ISPHER` flag.
    pub ispher: i64,
    /// Total charge (the reference field is called `icharg`).
    pub icharg: i64,
    /// Spin multiplicity.
    pub mult: i64,
    /// MP2 level (`0` or `2`).
    pub mplevl: i64,
}

/// `string_of_contrl`: render the complete `$CONTRL` block.
#[must_use]
pub fn string_of_contrl(contrl: &Contrl) -> String {
    format!(
        " $CONTRL
   EXETYP=RUN COORD=UNIQUE  UNITS=ANGS
   RUNTYP=ENERGY SCFTYP={} CITYP=NONE
   MAXIT={}
   ISPHER={}
   MULT={}
   ICHARG={}
   MPLEVL={}
 $END",
        string_of_scftyp(contrl.scftyp),
        contrl.maxit,
        contrl.ispher,
        contrl.mult,
        contrl.icharg,
        contrl.mplevl
    )
}

/// `make_contrl`: the reference defaults `maxit = 100`, `ispher = 1`,
/// `mplevl = 0`; the public fields can be changed afterwards.
#[must_use]
pub fn make_contrl(scftyp: Scftyp, mult: i64, charge: i64) -> Contrl {
    Contrl {
        scftyp,
        maxit: 100,
        ispher: 1,
        icharg: charge,
        mult,
        mplevl: 0,
    }
}

/// A `$VEC` block, the reference `vec_t`: the text between the `$VEC` and
/// `$END` lines, read from a GAMESS output by the reference's file helpers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoVec {
    /// Canonical orbitals (`--- OPTIMIZED MCSCF MO-S ---` or the closed/open
    /// shell headers).
    Canonical(String),
    /// Natural orbitals (`--- NATURAL ORBITALS OF MCSCF ---` or
    /// `MP2 NATURAL ORBITALS`).
    Natural(String),
}

impl MoVec {
    /// The stored `$VEC` text (`string_of_vec`, with the file read removed).
    #[must_use]
    pub fn text(&self) -> &str {
        match self {
            Self::Canonical(text) | Self::Natural(text) => text,
        }
    }
}

/// Initial-orbital guess, the reference `guess_t`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Guess {
    /// Hückel guess.
    Huckel,
    /// Core-Hamiltonian guess.
    Hcore,
    /// Read `norb` canonical orbitals from the stored text.
    Canonical(i64, String),
    /// Read `norb` natural orbitals from the stored text.
    Natural(i64, String),
}

/// `guess_of_string`: case-insensitive `"huckel"`/`"hcore"`.
///
/// # Errors
///
/// [`GamessError::BadGuess`] (reference `Invalid_argument "Bad MO guess"`)
/// for any other text.
pub fn guess_of_string(text: &str) -> Result<Guess, GamessError> {
    match text.to_lowercase().as_str() {
        "huckel" => Ok(Guess::Huckel),
        "hcore" => Ok(Guess::Hcore),
        _ => Err(GamessError::BadGuess),
    }
}

/// `string_of_guess`: the `$GUESS` block, followed by two newlines and the
/// `$VEC` text for the read-in guesses (and nothing for the analytic ones).
#[must_use]
pub fn string_of_guess(guess: &Guess) -> String {
    let mut result = String::from(" $GUESS\n  GUESS=");
    match guess {
        Guess::Hcore => result.push_str("HCORE\n"),
        Guess::Huckel => result.push_str("HUCKEL\n"),
        Guess::Canonical(norb, _) | Guess::Natural(norb, _) => {
            result.push_str(&format!("MOREAD\n  NORB={norb}\n"));
        }
    }
    result.push_str(" $END");
    match guess {
        Guess::Hcore | Guess::Huckel => {}
        Guess::Natural(_, text) | Guess::Canonical(_, text) => {
            result.push_str("\n\n");
            result.push_str(text);
        }
    }
    result
}

/// `string_of_basis`: the `$BASIS` block for a basis name.
#[must_use]
pub fn string_of_basis(basis: &str) -> String {
    format!(" $BASIS\n  GBASIS={basis}\n $END")
}

/// Geometry, the reference `coord_t`.
#[derive(Debug, Clone, PartialEq)]
pub enum Coord {
    /// A single atom at the origin.
    Atom(ElementSymbol),
    /// A homonuclear diatomic at `±r/2` on `z`.
    DiatomicHomo(ElementSymbol, f64),
    /// A heteronuclear diatomic.
    Diatomic(ElementSymbol, ElementSymbol, f64),
    /// An explicit list of `(element, x, y, z)`.
    Xyz(Vec<(ElementSymbol, f64, f64, f64)>),
}

/// Contents of a `$DATA` block, the reference `data_t`.
///
/// `xyz` is the pre-rendered coordinate text: the reference stores it as a
/// string and [`string_of_data`] concatenates it unchanged.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Data {
    /// Symmetry line, including its terminating newline (see the module
    /// docs for the missing `Sym.to_data`).
    pub sym: String,
    /// Title line.
    pub title: String,
    /// Coordinate lines, each already formatted.
    pub xyz: String,
    /// Total nuclear charge in electrons.
    pub nucl_charge: i64,
}

/// `data_of_atom`.
#[must_use]
pub fn data_of_atom(element: ElementSymbol) -> Data {
    let symbol = element.symbol();
    let charge = charge_of(element);
    Data {
        sym: "D4h\n".to_owned(),
        title: symbol.to_owned(),
        xyz: format!("{symbol}  {charge}.0  0. 0. 0."),
        nucl_charge: charge,
    }
}

/// `data_of_diatomic_homo`: the two nuclei are at `±r/2` on `z`, but only the
/// negative position is written (GAMESS rebuilds the other one from the
/// symmetry group).
///
/// # Errors
///
/// [`GamessError::NonPositiveDistance`] when `r <= 0` or `r` is NaN (the
/// reference asserts `r > 0.`).
pub fn data_of_diatomic_homo(element: ElementSymbol, r: f64) -> Result<Data, GamessError> {
    if r.is_nan() || r <= 0.0 {
        return Err(GamessError::NonPositiveDistance { r });
    }
    let symbol = element.symbol();
    let charge = charge_of(element);
    Ok(Data {
        sym: "D4h\n".to_owned(),
        title: format!("{symbol}2"),
        xyz: format!("{symbol}  {charge}.0  0. 0. {}", format_f(-r * 0.5, 0)),
        nucl_charge: 2 * charge,
    })
}

/// `data_of_diatomic`: first atom at the origin, second at `(0, 0, r)`.
///
/// # Errors
///
/// [`GamessError::NonPositiveDistance`] when `r <= 0` or `r` is NaN.
pub fn data_of_diatomic(
    element1: ElementSymbol,
    element2: ElementSymbol,
    r: f64,
) -> Result<Data, GamessError> {
    if r.is_nan() || r <= 0.0 {
        return Err(GamessError::NonPositiveDistance { r });
    }
    let atom1 = element1.symbol();
    let atom2 = element2.symbol();
    let charge1 = charge_of(element1);
    let charge2 = charge_of(element2);
    Ok(Data {
        sym: "C4v\n".to_owned(),
        title: format!("{atom1}{atom2}"),
        xyz: format!(
            "{atom1}  {charge1}.0  0. 0. 0.\n{atom2}  {charge2}.0  0. 0. {}",
            format_f(r, 0)
        ),
        nucl_charge: charge1 + charge2,
    })
}

/// `data_of_xyz`: the title is the literal `"..."`, and each coordinate line
/// is `"<symbol> <%f charge>  <x> <y> <z>"` (the charge is printed as a
/// float).
#[must_use]
pub fn data_of_xyz(atoms: &[(ElementSymbol, f64, f64, f64)]) -> Data {
    let xyz = atoms
        .iter()
        .map(|(element, x, y, z)| {
            format!(
                "{} {}  {} {} {}",
                element.symbol(),
                format_f(element.to_charge(), 0),
                format_f(*x, 0),
                format_f(*y, 0),
                format_f(*z, 0)
            )
        })
        .collect::<Vec<_>>()
        .join("\n");
    let nucl_charge = atoms
        .iter()
        .map(|(element, _, _, _)| element.to_charge() as i64)
        .sum();
    Data {
        sym: "C1\n".to_owned(),
        title: "...".to_owned(),
        xyz,
        nucl_charge,
    }
}

/// `make_data`: dispatch over the geometry representation.
///
/// # Errors
///
/// Propagates the distance check of the diatomic constructors.
pub fn make_data(coord: &Coord) -> Result<Data, GamessError> {
    match coord {
        Coord::Atom(element) => Ok(data_of_atom(*element)),
        Coord::DiatomicHomo(element, r) => data_of_diatomic_homo(*element, *r),
        Coord::Diatomic(element1, element2, r) => data_of_diatomic(*element1, *element2, *r),
        Coord::Xyz(atoms) => Ok(data_of_xyz(atoms)),
    }
}

/// `string_of_data`: `" $DATA\n<title>\n<sym><xyz>\n $END"`.
#[must_use]
pub fn string_of_data(data: &Data) -> String {
    format!(" $DATA\n{}\n{}{}\n $END", data.title, data.sym, data.xyz)
}

/// `string_of_gugdm2`: an empty string for one state.
///
/// # Errors
///
/// [`GamessError::ZeroStates`] (reference `Invalid_argument "Nstates must be
/// > 0"`) when `nstate < 1`.
pub fn string_of_gugdm2(nstate: i64) -> Result<String, GamessError> {
    match nstate {
        1 => Ok(String::new()),
        nstate if nstate < 1 => Err(GamessError::ZeroStates),
        nstate => {
            let mut weights = String::new();
            for state in 0..nstate {
                if state > 0 {
                    weights.push(',');
                }
                weights.push_str("1.");
            }
            Ok(format!("\n $GUGDM2\n   WSTATE(1)={weights}\n $END\n"))
        }
    }
}

/// `$GUGDIA` block, the reference `gugdia_t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Gugdia {
    /// Number of states.
    pub nstate: i64,
    /// Maximum number of iterations.
    pub itermx: i64,
}

/// `string_of_gugdia`.
#[must_use]
pub fn string_of_gugdia(gugdia: &Gugdia) -> String {
    format!(
        "
 $GUGDIA
  PRTTOL=0.0001
  NSTATE={}
  ITERMX={}
 $END
",
        gugdia.nstate, gugdia.itermx
    )
}

/// `make_gugdia`: the reference default is `itermx = 500`.
///
/// # Errors
///
/// [`GamessError::NonPositiveNstate`] when `nstate <= 0` and
/// [`GamessError::TooFewIterations`] when `itermx <= 1` (the reference
/// asserts both).
pub fn make_gugdia(nstate: i64) -> Result<Gugdia, GamessError> {
    if nstate <= 0 {
        return Err(GamessError::NonPositiveNstate { nstate });
    }
    Ok(Gugdia {
        nstate,
        itermx: 500,
    })
}

/// MCSCF algorithm, the reference `mcscf_t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mcscf {
    /// `FULLNR`.
    Fullnr,
    /// `SOSCF`.
    Soscf,
    /// `FOCAS`.
    Focas,
}

/// `string_of_mcscf`: the whole `$MCSCF` block for one algorithm.
#[must_use]
pub fn string_of_mcscf(mcscf: Mcscf) -> String {
    let flags = match mcscf {
        Mcscf::Focas => "   FOCAS=.T.    SOSCF=.F.   FULLNR=.F.",
        Mcscf::Soscf => "   FOCAS=.F.    SOSCF=.T.   FULLNR=.F.",
        Mcscf::Fullnr => "   FOCAS=.F.    SOSCF=.F.   FULLNR=.T.",
    };
    format!(
        " $MCSCF\n{flags}
   CISTEP=GUGA EKT=.F. QUAD=.F. JACOBI=.f.
   MAXIT=1000
 $END"
    )
}

/// `$DRT` block, the reference `drt_t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Drt {
    /// Number of MCC (core) orbitals.
    pub nmcc: i64,
    /// Number of doubly occupied orbitals.
    pub ndoc: i64,
    /// Number of active alpha orbitals.
    pub nalp: i64,
    /// Number of active virtual orbitals.
    pub nval: i64,
    /// State symmetry.
    pub istsym: i64,
}

/// `make_drt`: derive the DRT occupations from the electron counts.
///
/// `n_elec_alpha`/`n_elec_beta` are the numbers of alpha/beta electrons,
/// `n_e` the number of correlated electrons and `n_act` the number of active
/// orbitals. `istsym` defaults to 1, as in the reference; change the field
/// afterwards for another value.
#[must_use]
pub fn make_drt(n_elec_alpha: i64, n_elec_beta: i64, n_e: i64, n_act: i64) -> Drt {
    let n_elec_tot = n_elec_alpha + n_elec_beta;
    let nmcc = (n_elec_tot - n_e) / 2;
    let ndoc = n_elec_beta - nmcc;
    let nalp = n_elec_alpha - nmcc - ndoc;
    let nval = n_act - ndoc - nalp;
    Drt {
        nmcc,
        ndoc,
        nalp,
        nval,
        istsym: 1,
    }
}

/// `string_of_drt`: the `$DRT` block. The reference takes a `sym` argument
/// and never uses it; kept for signature parity.
#[must_use]
pub fn string_of_drt(drt: &Drt, _sym: &str) -> String {
    format!(
        " $DRT
  NMCC={}
  NDOC={}
  NALP={}
  NVAL={}
  NEXT=0
  ISTSYM={}
  FORS=.TRUE.
  GROUP=C1
  MXNINT= 600000
  NPRT=2
 $END",
        drt.nmcc, drt.ndoc, drt.nalp, drt.nval, drt.istsym
    )
}

/// `string_of_mp2`: the `$MP2` block (always the same text).
#[must_use]
pub fn string_of_mp2() -> String {
    " $MP2\n   MP2PRP=.TRUE.\n $END".to_owned()
}

/// Requested computation, the reference `computation`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Computation {
    /// Hartree-Fock.
    Hf,
    /// MP2.
    Mp2,
    /// CAS with `(n_e, n_a)`.
    Cas(i64, i64),
}

/// Molecular system, the reference `system`.
#[derive(Debug, Clone, PartialEq)]
pub struct System {
    /// Spin multiplicity.
    pub mult: i64,
    /// Total charge.
    pub charge: i64,
    /// Basis-set name for `$BASIS`.
    pub basis: String,
    /// Geometry.
    pub coord: Coord,
}

/// `n_elec`: `nucl_charge - charge`.
///
/// The nuclear charge is summed directly from the geometry, so this cannot
/// fail; the reference's `make_data` would assert on a non-positive diatomic
/// distance, but no distance is used here.
#[must_use]
pub fn n_elec(system: &System) -> i64 {
    nucl_charge(&system.coord) - system.charge
}

/// `n_elec_alpha_beta`: split `n_elec` electrons into `(alpha, beta)` for
/// multiplicity `mult`, with `alpha = (n + mult - 1) / 2` and
/// `beta = n - alpha` (integer division truncating toward zero).
#[must_use]
pub fn n_elec_alpha_beta(system: &System) -> (i64, i64) {
    let n = n_elec(system);
    let alpha = (n + system.mult - 1) / 2;
    (alpha, n - alpha)
}

/// `create_single_det_input`: the HF/MP2 deck.
///
/// `vecfile` is the extracted `$VEC` text; the empty string selects the
/// analytic guess.
///
/// # Errors
///
/// Propagates [`make_data`] and, for MP2, the state-count checks of the
/// callers (none here).
pub fn create_single_det_input(
    mp2: bool,
    guess: &Guess,
    vecfile: &str,
    system: &System,
) -> Result<String, GamessError> {
    let scftyp = if system.mult == 1 {
        Scftyp::Rhf
    } else {
        Scftyp::Rohf
    };
    let (n_elec_alpha, _) = n_elec_alpha_beta(system);
    let mut contrl = make_contrl(scftyp, system.mult, system.charge);
    contrl.mplevl = if mp2 { 2 } else { 0 };
    let guess = if vecfile.is_empty() {
        string_of_guess(guess)
    } else {
        string_of_guess(&Guess::Canonical(n_elec_alpha, vecfile.to_owned()))
    };
    let parts = [
        string_of_contrl(&contrl),
        guess,
        string_of_basis(&system.basis),
        if mp2 { string_of_mp2() } else { String::new() },
        string_of_data(&make_data(&system.coord)?),
    ];
    Ok(parts.join("\n\n"))
}

/// `create_hf_input`.
///
/// # Errors
///
/// Propagates [`create_single_det_input`].
pub fn create_hf_input(
    guess: &Guess,
    vecfile: &str,
    system: &System,
) -> Result<String, GamessError> {
    create_single_det_input(false, guess, vecfile, system)
}

/// `create_mp2_input`.
///
/// # Errors
///
/// Propagates [`create_single_det_input`].
pub fn create_mp2_input(
    guess: &Guess,
    vecfile: &str,
    system: &System,
) -> Result<String, GamessError> {
    create_single_det_input(true, guess, vecfile, system)
}

/// `create_cas_input`: the MCSCF/CAS deck.
///
/// # Errors
///
/// [`GamessError`] from [`make_data`] and from the state-count checks
/// (`nstate`).
pub fn create_cas_input(
    guess: &Guess,
    vecfile: &str,
    nstate: i64,
    system: &System,
    n_e: i64,
    n_a: i64,
) -> Result<String, GamessError> {
    let (n_elec_alpha, n_elec_beta) = n_elec_alpha_beta(system);
    let drt = make_drt(n_elec_alpha, n_elec_beta, n_e, n_a);
    let data = make_data(&system.coord)?;
    let mut contrl = make_contrl(Scftyp::Mcscf, system.mult, system.charge);
    contrl.mplevl = 0;
    let guess = if vecfile.is_empty() {
        string_of_guess(guess)
    } else {
        let norb = drt.nmcc + drt.ndoc + drt.nval + drt.nalp;
        string_of_guess(&Guess::Natural(norb, vecfile.to_owned()))
    };
    let parts = [
        string_of_contrl(&contrl),
        guess,
        string_of_basis(&system.basis),
        string_of_mcscf(Mcscf::Fullnr),
        string_of_drt(&drt, &data.sym),
        string_of_gugdia(&make_gugdia(nstate)?),
        string_of_gugdm2(nstate)?,
        string_of_data(&data),
    ];
    Ok(parts.join("\n\n"))
}

/// `create_input`: dispatch on the requested computation.
///
/// # Errors
///
/// Propagates the constructor of the selected computation.
pub fn create_input(
    computation: Computation,
    system: &System,
    nstate: i64,
    guess: &Guess,
    vecfile: &str,
) -> Result<String, GamessError> {
    match computation {
        Computation::Hf => create_hf_input(guess, vecfile, system),
        Computation::Mp2 => create_mp2_input(guess, vecfile, system),
        Computation::Cas(n_e, n_a) => create_cas_input(guess, vecfile, nstate, system, n_e, n_a),
    }
}

/// Failure modes of the GAMESS writers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GamessError {
    /// Reference `Invalid_argument "Bad MO guess"`.
    BadGuess,
    /// Reference `assert (r > 0.)` in a diatomic constructor (no message).
    NonPositiveDistance {
        /// The rejected bond length.
        r: f64,
    },
    /// Reference `Invalid_argument "Nstates must be > 0"`.
    ZeroStates,
    /// Reference `assert (nstate > 0)` in `make_gugdia` (no message).
    NonPositiveNstate {
        /// The rejected state count.
        nstate: i64,
    },
    /// Reference `assert (itermx > 1)` in `make_gugdia` (no message).
    TooFewIterations {
        /// The rejected iteration count.
        itermx: i64,
    },
}

impl fmt::Display for GamessError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::BadGuess => f.write_str("Bad MO guess"),
            Self::NonPositiveDistance { r } => {
                write!(f, "Gamess: the bond length must be > 0, got {r}")
            }
            Self::ZeroStates => f.write_str("Nstates must be > 0"),
            Self::NonPositiveNstate { nstate } => {
                write!(f, "Gamess: Nstate must be > 0, got {nstate}")
            }
            Self::TooFewIterations { itermx } => {
                write!(f, "Gamess: ITERMX must be > 1, got {itermx}")
            }
        }
    }
}

impl std::error::Error for GamessError {}

/// `Element.to_charge |> Charge.to_int`: the integer nuclear charge.
fn charge_of(element: ElementSymbol) -> i64 {
    element.to_charge() as i64
}

/// Total nuclear charge of a geometry, summed the way `data_of_*` does it.
fn nucl_charge(coord: &Coord) -> i64 {
    match coord {
        Coord::Atom(element) => charge_of(*element),
        Coord::DiatomicHomo(element, _) => 2 * charge_of(*element),
        Coord::Diatomic(element1, element2, _) => charge_of(*element1) + charge_of(*element2),
        Coord::Xyz(atoms) => atoms
            .iter()
            .map(|(element, _, _, _)| charge_of(*element))
            .sum(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn system() -> System {
        let h = ElementSymbol::of_string("H").unwrap();
        System {
            mult: 1,
            charge: 0,
            basis: "sto-3g".to_owned(),
            coord: Coord::Xyz(vec![(h, 0.0, 0.0, 0.0), (h, 0.0, 0.0, 0.7)]),
        }
    }

    #[test]
    fn contrl_block_matches_the_reference_format() {
        let contrl = make_contrl(Scftyp::Rhf, 1, 0);
        assert_eq!(
            string_of_contrl(&contrl),
            " $CONTRL\n   EXETYP=RUN COORD=UNIQUE  UNITS=ANGS\n   \
             RUNTYP=ENERGY SCFTYP=RHF CITYP=NONE\n   MAXIT=100\n   ISPHER=1\n   \
             MULT=1\n   ICHARG=0\n   MPLEVL=0\n $END"
        );
    }

    #[test]
    fn guesses_and_basis() {
        assert_eq!(guess_of_string("HCore"), Ok(Guess::Hcore));
        assert_eq!(guess_of_string("huckel"), Ok(Guess::Huckel));
        assert_eq!(guess_of_string("other"), Err(GamessError::BadGuess));
        assert_eq!(
            string_of_guess(&Guess::Hcore),
            " $GUESS\n  GUESS=HCORE\n $END"
        );
        assert_eq!(
            string_of_guess(&Guess::Canonical(4, "$VEC".to_owned())),
            " $GUESS\n  GUESS=MOREAD\n  NORB=4\n $END\n\n$VEC"
        );
        assert_eq!(
            string_of_basis("cc-pvdz"),
            " $BASIS\n  GBASIS=cc-pvdz\n $END"
        );
    }

    #[test]
    fn data_blocks() {
        let h = ElementSymbol::of_string("H").unwrap();
        let data = data_of_atom(h);
        assert_eq!(data.title, "H");
        assert_eq!(data.xyz, "H  1.0  0. 0. 0.");
        assert_eq!(data.nucl_charge, 1);
        assert_eq!(
            string_of_data(&data),
            " $DATA\nH\nD4h\nH  1.0  0. 0. 0.\n $END"
        );

        let homo = data_of_diatomic_homo(h, 1.4).unwrap();
        assert_eq!(homo.title, "H2");
        assert_eq!(homo.xyz, "H  1.0  0. 0. -0.700000");
        assert_eq!(homo.nucl_charge, 2);

        let o = ElementSymbol::of_string("O").unwrap();
        let hetero = data_of_diatomic(h, o, 1.8).unwrap();
        assert_eq!(hetero.title, "HO");
        assert_eq!(hetero.xyz, "H  1.0  0. 0. 0.\nO  8.0  0. 0. 1.800000");
        assert_eq!(hetero.nucl_charge, 9);
        assert!(matches!(
            data_of_diatomic(h, o, 0.0),
            Err(GamessError::NonPositiveDistance { .. })
        ));

        let xyz = data_of_xyz(&[(o, 0.0, 0.0, 0.0), (h, 0.0, 0.0, 0.96)]);
        assert_eq!(
            xyz.xyz,
            "O 8.000000  0.000000 0.000000 0.000000\nH 1.000000  0.000000 0.000000 0.960000"
        );
        assert_eq!(xyz.nucl_charge, 9);
    }

    #[test]
    fn state_blocks() {
        assert_eq!(string_of_gugdm2(1), Ok(String::new()));
        assert_eq!(
            string_of_gugdm2(3),
            Ok("\n $GUGDM2\n   WSTATE(1)=1.,1.,1.\n $END\n".to_owned())
        );
        assert_eq!(string_of_gugdm2(0), Err(GamessError::ZeroStates));
        assert_eq!(
            string_of_gugdia(&make_gugdia(3).unwrap()),
            "\n $GUGDIA\n  PRTTOL=0.0001\n  NSTATE=3\n  ITERMX=500\n $END\n"
        );
        assert!(matches!(
            make_gugdia(0),
            Err(GamessError::NonPositiveNstate { .. })
        ));
    }

    #[test]
    fn mcscf_and_drt_blocks() {
        assert_eq!(
            string_of_mcscf(Mcscf::Fullnr),
            " $MCSCF\n   FOCAS=.F.    SOSCF=.F.   FULLNR=.T.\n   \
             CISTEP=GUGA EKT=.F. QUAD=.F. JACOBI=.f.\n   MAXIT=1000\n $END"
        );
        let drt = make_drt(5, 5, 4, 4);
        assert_eq!(
            (drt.nmcc, drt.ndoc, drt.nalp, drt.nval, drt.istsym),
            (3, 2, 0, 2, 1)
        );
        let text = string_of_drt(&drt, "ignored");
        assert!(text.contains("NMCC=3\n  NDOC=2\n  NALP=0\n  NVAL=2"));
        assert!(text.contains("ISTSYM=1"));
    }

    #[test]
    fn electron_counts() {
        let system = system();
        assert_eq!(n_elec(&system), 2);
        assert_eq!(n_elec_alpha_beta(&system), (1, 1));
        let mut triplet = system.clone();
        triplet.mult = 3;
        assert_eq!(n_elec_alpha_beta(&triplet), (2, 0));
    }

    #[test]
    fn complete_decks_are_joined_with_blank_lines() {
        let system = system();
        let deck = create_hf_input(&Guess::Huckel, "", &system).unwrap();
        assert!(deck.starts_with(" $CONTRL\n"));
        assert!(deck.contains("\n\n $GUESS\n"));
        assert!(deck.contains("\n\n $BASIS\n"));
        assert!(!deck.contains("$MP2"));
        assert!(deck.ends_with(" $END"));
        let mp2 = create_mp2_input(&Guess::Hcore, "", &system).unwrap();
        assert!(mp2.contains("$MP2"));
        let cas = create_cas_input(&Guess::Hcore, "", 2, &system, 2, 2).unwrap();
        assert!(cas.contains("$MCSCF"));
        assert!(cas.contains("$GUGDIA"));
        let _ = create_input(Computation::Hf, &system, 1, &Guess::Huckel, "");
        let _ = create_input(Computation::Mp2, &system, 1, &Guess::Huckel, "");
        let _ = create_input(Computation::Cas(2, 2), &system, 2, &Guess::Huckel, "");
    }
}
