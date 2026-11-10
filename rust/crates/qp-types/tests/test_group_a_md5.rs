//! Golden md5 test over `tests/golden/md5/basis-md5.jsonl` (the P0 corpus of
//! `CONTRACTS.md` §3), plus the ports of `ocaml/tests/test_basis.ml` and
//! `ocaml/tests/test_gto.ml`, which read the installed basis library.
//!
//! # What each row means
//!
//! One row is `(basis-library file, element)`, and the object that was hashed
//! is the two-nucleus corpus `qp_create_ezfio` produces: the element's shells
//! read twice, tagged `Nucl_number` 1 then 2, with every exponent and
//! coefficient put through the EZFIO `%24.15e` → `float_of_string` round trip
//! ([`corpus_basis`] does all of that). `ao_md5` is
//! `Basis.to_md5` of that object; `n_ao` is its Cartesian function count;
//! `ao_extra_md5` is `MD5("None")` because the corpus has no extra basis.
//!
//! The corpus is **not** modified by this test and is read read-only.
//!
//! # Locating the basis library
//!
//! `$QP_BASIS_DIR`, else `$QP_ROOT/data/basis`, else `$HOME/qp2/data/basis`.
//! When none exists the tests fail with the list of paths they tried: a
//! missing basis library must not turn into a green run.
//!
//! `basis_file_sha256` and `pseudo` are parsed (so a schema change is noticed)
//! but not used: this crate has no SHA-256 and no pseudo-potential support,
//! and for these rows the pseudo does not enter `ao_md5` (the prototype that
//! produced the corpus reproduces the four ECP rows without it).

use std::collections::BTreeSet;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use qp_types::gto::{Angmom, Gto};
use qp_types::{Basis, Element, Fmt, LongBasis, MD5_OF_NONE, basis, element};

// ---------------------------------------------------------------------
// Locating the corpus and the basis library
// ---------------------------------------------------------------------

/// `rust/tests/golden/md5/basis-md5.jsonl`, from the crate directory.
fn corpus_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/md5/basis-md5.jsonl")
}

/// The installed basis-library directory, or a panic listing what was tried.
fn basis_dir() -> PathBuf {
    let mut tried = Vec::new();
    let mut candidates: Vec<PathBuf> = Vec::new();
    if let Some(d) = std::env::var_os("QP_BASIS_DIR") {
        candidates.push(PathBuf::from(d));
    }
    if let Some(root) = std::env::var_os("QP_ROOT") {
        candidates.push(PathBuf::from(root).join("data").join("basis"));
    }
    if let Some(home) = std::env::var_os("HOME") {
        candidates.push(PathBuf::from(home).join("qp2").join("data").join("basis"));
    }
    for candidate in candidates {
        if candidate.is_dir() {
            return candidate;
        }
        tried.push(candidate);
    }
    panic!(
        "basis library not found; set QP_BASIS_DIR. Tried: {}",
        tried
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );
}

// ---------------------------------------------------------------------
// The corpus: a flat JSON object per line, parsed by hand
// ---------------------------------------------------------------------

/// The subset of JSON the corpus uses: a flat object whose values are
/// strings, integers or `null`.
mod json {
    /// One JSON value.
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum Value {
        /// A double-quoted string.
        Str(String),
        /// An integer (the corpus only has non-negative counts).
        Int(i64),
        /// `null`.
        Null,
    }

    /// A byte cursor over one line.
    pub struct Parser<'a> {
        input: &'a [u8],
        pos: usize,
    }

    impl<'a> Parser<'a> {
        /// Start at the beginning of `s`.
        pub fn new(s: &'a str) -> Parser<'a> {
            Parser {
                input: s.as_bytes(),
                pos: 0,
            }
        }

        fn skip_ws(&mut self) {
            while let Some(&c) = self.input.get(self.pos) {
                if matches!(c, b' ' | b'\t' | b'\r' | b'\n') {
                    self.pos += 1;
                } else {
                    break;
                }
            }
        }

        fn eat(&mut self, c: u8) -> bool {
            self.skip_ws();
            if self.input.get(self.pos) == Some(&c) {
                self.pos += 1;
                true
            } else {
                false
            }
        }

        fn expect(&mut self, c: u8) -> Result<(), String> {
            if self.eat(c) {
                Ok(())
            } else {
                Err(format!("expected {:?} at byte {}", char::from(c), self.pos))
            }
        }

        fn string(&mut self) -> Result<String, String> {
            self.expect(b'"')?;
            let mut out = String::new();
            loop {
                let Some(&c) = self.input.get(self.pos) else {
                    return Err("unterminated string".to_owned());
                };
                self.pos += 1;
                match c {
                    b'"' => return Ok(out),
                    b'\\' => {
                        let Some(&e) = self.input.get(self.pos) else {
                            return Err("unterminated escape".to_owned());
                        };
                        self.pos += 1;
                        out.push(match e {
                            b'"' => '"',
                            b'\\' => '\\',
                            b'/' => '/',
                            b'b' => '\u{8}',
                            b'f' => '\u{c}',
                            b'n' => '\n',
                            b'r' => '\r',
                            b't' => '\t',
                            other => {
                                return Err(format!("unsupported escape \\{}", char::from(other)));
                            }
                        });
                    }
                    0x00..=0x7f => out.push(char::from(c)),
                    _ => return Err("non-ASCII byte inside a corpus string".to_owned()),
                }
            }
        }

        fn integer(&mut self) -> Result<i64, String> {
            self.skip_ws();
            let start = self.pos;
            if self.input.get(self.pos) == Some(&b'-') {
                self.pos += 1;
            }
            while matches!(self.input.get(self.pos), Some(c) if c.is_ascii_digit()) {
                self.pos += 1;
            }
            let text =
                std::str::from_utf8(&self.input[start..self.pos]).map_err(|e| e.to_string())?;
            text.parse()
                .map_err(|_| format!("not an integer: {text:?}"))
        }

        fn literal(&mut self, word: &[u8]) -> Result<(), String> {
            self.skip_ws();
            if self.input[self.pos..].starts_with(word) {
                self.pos += word.len();
                Ok(())
            } else {
                Err(format!("expected {}", String::from_utf8_lossy(word)))
            }
        }

        fn value(&mut self) -> Result<Value, String> {
            self.skip_ws();
            match self.input.get(self.pos) {
                Some(b'"') => Ok(Value::Str(self.string()?)),
                Some(b'n') => {
                    self.literal(b"null")?;
                    Ok(Value::Null)
                }
                Some(_) => Ok(Value::Int(self.integer()?)),
                None => Err("unexpected end of line".to_owned()),
            }
        }

        /// Parse one whole flat object; trailing bytes are an error.
        pub fn object(&mut self) -> Result<Vec<(String, Value)>, String> {
            self.expect(b'{')?;
            let mut fields = Vec::new();
            if self.eat(b'}') {
                return Ok(fields);
            }
            loop {
                let key = self.string()?;
                self.expect(b':')?;
                let value = self.value()?;
                fields.push((key, value));
                if self.eat(b',') {
                    continue;
                }
                self.expect(b'}')?;
                break;
            }
            self.skip_ws();
            if self.pos != self.input.len() {
                return Err(format!("trailing data at byte {}", self.pos));
            }
            Ok(fields)
        }
    }
}

/// One row of the corpus.
#[derive(Debug, Clone, Default)]
struct Row {
    ao_md5: String,
    ao_extra_md5: String,
    basis: String,
    #[allow(dead_code)] // parsed for schema checking; no SHA-256 in this crate
    basis_file_sha256: String,
    element: String,
    n_ao: usize,
    #[allow(dead_code)] // parsed for schema checking; pseudo support is group B
    pseudo: Option<String>,
}

impl Row {
    /// The seven keys of the corpus schema, all required, no duplicates.
    const KEYS: [&'static str; 7] = [
        "ao_extra_md5",
        "ao_md5",
        "basis",
        "basis_file_sha256",
        "element",
        "n_ao",
        "pseudo",
    ];

    fn parse(line: &str) -> Result<Row, String> {
        let fields = json::Parser::new(line).object()?;
        let mut row = Row::default();
        let mut seen = BTreeSet::new();
        for (key, value) in fields {
            if !seen.insert(key.clone()) {
                return Err(format!("duplicate key {key:?}"));
            }
            match (key.as_str(), value) {
                ("ao_md5", json::Value::Str(s)) => row.ao_md5 = s,
                ("ao_extra_md5", json::Value::Str(s)) => row.ao_extra_md5 = s,
                ("basis", json::Value::Str(s)) => row.basis = s,
                ("basis_file_sha256", json::Value::Str(s)) => row.basis_file_sha256 = s,
                ("element", json::Value::Str(s)) => row.element = s,
                ("n_ao", json::Value::Int(n)) => {
                    row.n_ao = usize::try_from(n).map_err(|_| format!("n_ao {n} is negative"))?;
                }
                ("pseudo", json::Value::Null) => row.pseudo = None,
                ("pseudo", json::Value::Str(s)) => row.pseudo = Some(s),
                (k, v) => return Err(format!("unexpected value {v:?} for key {k:?}")),
            }
        }
        for key in Row::KEYS {
            if !seen.contains(key) {
                return Err(format!("missing key {key:?}"));
            }
        }
        if row.basis_file_sha256.len() != 64 {
            return Err("basis_file_sha256 is not 64 hex characters".to_owned());
        }
        Ok(row)
    }
}

/// All rows, in file order. Blank lines are skipped; every other line must
/// parse.
fn load_rows() -> Vec<Row> {
    let path = corpus_path();
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read the corpus {}: {e}", path.display()));
    text.lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty())
        .map(|(i, line)| Row::parse(line).unwrap_or_else(|e| panic!("corpus line {}: {e}", i + 1)))
        .collect()
}

/// `Element.of_string` on the corpus spelling (the long name, uppercase).
fn element_of(row: &Row) -> Element {
    element::of_string(&row.element)
        .unwrap_or_else(|e| panic!("row {}/{}: {e}", row.basis, row.element))
}

/// `float_of_string (Printf.sprintf "%24.15e" x)`: the EZFIO real-array
/// round trip every `ao_expo`/`ao_coef` went through before it was hashed.
fn ezfio_round_trip(x: f64) -> f64 {
    format!("{x:.15e}")
        .parse::<f64>()
        .expect("internal invariant: Rust's `e` format always re-parses")
}

/// The `qp_create_ezfio` sequence the corpus digest is taken over: the
/// element's section read twice (nucleus 1 then nucleus 2, the stream
/// rewound by [`basis::find`] between the reads) with every exponent and
/// coefficient put through the EZFIO `%24.15e` → `float_of_string` round
/// trip.
fn corpus_basis(path: &Path, element: Element) -> anyhow::Result<Basis> {
    let file = File::open(path).map_err(|e| anyhow::anyhow!("opening {}: {e}", path.display()))?;
    let mut r = BufReader::new(file);
    let mut b = basis::read_element(&mut r, 1, element)
        .map_err(|e| anyhow::anyhow!("reading basis section {element}: {e:#}"))?;
    b.extend(
        basis::read_element(&mut r, 2, element)
            .map_err(|e| anyhow::anyhow!("reading basis section {element}: {e:#}"))?,
    );
    for (gto, _) in &mut b {
        for (primitive, coef) in &mut gto.lc {
            primitive.expo = ezfio_round_trip(primitive.expo);
            *coef = ezfio_round_trip(*coef);
        }
    }
    Ok(b)
}

// ---------------------------------------------------------------------
// The golden test
// ---------------------------------------------------------------------

#[test]
fn every_golden_row_reproduces_ao_md5_and_n_ao() {
    let dir = basis_dir();
    let rows = load_rows();
    assert_eq!(rows.len(), 27, "the P0 corpus has 27 rows");

    let mut failures = Vec::new();
    for row in &rows {
        let element = element_of(row);
        let path = dir.join(&row.basis);
        let b = match corpus_basis(&path, element) {
            Ok(b) => b,
            Err(e) => {
                failures.push(format!("{}/{}: {e:#}", row.basis, row.element));
                continue;
            }
        };
        let md5 = basis::to_md5(&b);
        if md5 != row.ao_md5 {
            failures.push(format!(
                "{}/{}: ao_md5 {md5} != {}",
                row.basis, row.element, row.ao_md5
            ));
        }
        let n_ao = LongBasis::of_basis(&b).len();
        if n_ao != row.n_ao {
            failures.push(format!(
                "{}/{}: n_ao {} != {}",
                row.basis, row.element, n_ao, row.n_ao
            ));
        }
        // Corpus property, not a computation of ours: no row carries an
        // extra basis, so the reference wrote MD5("None") for all of them.
        if row.ao_extra_md5 != MD5_OF_NONE {
            failures.push(format!(
                "{}/{}: ao_extra_md5 {} is not MD5(\"None\")",
                row.basis, row.element, row.ao_extra_md5
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} rows failed:\n{}",
        failures.len(),
        rows.len(),
        failures.join("\n")
    );
}

#[test]
fn the_two_documented_collisions_hold_on_the_computed_digests() {
    // `CONTRACTS.md` §3 predicts these by construction: the extra
    // polarization functions change `n_ao` but not the short basis, and the
    // digest is taken over the short basis. Checking it on two independently
    // computed digests does not trust the corpus.
    let dir = basis_dir();
    let rows = load_rows();
    for (element, a, b) in [
        (Element::Kr, "6-31g_star", "6-31g_star_star"),
        (Element::Rn, "def2-svp", "def2-tzvpp"),
    ] {
        let read = |basis: &str| {
            let b = corpus_basis(&dir.join(basis), element)
                .unwrap_or_else(|e| panic!("{basis}/{element}: {e:#}"));
            qp_types::basis::to_md5(&b)
        };
        let ma = read(a);
        let mb = read(b);
        assert_eq!(ma, mb, "{a} and {b} on {element} must collide");
        for (basis, md5) in [(a, &ma), (b, &mb)] {
            let row = rows
                .iter()
                .find(|r| r.basis == basis && element_of(r) == element)
                .unwrap_or_else(|| panic!("row {basis}/{element} is missing"));
            assert_eq!(md5, &row.ao_md5, "{basis}/{element}");
        }
    }
}

// ---------------------------------------------------------------------
// Ports of ocaml/tests/test_basis.ml and ocaml/tests/test_gto.ml
// ---------------------------------------------------------------------

/// `test_basis.ml`: `cc-pvdz` fluorine, nuclei 1 and 2, plus the long basis
/// round trip, the two writers and the digest.
///
/// The reference values (`nbasis = 12`, `nlong = 30`,
/// `md5 = 227bd952313a812f889a6aa83194c161`) come from running the OCaml
/// test against the same basis library.
#[test]
fn cc_pvdz_fluorine_two_nuclei_matches_the_reference_test() {
    let path = basis_dir().join("cc-pvdz");
    let file = File::open(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut r = BufReader::new(file);

    let basis: Basis = basis::read_element(&mut r, 1, Element::F)
        .expect("fluorine, nucleus 1")
        .into_iter()
        .chain(basis::read_element(&mut r, 2, Element::F).expect("fluorine, nucleus 2"))
        .collect();

    // Six shells per nucleus: 2 x (S 9, S 9, S 1, P 4, P 1, D 1).
    assert_eq!(basis.len(), 12);
    assert_eq!(basis.iter().filter(|(_, n)| *n == 1).count(), 6);
    assert_eq!(basis.iter().filter(|(_, n)| *n == 2).count(), 6);

    let long = LongBasis::of_basis(&basis);
    // 2 x (1 + 1 + 1 + 3 + 3 + 6) Cartesian functions.
    assert_eq!(long.len(), 30);

    let short = long.to_basis();
    assert_eq!(short, basis, "(short_basis <> basis)");

    assert_eq!(basis::to_md5(&basis), "227bd952313a812f889a6aa83194c161");

    // `Basis.to_string` (the Gamess layout, no element array): the reference
    // conses the separator before the header, so the blank line sits
    // between the last shell of one block and `Atom 2`.
    let gamess = basis::to_string(Fmt::Gamess, &basis, None);
    assert!(
        gamess.starts_with("Atom 1\nS   9\n  1     14710.000000          0.000721\n"),
        "{gamess}"
    );
    assert!(gamess.contains("\n\nAtom 2\nS   9\n"), "{gamess}");
    assert!(!gamess.ends_with('\n'), "{gamess}");

    // With an element array the header is the symbol, not `Atom n`.
    let with_elements = basis::to_string_gamess(&basis, Some(&[Element::F, Element::F]));
    assert!(with_elements.starts_with("F\nS   9\n"), "{with_elements}");
    assert!(with_elements.contains("\n\nF\nS   9\n"), "{with_elements}");

    // The Gaussian layout separates blocks with `****` and ends with one.
    let gaussian = basis::to_string_gaussian(&basis, Some(&[Element::F, Element::F]));
    assert!(gaussian.starts_with("F\nS   9   1.00\n"), "{gaussian}");
    assert!(gaussian.contains("\n****\nF\n"), "{gaussian}");
    assert!(gaussian.ends_with("\n****"), "{gaussian}");
}

/// `test_gto.ml`: `Gto.read_one` reads the hydrogen section shell by shell,
/// and re-reading from the top gives the same shell.
///
/// The reference test seeks to offset 0, skips one line and re-reads; the
/// portable statement of what it checks is that the reader is deterministic
/// and that the file's first two hydrogen shells are the ones asserted here.
#[test]
fn hydrogen_shells_are_read_deterministically() {
    let path = basis_dir().join("cc-pvdz");

    let read_first_two = |path: &Path| {
        let mut r = BufReader::new(File::open(path).expect("cc-pvdz is installed"));
        basis::find(&mut r, Element::H).expect("hydrogen section");
        let first = Gto::read_one(&mut r).expect("first hydrogen shell");
        let second = Gto::read_one(&mut r).expect("second hydrogen shell");
        (first, second)
    };

    let (first, second) = read_first_two(&path);
    let (again, _) = read_first_two(&path);
    assert_eq!(again, first, "the reader must be deterministic");

    // The first shell is `S   4` with the cc-pVDZ hydrogen exponents.
    assert_eq!(first.sym, Angmom::S);
    assert_eq!(first.lc.len(), 4);
    let expo: Vec<f64> = first.lc.iter().map(|(p, _)| p.expo).collect();
    assert_eq!(expo, [13.01, 1.962, 0.4446, 0.122]);
    let coef: Vec<f64> = first.lc.iter().map(|(_, c)| *c).collect();
    assert_eq!(coef, [0.019685, 0.137977, 0.478148, 0.50124]);
    assert_eq!(first.to_string_gamess().lines().count(), 5);

    // The second is the single-exponent `S   1`, so `gto2 = gto` from the
    // reference test is (correctly) false for today's file.
    assert_ne!(second, first);
    assert_eq!(second.sym, Angmom::S);
    assert_eq!(second.lc.len(), 1);
    assert_eq!(second.lc[0].0.expo, 0.122);

    // `test_gto_2.ml`: `Basis.read` after the header stops at the blank line.
    let mut r = BufReader::new(File::open(&path).expect("cc-pvdz is installed"));
    basis::find(&mut r, Element::H).expect("hydrogen section");
    let section = basis::read(&mut r, 1).expect("hydrogen section");
    assert_eq!(section.len(), 3, "S 4, S 1, P 1");
    assert_eq!(section[2].0.sym, Angmom::P);
    assert_eq!(section[2].0.lc[0].0.expo, 0.727);
    assert!(section.iter().all(|(_, n)| *n == 1));
}
