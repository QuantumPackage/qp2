//! Golden byte tests against payloads produced by the reference clients.
//!
//! Every case below is checked twice, in both directions:
//!
//! 1. **Read**: the committed decompressed payload is installed verbatim in a
//!    temporary tree (gzipped for arrays), then read through the public API and
//!    compared with hand-written expected values.
//! 2. **Write**: the same values are written through [`EzFio::set_array`] /
//!    [`EzFio::set`] with the producer's [`ByteLayout`], the resulting file is
//!    decompressed, and compared **byte for byte** with the fixture.
//!
//! The expected values are the literals the fixture generators were given
//! (`tests/fixtures/src/gen_fortran.f90`, `gen_ocaml.ml`); they do not depend
//! on the formatting code under test.
//!
//! Cases the crate deliberately does not reproduce are marked
//! `reproduced: false` and pinned by a hand-derived byte string in the test
//! that owns them, so a silent change of the deviation is still a failure.
//! `tests/fixtures/README.md` lists the deviations.

mod common;

use common::*;
use qp_ezfio::{Array, Element};

/// One reference-produced array payload.
struct ArrayCase<T> {
    producer: Producer,
    group: &'static str,
    key: &'static str,
    dims: &'static [usize],
    /// Values handed to the fixture generator, i.e. what `set_array` writes.
    values: Vec<T>,
    /// Values the reference payload decodes to through the reader.
    ///
    /// Differs from [`Self::values`] only where the format is lossy: `real`
    /// fields carry 15 (Fortran) or 16 (OCaml) significant digits.
    decoded: Vec<T>,
    /// Whether `set_array` must reproduce the reference bytes exactly.
    reproduced: bool,
}

/// One reference-produced scalar payload.
struct ScalarCase<T> {
    producer: Producer,
    group: &'static str,
    key: &'static str,
    value: T,
    reproduced: bool,
}

fn array_case<T: Clone>(
    producer: Producer,
    group: &'static str,
    key: &'static str,
    dims: &'static [usize],
    values: Vec<T>,
) -> ArrayCase<T> {
    ArrayCase {
        producer,
        group,
        key,
        dims,
        decoded: values.clone(),
        values,
        reproduced: true,
    }
}

fn scalar_case<T>(
    producer: Producer,
    group: &'static str,
    key: &'static str,
    value: T,
) -> ScalarCase<T> {
    ScalarCase {
        producer,
        group,
        key,
        value,
        reproduced: true,
    }
}

/// Install the reference payload, read it back, and compare dims and values.
fn read_array_case<T: Element>(case: &ArrayCase<T>) -> Array<T> {
    let payload = fixture(case.producer, case.group, case.key);
    let (_dir, ezfio) = tree();
    install_array(&ezfio, case.group, case.key, &payload);
    let got = ezfio
        .get_array::<T>(case.group, case.key)
        .unwrap_or_else(|e| {
            panic!(
                "read {}: {e}",
                fixture_path(case.producer, case.group, case.key).display()
            )
        });
    assert_eq!(got.dims(), case.dims, "{}:{} dims", case.group, case.key);
    assert_eq!(
        got.data(),
        case.decoded.as_slice(),
        "{}:{} values",
        case.group,
        case.key
    );
    got
}

/// Write the case with its producer's layout and return the decompressed bytes.
fn write_array_case<T: Element>(case: &ArrayCase<T>) -> Vec<u8> {
    let (_dir, ezfio) = tree();
    let ezfio = ezfio.with_layout(case.producer.layout());
    let value = array(case.dims, &case.values);
    ezfio
        .set_array(case.group, case.key, &value)
        .unwrap_or_else(|e| panic!("write {}:{}: {e}", case.group, case.key));
    written_array(&ezfio, case.group, case.key)
}

fn check_arrays<T: Element>(cases: &[ArrayCase<T>]) {
    for case in cases {
        read_array_case(case);
        if case.reproduced {
            let written = write_array_case(case);
            let expected = fixture(case.producer, case.group, case.key);
            assert_eq!(
                written, expected,
                "{}:{} ({:?}) bytes differ from the reference payload",
                case.group, case.key, case.producer
            );
        }
    }
}

fn read_scalar_case<T: Element>(case: &ScalarCase<T>) -> T {
    let payload = fixture(case.producer, case.group, case.key);
    let (_dir, ezfio) = tree();
    install_scalar(&ezfio, case.group, case.key, &payload);
    ezfio
        .get::<T>(case.group, case.key)
        .unwrap_or_else(|e| panic!("read {}:{}: {e}", case.group, case.key))
}

fn write_scalar_case<T: Element>(case: &ScalarCase<T>) -> Vec<u8> {
    let (_dir, ezfio) = tree();
    let ezfio = ezfio.with_layout(case.producer.layout());
    ezfio
        .set(case.group, case.key, &case.value)
        .unwrap_or_else(|e| panic!("write {}:{}: {e}", case.group, case.key));
    written_scalar(&ezfio, case.group, case.key)
}

fn check_scalars<T: Element>(cases: &[ScalarCase<T>]) {
    for case in cases {
        let got = read_scalar_case(case);
        assert_eq!(got, case.value, "{}:{} value", case.group, case.key);
        if case.reproduced {
            let written = write_scalar_case(case);
            let expected = fixture(case.producer, case.group, case.key);
            assert_eq!(
                written, expected,
                "{}:{} ({:?}) bytes differ from the reference payload",
                case.group, case.key, case.producer
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Real arrays
// ---------------------------------------------------------------------------

/// `1.0, -0.5, 1e-3, 1234567.8901234567, 0.0, -0.0`.
fn real_values() -> Vec<f64> {
    vec![1.0, -0.5, 1.0e-3, 1234567.8901234567, 0.0, -0.0]
}

/// `1e-100, 1e100, -1e-5, pi, 1e300, -2.5e-308`.
fn real_edge_values() -> Vec<f64> {
    vec![
        1.0e-100,
        1.0e100,
        -1.0e-5,
        std::f64::consts::PI,
        1.0e300,
        -2.5e-308,
    ]
}

fn real_cases(producer: Producer) -> Vec<ArrayCase<f64>> {
    // `real` fields are lossy: `E24.15` keeps 15 significant digits, the OCaml
    // `%24.15e` keeps 16. The long literals below therefore do not round-trip,
    // and the read side is checked against the decoded values, not the inputs.
    let mut decoded = real_values();
    decoded[3] = match producer {
        Producer::Ocaml => 1_234_567.890_123_457,
        Producer::Fortran => 1_234_567.890_123_46,
    };
    let mut decoded_edge = real_edge_values();
    if producer == Producer::Fortran {
        // 15 significant digits: `pi` loses its last two. Written as the
        // decimal the field decodes to, parsed so it is not mistaken for
        // `f64::consts::PI` (which is the 16-digit value, element 3 of the
        // OCaml payload).
        decoded_edge[3] = "3.14159265358979".parse().expect("decimal literal");
    }
    vec![
        ArrayCase {
            decoded,
            ..array_case(producer, "oracle", "reals", &[6], real_values())
        },
        ArrayCase {
            decoded: decoded_edge,
            ..array_case(producer, "oracle", "redge", &[6], real_edge_values())
        },
        array_case(
            producer,
            "oracle",
            "r2d",
            &[2, 3],
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        ),
        array_case(producer, "oracle", "empty", &[3, 0], vec![]),
    ]
}

#[test]
fn fortran_real_arrays_match_reference_bytes() {
    let cases = real_cases(Producer::Fortran);
    check_arrays(&cases);
    // `assert_eq!` on f64 treats -0.0 == 0.0; the reference writes the sign,
    // and so must the reader.
    let (_dir, ezfio) = tree();
    let payload = fixture(Producer::Fortran, "oracle", "reals");
    install_array(&ezfio, "oracle", "reals", &payload);
    let got = ezfio.get_array::<f64>("oracle", "reals").unwrap();
    assert!(got.data()[4].is_sign_positive());
    assert!(got.data()[5].is_sign_negative());
}

#[test]
fn ocaml_real_arrays_match_reference_bytes() {
    let cases = real_cases(Producer::Ocaml);
    check_arrays(&cases);
    let (_dir, ezfio) = tree();
    let payload = fixture(Producer::Ocaml, "oracle", "reals");
    install_array(&ezfio, "oracle", "reals", &payload);
    let got = ezfio.get_array::<f64>("oracle", "reals").unwrap();
    assert!(got.data()[4].is_sign_positive());
    assert!(got.data()[5].is_sign_negative());
}

// ---------------------------------------------------------------------------
// Integer arrays
// ---------------------------------------------------------------------------

/// Both producers write `integer` and `integer*8` with the same `I20` field;
/// the Rust side stores both in `i64`.
fn integer_cases(producer: Producer) -> Vec<ArrayCase<i64>> {
    vec![
        array_case(producer, "oracle", "ints", &[4], vec![1, -42, 1000000, 0]),
        array_case(
            producer,
            "oracle",
            "i3d",
            &[1, 2, 3],
            vec![1, 2, 3, 4, 5, 6],
        ),
        array_case(
            producer,
            "oracle",
            "i8s",
            &[3],
            vec![1, -9007199254740993, 123456789012345],
        ),
    ]
}

#[test]
fn fortran_integer_arrays_match_reference_bytes() {
    check_arrays(&integer_cases(Producer::Fortran));
}

#[test]
fn ocaml_integer_arrays_match_reference_bytes() {
    check_arrays(&integer_cases(Producer::Ocaml));
}

#[test]
fn integer_scalars_match_reference_bytes() {
    check_scalars(&[scalar_case(Producer::Fortran, "oracle", "iscalar", 7i64)]);
    check_scalars(&[scalar_case(Producer::Ocaml, "oracle", "iscalar", 7i64)]);
}

// ---------------------------------------------------------------------------
// Logical arrays
// ---------------------------------------------------------------------------

#[test]
fn fortran_logical_arrays_match_reference_bytes() {
    check_arrays(&[array_case(
        Producer::Fortran,
        "oracle",
        "logicals",
        &[3],
        vec![true, false, true],
    )]);
    check_scalars(&[scalar_case(Producer::Fortran, "oracle", "lscalar", true)]);
}

#[test]
fn ocaml_logical_scalar_matches_reference_bytes() {
    // The reference OCaml client's `print_bool` emits `T\n` and the write
    // format appends a second `\n`: the scalar file is 3 bytes (`T\n\n` /
    // `F\n\n`), and every OCaml-client-written EZFIO in the golden corpus
    // carries that form (e.g. `work/empty`). The reader takes the first line.
    check_scalars(&[scalar_case(Producer::Ocaml, "oracle", "lscalar", true)]);
    // The `false` variant has no oracle fixture; pin the bytes directly.
    let case = scalar_case(Producer::Ocaml, "oracle", "lscalar", false);
    assert_eq!(write_scalar_case(&case), b"F\n\n");
}

// ---------------------------------------------------------------------------
// Character arrays
// ---------------------------------------------------------------------------

/// `"abc"`, `"  spaced"` (leading blanks are data) and `"x"`.
fn char_values() -> Vec<String> {
    vec!["abc".to_string(), "  spaced".to_string(), "x".to_string()]
}

#[test]
fn fortran_character_arrays_match_reference_bytes() {
    // `nuclei/nucl_label` is declared `character*(32)`, so the crate knows the
    // field width and must reproduce the blank padding exactly.
    check_arrays(&[array_case(
        Producer::Fortran,
        "nuclei",
        "nucl_label",
        &[3],
        vec!["H".to_string(), "O".to_string(), "O".to_string()],
    )]);

    // `oracle/chars` was written from a `character*(8)` variable, but the key
    // is unknown to the registry, so the crate writes it raw (field width 0)
    // and cannot reproduce the padding. Read it, and pin our bytes.
    let case = ArrayCase {
        reproduced: false,
        ..array_case(Producer::Fortran, "oracle", "chars", &[3], char_values())
    };
    read_array_case(&case);
    let written = write_array_case(&case);
    let expected = concat!(
        "  1\n",
        "                   3\n",
        "abc\n",
        "  spaced\n",
        "x\n",
    );
    assert_eq!(String::from_utf8(written).unwrap(), expected);
}

#[test]
fn ocaml_character_arrays_match_reference_bytes() {
    check_arrays(&[
        array_case(Producer::Ocaml, "oracle", "chars", &[3], char_values()),
        array_case(
            Producer::Ocaml,
            "nuclei",
            "nucl_label",
            &[3],
            vec!["H".to_string(), "O".to_string(), "O".to_string()],
        ),
    ]);

    // Fortran-written `hello` is a literal, so it is not padded either.
    check_scalars(&[scalar_case(
        Producer::Ocaml,
        "oracle",
        "cscalar",
        "hello".to_string(),
    )]);
}

#[test]
fn fortran_padded_character_scalar_is_read_but_not_reproduced() {
    // `character*(32)` scalar: `short` + 27 blanks. The reference reader trims;
    // so does ours. The crate writes scalars verbatim, so the padding is not
    // reproduced (it is a property of the Fortran variable, not of the key).
    let case = ScalarCase {
        reproduced: false,
        ..scalar_case(
            Producer::Fortran,
            "oracle3",
            "short_padded",
            "short".to_string(),
        )
    };
    assert_eq!(read_scalar_case(&case), "short");
    assert_eq!(write_scalar_case(&case), b"short\n");

    // Unpadded Fortran scalar: byte-identical.
    check_scalars(&[scalar_case(
        Producer::Fortran,
        "oracle",
        "cscalar",
        "hello".to_string(),
    )]);
}

// ---------------------------------------------------------------------------
// Real scalars
// ---------------------------------------------------------------------------

#[test]
fn real_scalars_match_reference_bytes() {
    check_scalars(&[scalar_case(Producer::Fortran, "oracle", "rscalar", 2.5f64)]);
    check_scalars(&[scalar_case(Producer::Ocaml, "oracle", "rscalar", 2.5f64)]);
}

// ---------------------------------------------------------------------------
// The OCaml logical-array bug
// ---------------------------------------------------------------------------

#[test]
fn reference_ocaml_logical_array_is_rejected() {
    // `write_bool_array` emits a blank line after every element, which the
    // reference OCaml reader itself cannot parse past the first element. The
    // crate must not silently accept it: reading the *second* logical fails.
    let payload = fixture_raw("reference-bugs/ocaml/logicals.txt");
    assert_eq!(
        String::from_utf8(payload.clone()).unwrap(),
        "  1\n                   3 \nT\n\nF\n\nT\n\n"
    );
    let (_dir, ezfio) = tree();
    install_array(&ezfio, "oracle", "logicals", &payload);
    let err = ezfio
        .get_array::<bool>("oracle", "logicals")
        .expect_err("a blank line is not a logical field");
    assert!(
        matches!(err, qp_ezfio::Error::Malformed { .. }),
        "unexpected error: {err}"
    );
    assert!(format!("{err}").contains("line 4"), "{err}");
}

#[test]
fn our_ocaml_logical_array_writes_one_field_per_line() {
    // Deliberate deviation: one field per line, readable by both readers.
    let case = ArrayCase {
        reproduced: false,
        ..array_case(
            Producer::Ocaml,
            "oracle",
            "logicals",
            &[3],
            vec![true, false, true],
        )
    };
    let written = write_array_case(&case);
    assert_eq!(
        String::from_utf8(written).unwrap(),
        "  1\n                   3 \nT\nF\nT\n"
    );
    // And it round-trips through our own reader.
    let (_dir, ezfio) = tree();
    let ezfio = ezfio.with_layout(qp_ezfio::ByteLayout::Ocaml);
    ezfio
        .set_array("oracle", "logicals", &array(&[3], &[true, false, true]))
        .unwrap();
    assert_eq!(
        ezfio
            .get_array::<bool>("oracle", "logicals")
            .unwrap()
            .data(),
        &[true, false, true]
    );
}
