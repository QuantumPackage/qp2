//! Round-trip, interop, and runtime-dimension tests for `qp-ezfio`.
//!
//! These tests do not compare against whole reference files (that is
//! `tests/golden.rs`); they check the properties the crate promises:
//!
//! * every element type survives write -> read, including rank 1/2/3,
//!   single-element and zero-extent arrays;
//! * a Fortran-produced payload, read by the crate and rewritten, reproduces
//!   the OCaml producer's fixture (the crate has one writer, deviation E5),
//!   except for the digits the lossy `real` format cannot carry across (see
//!   `fortran_payload_rewritten_as_ocaml_is_exact_except_lost_real_digits`);
//! * the registry's dimension expressions are resolved at call time from the
//!   scalars currently in the tree, not at start-up.

mod common;

use common::*;
use qp_ezfio::{Array, Error};

/// Write `value`, read it back, and return the reread array.
fn roundtrip<T: qp_ezfio::Element>(group: &str, key: &str, value: &Array<T>) -> Array<T> {
    let (_dir, ezfio) = tree();
    ezfio.set_array(group, key, value).unwrap();
    ezfio.get_array::<T>(group, key).unwrap()
}

#[test]
fn roundtrip_every_type() {
    // Real, rank 1 and rank 2, including both signed zeros.
    let reals = array(&[2, 3], &[1.0f64, -0.5, 1.0e-3, 0.0, -0.0, 1.0e300]);
    let got = roundtrip("g", "reals", &reals);
    assert_eq!(got.dims(), &[2, 3]);
    assert_eq!(got.data(), reals.data());
    assert!(got.data()[4].is_sign_negative());

    // Integer, rank 3 with a leading extent of 1.
    let ints = array(&[1, 2, 3], &[1i64, -42, 1000000, 0, i64::MAX, i64::MIN]);
    let got = roundtrip("g", "ints", &ints);
    assert_eq!(got.dims(), &[1, 2, 3]);
    assert_eq!(got.data(), ints.data());

    // Logical.
    let logicals = array(&[3], &[true, false, true]);
    let got = roundtrip("g", "logicals", &logicals);
    assert_eq!(got.data(), logicals.data());

    // Character.
    let chars: Array<String> = array(
        &[3],
        &["abc".to_string(), "spaced".to_string(), "x".to_string()],
    );
    let got = roundtrip("g", "chars", &chars);
    assert_eq!(got.data(), chars.data());

    // Single element and zero extent.
    let one = array(&[1], &[42.0f64]);
    assert_eq!(roundtrip("g", "one", &one).data(), one.data());
    let empty: Array<f64> = array(&[3, 0], &[]);
    let got = roundtrip("g", "empty", &empty);
    assert_eq!(got.dims(), &[3, 0]);
    assert!(got.is_empty());
}

/// Character payloads are written verbatim and read trimmed.
///
/// `ezfio_core`'s writer keeps the blanks of `"  spaced"` (leading blanks are
/// data, reference OCaml semantics), but its reader applies `String.trim` on
/// both ends. The previous local reader kept the leading blanks; this is a
/// documented read delta.
#[test]
fn character_leading_blanks_are_written_but_trimmed_on_read() {
    let (_dir, ezfio) = tree();
    let chars: Array<String> = array(&[2], &["  spaced".to_string(), "trail  ".to_string()]);
    ezfio.set_array("g", "chars", &chars).unwrap();
    assert_eq!(
        String::from_utf8(written_array(&ezfio, "g", "chars")).unwrap(),
        "  1\n                   2 \n  spaced\ntrail  \n"
    );
    assert_eq!(
        ezfio.get_array::<String>("g", "chars").unwrap().data(),
        &["spaced".to_string(), "trail".to_string()]
    );
}

#[test]
fn roundtrip_scalars_every_type() {
    let (_dir, ezfio) = tree();
    ezfio.set("g", "r", &2.5f64).unwrap();
    ezfio.set("g", "i", &-7i64).unwrap();
    ezfio.set("g", "l", &true).unwrap();
    ezfio.set("g", "c", &"hello".to_string()).unwrap();
    assert_eq!(ezfio.get::<f64>("g", "r").unwrap(), 2.5);
    assert_eq!(ezfio.get::<i64>("g", "i").unwrap(), -7);
    assert!(ezfio.get::<bool>("g", "l").unwrap());
    assert_eq!(ezfio.get::<String>("g", "c").unwrap(), "hello");
}

/// Read the Fortran fixture, rewrite it with the crate's writer (the OCaml
/// dialect of `ezfio_core`, deviation E5), and compare with the OCaml fixture.
///
/// This cannot be byte-identical for `real` data. The Fortran `E24.15` payload
/// carries 15 significant digits and the OCaml `%24.15e` writer emits 16: the
/// digit the Fortran client rounded away is gone, so a value written from a
/// 16-digit literal (`1234567.8901234567`, `pi`) cannot come back. Every other
/// line must still match the OCaml fixture exactly, and exactly one line per
/// payload is allowed to differ.
///
/// (The opposite direction, OCaml bytes rewritten as Fortran, has no test by
/// design: the Fortran-layout writer was removed with `ByteLayout`.)
#[test]
fn fortran_payload_rewritten_as_ocaml_is_exact_except_lost_real_digits() {
    let (_dir, base) = tree();

    // Values that both layouts hold exactly: byte-identical conversion.
    for (key, dims) in [("r2d", vec![2, 3]), ("empty", vec![3, 0])] {
        let payload = fixture(Producer::Fortran, "oracle", key);
        install_array(&base, "oracle", key, &payload);
        let read = base.get_array::<f64>("oracle", key).unwrap();
        assert_eq!(read.dims(), dims.as_slice(), "{key}");
        base.set_array("oracle", key, &read).unwrap();
        assert_eq!(
            written_array(&base, "oracle", key),
            fixture(Producer::Ocaml, "oracle", key),
            "{key}: Fortran -> OCaml conversion is not byte-identical"
        );
    }

    // `lost` is the index of the single element whose 16-digit literal was
    // rounded to 15 digits by the Fortran writer.
    for (key, lost) in [("reals", 3usize), ("redge", 3usize)] {
        let payload = fixture(Producer::Fortran, "oracle", key);
        install_array(&base, "oracle", key, &payload);
        let read = base.get_array::<f64>("oracle", key).unwrap();
        base.set_array("oracle", key, &read).unwrap();

        let written = written_array(&base, "oracle", key);
        let reference = fixture(Producer::Ocaml, "oracle", key);
        let got: Vec<&[u8]> = written.split(|&b| b == b'\n').collect();
        let want: Vec<&[u8]> = reference.split(|&b| b == b'\n').collect();
        assert_eq!(got.len(), want.len(), "{key}: line count");
        for (line, (a, b)) in got.iter().zip(&want).enumerate() {
            let element = line.checked_sub(2);
            if element == Some(lost) {
                assert_ne!(a, b, "{key}: element {lost} should have lost digits");
            } else {
                assert_eq!(a, b, "{key}: line {line}");
            }
        }
    }

    // Integers are written identically by both clients, so the conversion is
    // trivially exact; assert it anyway because it pins the dims-line change.
    let payload = fixture(Producer::Fortran, "oracle", "ints");
    install_array(&base, "oracle", "ints", &payload);
    let read = base.get_array::<i64>("oracle", "ints").unwrap();
    base.set_array("oracle", "ints", &read).unwrap();
    assert_eq!(
        written_array(&base, "oracle", "ints"),
        fixture(Producer::Ocaml, "oracle", "ints")
    );
}

/// `determinants/psi_det` is declared with dims
/// `(n_int*bit_kind/8, 2, n_det)`: the shape is a function of scalars that are
/// written after the array in a real run. The write path resolves them at call
/// time to check the array it is handed; the read path never consults them
/// (the reference getters do not either), so a stale tree stays readable.
#[test]
fn dimension_expressions_resolve_at_call_time() {
    let (_dir, ezfio) = tree();
    ezfio.set("determinants", "n_int", &1i64).unwrap();
    ezfio.set("determinants", "bit_kind", &64i64).unwrap();
    ezfio.set("determinants", "n_det", &2i64).unwrap();

    // n_int * bit_kind / 8 = 8.
    let det: Array<i64> = array(&[8, 2, 2], &(0..32).map(i64::from).collect::<Vec<_>>());
    ezfio.set_array("determinants", "psi_det", &det).unwrap();
    assert_eq!(
        ezfio
            .get_array::<i64>("determinants", "psi_det")
            .unwrap()
            .dims(),
        &[8, 2, 2]
    );

    // Wrong first extent against the current n_int.
    let wrong: Array<i64> = array(&[7, 2, 2], &vec![0i64; 28]);
    let err = ezfio
        .set_array("determinants", "psi_det", &wrong)
        .unwrap_err();
    assert!(matches!(err, Error::Shape { .. }), "{err}");

    // Growing n_int re-interprets the same file: the stored shape no longer
    // matches the declared dims.  The *read* path must stay lenient (the
    // reference getters never consult `size:`), so the array comes back with
    // its on-disk dims and data; only the write path rejects the mismatch.
    ezfio.set("determinants", "n_int", &2i64).unwrap();
    let stale = ezfio
        .get_array::<i64>("determinants", "psi_det")
        .expect("a stale declared shape must not fail a read");
    assert_eq!(stale.dims(), &[8, 2, 2]);
    assert_eq!(stale.data(), det.data());

    // 2*64/8 = 16, so the same values can be rewritten at the new shape.
    let det: Array<i64> = array(&[16, 2, 2], &(0..64).map(i64::from).collect::<Vec<_>>());
    ezfio.set_array("determinants", "psi_det", &det).unwrap();
    assert_eq!(
        ezfio
            .get_array::<i64>("determinants", "psi_det")
            .unwrap()
            .dims(),
        &[16, 2, 2]
    );
}

/// A key the registry does not know (a plugin attribute) is read and written
/// without type or shape checks, like both reference clients.
#[test]
fn unknown_keys_are_unchecked() {
    let (_dir, ezfio) = tree();
    let data: Array<i64> = array(&[7, 3], &(0..21).map(i64::from).collect::<Vec<_>>());
    ezfio.set_array("plugin", "cube", &data).unwrap();
    assert_eq!(
        ezfio.get_array::<i64>("plugin", "cube").unwrap().dims(),
        &[7, 3]
    );

    // The registry would reject this type for a declared key; for an unknown
    // one there is nothing to contradict.
    ezfio.set("plugin", "scale", &1.5f64).unwrap();
    assert_eq!(ezfio.get::<f64>("plugin", "scale").unwrap(), 1.5);
}

/// `EzFio::create` writes only the root marker; a group appears when a key is
/// written and disappears from `groups()` when its marker is removed.
#[test]
fn create_then_list_groups_and_keys() {
    let (_dir, ezfio) = tree();
    assert_eq!(ezfio.groups().unwrap(), vec!["nuclei".to_string()]);
    ezfio.set("zzz", "k", &1i64).unwrap();
    assert_eq!(
        ezfio.groups().unwrap(),
        vec!["nuclei".to_string(), "zzz".to_string()]
    );
    assert_eq!(ezfio.keys("zzz").unwrap(), vec!["k".to_string()]);
    ezfio.remove("zzz", "k").unwrap();
    assert!(ezfio.keys("zzz").unwrap().is_empty());
}
