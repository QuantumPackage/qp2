//! Shared helpers for the `qp-ezfio` integration tests.
//!
//! The fixtures live in `tests/fixtures/<producer>/<group>/<key>.txt` and hold
//! the **decompressed** payloads produced by the reference clients (decision
//! D6: raw `.gz` bytes depend on the gzip implementation and are not
//! comparable). This module turns a fixture payload
//! into the files an EZFIO tree would contain, and reads back the payload the
//! crate wrote, so that a test can compare bytes and values separately.
//!
//! Nothing here goes through the crate's own writer when *installing* a
//! reference payload: `install_array` writes the gzip file itself, so a bug in
//! `set_array` cannot corrupt the input of a read test.

#![allow(dead_code)] // each integration test binary uses a different subset

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use qp_ezfio::{Array, ByteLayout, EzFio};
use tempfile::TempDir;

/// Root of the committed fixture corpus.
pub const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures");

/// Which reference client produced a payload.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Producer {
    /// `gfortran` + `libezfio_irp.a`, `E24.15` / `I20` / `L1` fields.
    Fortran,
    /// Reference OCaml client, `%24.15e` / `%20d` / raw strings.
    Ocaml,
}

impl Producer {
    /// Directory name under `tests/fixtures/`.
    pub fn dir(self) -> &'static str {
        match self {
            Producer::Fortran => "fortran",
            Producer::Ocaml => "ocaml",
        }
    }

    /// Byte layout the crate must emit to reproduce this producer.
    pub fn layout(self) -> ByteLayout {
        match self {
            Producer::Fortran => ByteLayout::Fortran,
            Producer::Ocaml => ByteLayout::Ocaml,
        }
    }
}

/// Read the decompressed reference payload `fixtures/<producer>/<group>/<key>.txt`.
pub fn fixture(producer: Producer, group: &str, key: &str) -> Vec<u8> {
    let path = Path::new(FIXTURES)
        .join(producer.dir())
        .join(group)
        .join(format!("{key}.txt"));
    fs::read(&path).unwrap_or_else(|e| panic!("cannot read fixture {}: {e}", path.display()))
}

/// Read a file under `tests/fixtures/` by relative path (reference-bug cases).
pub fn fixture_raw(relative: &str) -> Vec<u8> {
    let path = Path::new(FIXTURES).join(relative);
    fs::read(&path).unwrap_or_else(|e| panic!("cannot read fixture {}: {e}", path.display()))
}

/// Gzip `payload`. Only used to install reference bytes in a tree; the crate's
/// own writer is tested through `EzFio::set_array`.
pub fn gzip(payload: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(payload).expect("gzip write");
    encoder.finish().expect("gzip finish")
}

/// Decompress a `.gz` payload read from the tree.
pub fn gunzip(compressed: &[u8]) -> Vec<u8> {
    let mut decoder = GzDecoder::new(compressed);
    let mut out = Vec::new();
    decoder.read_to_end(&mut out).expect("gunzip");
    out
}

/// Create an empty EZFIO tree in a temporary directory.
///
/// `nuclei/nucl_num = 3` is set so that the registry-known fixtures
/// (`nuclei/nucl_label`, dims `(nucl_num)`) resolve and their shape check is
/// actually exercised.
pub fn tree() -> (TempDir, EzFio) {
    let dir = tempfile::tempdir().expect("tempdir");
    let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create ezfio");
    ezfio
        .set("nuclei", "nucl_num", &3i64)
        .expect("set nucl_num");
    (dir, ezfio)
}

/// Install a reference array payload verbatim as `<root>/<group>/<key>.gz`.
pub fn install_array(ezfio: &EzFio, group: &str, key: &str, payload: &[u8]) {
    let dir = ezfio.root().join(group);
    fs::create_dir_all(&dir).expect("create group dir");
    fs::write(dir.join(format!("{key}.gz")), gzip(payload)).expect("write reference array");
}

/// Install a reference scalar payload verbatim as `<root>/<group>/<key>`.
pub fn install_scalar(ezfio: &EzFio, group: &str, key: &str, payload: &[u8]) {
    let dir = ezfio.root().join(group);
    fs::create_dir_all(&dir).expect("create group dir");
    fs::write(dir.join(key), payload).expect("write reference scalar");
}

/// Read back the *decompressed* payload of the array the crate wrote.
pub fn written_array(ezfio: &EzFio, group: &str, key: &str) -> Vec<u8> {
    let path = ezfio.root().join(group).join(format!("{key}.gz"));
    let compressed =
        fs::read(&path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    gunzip(&compressed)
}

/// Read back the scalar file the crate wrote.
pub fn written_scalar(ezfio: &EzFio, group: &str, key: &str) -> Vec<u8> {
    let path = ezfio.root().join(group).join(key);
    fs::read(&path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
}

/// Build an `Array<T>` from extents and flattened data.
pub fn array<T: Clone + PartialEq + std::fmt::Debug>(dims: &[usize], data: &[T]) -> Array<T> {
    Array::new(dims.to_vec(), data.to_vec()).expect("consistent array")
}

/// Path of a fixture inside the corpus (for error messages).
pub fn fixture_path(producer: Producer, group: &str, key: &str) -> PathBuf {
    Path::new(FIXTURES)
        .join(producer.dir())
        .join(group)
        .join(format!("{key}.txt"))
}
