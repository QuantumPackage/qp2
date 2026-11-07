//! `Bitmasks`: the determinant packing parameters derived from `mo_num`.
//!
//! Port of `ocaml/Input_bitmasks.ml`.  The group is `bitmasks` and the record
//! is two integers:
//!
//! | field | EZFIO key | derivation |
//! |---|---|---|
//! | `n_int` | `bitmasks/n_int` | `(mo_num - 1) / 64 + 1`, written when absent |
//! | `bit_kind` | `bitmasks/bit_kind` | 8 (bytes per Fortran integer), written when absent |
//!
//! [`read`] returns `None` unless `mo_basis/mo_num` exists: without MOs there
//! is nothing to pack, and the section does not appear in the editor.
//!
//! # Not a section, and not called
//!
//! `qp_edit` never mentions `Bitmasks`.  `ocaml/qp_edit.ml`'s task list has no
//! such entry, and `Input_determinants_by_hand` carries its own private
//! `read_n_int`/`read_bit_kind` with the same body (writing to
//! `determinants/*` instead of `bitmasks/*`).  The module is reached only
//! through `Input.ml`, for the benefit of other tools.  The port keeps it
//! faithful and unit-tested, and marks it `#[allow(dead_code)]` so that
//! `clippy -D warnings` records the absence of a caller instead of failing on
//! it.
//!
//! Units: `n_int` is a number of 64-bit words per spin string; `bit_kind` is
//! the byte size of the Fortran integer holding one word.
#![allow(dead_code)]

use anyhow::Result;
use qp_ezfio::EzFio;

use super::qptypes;

/// The EZFIO group holding the two integers.
const GROUP: &str = "bitmasks";

/// The `Bitmasks.t` record of `Input_bitmasks.ml`, in declaration order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Bitmasks {
    /// Words per spin string; `N_int_number`, at least 1.
    pub n_int: i64,
    /// Bytes per word; `Bit_kind`, one of 1, 2, 4, 8.
    pub bit_kind: i64,
}

/// `Bitmasks.read_n_int`: materialise `bitmasks/n_int` from `mo_num`, then
/// validate it through `N_int_number.of_int`.
///
/// The materialised value goes through `Bitlist.n_int_of_mo_num`, i.e. through
/// `N_int_number.of_int` with the *current* maximum (`determinants/n_int`, or
/// 128 when absent) — so this can fail before the value is even written.
fn read_n_int(ezfio: &EzFio) -> Result<i64> {
    if !ezfio.has(GROUP, "n_int")? {
        let mo_num = ezfio.get::<i64>("mo_basis", "mo_num")?;
        let n_int = qptypes::n_int_of_mo_num(ezfio, mo_num)?;
        ezfio.set::<i64>(GROUP, "n_int", &n_int)?;
    }
    qptypes::n_int_of_int(ezfio, ezfio.get::<i64>(GROUP, "n_int")?)
}

/// `Bitmasks.read_bit_kind`: materialise 8 (from the build configuration),
/// then validate it through `Bit_kind.of_int`.
fn read_bit_kind(ezfio: &EzFio) -> Result<i64> {
    if !ezfio.has(GROUP, "bit_kind")? {
        ezfio.set::<i64>(GROUP, "bit_kind", &qptypes::BIT_KIND)?;
    }
    qptypes::bit_kind_of_int(ezfio.get::<i64>(GROUP, "bit_kind")?)
}

/// `Bitmasks.read`: `None` without `mo_basis/mo_num`.
///
/// The record fields are evaluated in declaration order (`n_int` first), which
/// is the order in which a fresh tree receives the two files.
pub(crate) fn read(ezfio: &EzFio) -> Result<Option<Bitmasks>> {
    if !ezfio.has("mo_basis", "mo_num")? {
        return Ok(None);
    }
    let n_int = read_n_int(ezfio)?;
    let bit_kind = read_bit_kind(ezfio)?;
    Ok(Some(Bitmasks { n_int, bit_kind }))
}

/// `Bitmasks.to_string`: the derived pretty-printer, no caller in `qp_edit`.
pub(crate) fn to_string(data: &Bitmasks) -> String {
    format!(
        "\nn_int              = {}\nbit_kind           = {}\n",
        data.n_int, data.bit_kind
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch() -> (tempfile::TempDir, EzFio) {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        (dir, ezfio)
    }

    #[test]
    fn read_is_none_without_the_mo_basis() {
        let (_dir, ezfio) = scratch();
        assert_eq!(read(&ezfio).expect("read"), None);
        // Nothing was materialised: the guard comes first.
        assert!(!ezfio.has(GROUP, "n_int").expect("has n_int"));
        assert!(!ezfio.has(GROUP, "bit_kind").expect("has bit_kind"));
    }

    #[test]
    fn read_materialises_both_files_from_mo_num() {
        let (_dir, ezfio) = scratch();
        ezfio.set::<i64>("mo_basis", "mo_num", &55).expect("set");
        let data = read(&ezfio).expect("read").expect("some");
        // (55 - 1) / 64 + 1 = 1, and the build's bit_kind is 8.
        assert_eq!(
            data,
            Bitmasks {
                n_int: 1,
                bit_kind: 8
            }
        );
        assert_eq!(ezfio.get::<i64>(GROUP, "n_int").expect("n_int"), 1);
        assert_eq!(ezfio.get::<i64>(GROUP, "bit_kind").expect("bit_kind"), 8);
    }

    #[test]
    fn read_does_not_overwrite_a_stored_value() {
        let (_dir, ezfio) = scratch();
        ezfio.set::<i64>("mo_basis", "mo_num", &55).expect("set");
        ezfio.set::<i64>(GROUP, "n_int", &2).expect("set n_int");
        ezfio
            .set::<i64>(GROUP, "bit_kind", &4)
            .expect("set bit_kind");
        assert_eq!(
            read(&ezfio).expect("read").expect("some"),
            Bitmasks {
                n_int: 2,
                bit_kind: 4
            }
        );
    }

    #[test]
    fn read_rejects_a_stored_bit_kind_outside_the_four_widths() {
        let (_dir, ezfio) = scratch();
        ezfio.set::<i64>("mo_basis", "mo_num", &55).expect("set");
        ezfio.set::<i64>(GROUP, "bit_kind", &3).expect("set");
        let error = read(&ezfio).expect_err("3 is not a width");
        assert_eq!(error.to_string(), "Bit_kind should be (1|2|4|8).");
    }

    #[test]
    fn to_string_is_the_ocaml_layout() {
        let data = Bitmasks {
            n_int: 3,
            bit_kind: 8,
        };
        assert_eq!(
            to_string(&data),
            "\nn_int              = 3\nbit_kind           = 8\n"
        );
    }
}
