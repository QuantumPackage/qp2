//! MD5 digests of serialized objects (`ocaml/To_md5.ml`, `ocaml/MD5.ml`).
//!
//! The QP2 tool layer identifies a basis set by the MD5 of its S-expression
//! serialization: `To_md5.ml` is
//!
//! ```text
//! sexp_of_t t |> Sexp.to_string |> Digest.string |> Digest.to_hex |> MD5.of_string
//! ```
//!
//! `Digest.to_hex` is OCaml's MD5 of the *bytes* of the string (no
//! normalization, no trailing newline), lowercase hex, 32 characters.
//! The string is UTF-8 here because every atom the printers can produce is
//! ASCII (see [`crate::sexplib_compat`]), so `as_bytes` is the byte string
//! OCaml hashes.
//!
//! [`MD5_OF_NONE`] is the digest of the four bytes `None`, written when a
//! basis section cannot be read (`CONTRACTS.md` §3): it is the value of
//! both `ao_extra_md5` and `mo_basis/ao_md5` in that case.

/// `MD5("None")`, the sentinel digest for "no basis".
///
/// `Digest.string "None"` is `6adf97f83acf6453d4a6a4b1070f3754`; the
/// golden corpus records it in every `ao_extra_md5` field.
pub const MD5_OF_NONE: &str = "6adf97f83acf6453d4a6a4b1070f3754";

/// The MD5 of `s`, as 32 lowercase hex characters.
///
/// This is `Digest.to_hex (Digest.string s)`: the digest of the UTF-8 bytes
/// of `s`.
///
/// Cost: `O(len(s))` time and one 16-byte state; the hex formatting
/// allocates one 32-byte `String`.
pub fn md5_hex(s: &str) -> String {
    format!("{:x}", md5::compute(s.as_bytes()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn none_sentinel_matches_the_recorded_constant() {
        // Independently known: MD5 of the ASCII bytes "None".
        assert_eq!(md5_hex("None"), MD5_OF_NONE);
    }

    #[test]
    fn empty_string_is_the_well_known_digest() {
        assert_eq!(md5_hex(""), "d41d8cd98f00b204e9800998ecf8427e");
    }

    #[test]
    fn abc_is_the_rfc_1321_example() {
        assert_eq!(md5_hex("abc"), "900150983cd24fb0d6963f7d28e17f72");
    }
}
