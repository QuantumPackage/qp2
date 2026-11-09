//! `Electrons`: the electron counts and the spin they imply.
//!
//! Port of `ocaml/Input_electrons.ml`.  The section is two scalars,
//! `elec_alpha_num` and `elec_beta_num`; everything else is derived:
//!
//! * `read` returns `None` when `elec_alpha_num` is absent — the section then
//!   disappears from the editor;
//! * `to_rst` prints the multiplicity implied by the two counts;
//! * `write` writes both counts back, in record order;
//! * `of_rst` uses the generic sexp-of-RST path (`Generic_input_of_rst`),
//!   because the generated record is exactly two integers.
//!
//! The `elec_num` provider in `EZFIO.cfg` is **not** part of the section: the
//! hand-written record does not carry it, and the OCaml `read_elec_num` only
//! adds the two counts for `to_string`, which `qp_edit` never calls.  Nothing
//! here writes `electrons/elec_num`.
//!
//! Units: none; electron counts are integers.

use anyhow::Result;
use qp_ezfio::EzFio;
use qp_types::Multiplicity;

use super::qptypes;
use super::{ParseResult, apply_parsed};
use crate::engine::make_header;
use crate::ocaml_compat::assertion;
use crate::sexp;

/// The EZFIO group holding the counters.
const GROUP: &str = "electrons";

/// The `Electrons.t` record of `Input_electrons.ml`, in declaration order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Electrons {
    /// Number of alpha electrons; `Elec_alpha_number`, strictly positive.
    pub elec_alpha_num: i64,
    /// Number of beta electrons; `Elec_beta_number`, non-negative.
    pub elec_beta_num: i64,
}

/// `Electrons.read`: `None` when `electrons/elec_alpha_num` is absent.
///
/// The positivity checks are the `of_int` conversions of
/// `Elec_alpha_number` / `Elec_beta_number`.  OCaml raises `Assert_failure`
/// there; this port returns the same classification
/// ([`crate::ocaml_compat::AssertFailure`]) so the driver exits 2 like the
/// reference instead of panicking on data.
pub fn read(ezfio: &EzFio) -> Result<Option<Electrons>> {
    if !ezfio.has(GROUP, "elec_alpha_num")? {
        return Ok(None);
    }
    let elec_alpha_num = ezfio.get::<i64>(GROUP, "elec_alpha_num")?;
    let elec_beta_num = ezfio.get::<i64>(GROUP, "elec_beta_num")?;
    qptypes::elec_alpha_of_int(elec_alpha_num)?;
    qptypes::elec_beta_of_int(elec_beta_num)?;
    Ok(Some(Electrons {
        elec_alpha_num,
        elec_beta_num,
    }))
}

/// `Electrons.write`: the two counters, `elec_alpha_num` first.
pub(crate) fn write(ezfio: &EzFio, data: &Electrons) -> Result<()> {
    ezfio.set::<i64>(GROUP, "elec_alpha_num", &data.elec_alpha_num)?;
    ezfio.set::<i64>(GROUP, "elec_beta_num", &data.elec_beta_num)?;
    Ok(())
}

/// `Electrons.to_rst`, byte for byte.
///
/// The multiplicity line is `Multiplicity.of_alpha_beta`, whose `a >= b`
/// assertion (an `Assert_failure` in OCaml) is reported as an error here.
pub(crate) fn to_rst(data: &Electrons) -> Result<String> {
    let alpha = i32::try_from(data.elec_alpha_num).map_err(|_| {
        anyhow::anyhow!(
            "electron count {} does not fit in an int",
            data.elec_alpha_num
        )
    })?;
    let beta = i32::try_from(data.elec_beta_num).map_err(|_| {
        anyhow::anyhow!(
            "electron count {} does not fit in an int",
            data.elec_beta_num
        )
    })?;
    if alpha < beta {
        return Err(assertion(format!(
            "Multiplicity.of_alpha_beta: assertion (a >= b) failed for a = {alpha}, b = {beta}"
        )));
    }
    let multiplicity = Multiplicity::of_alpha_beta(alpha, beta);
    Ok(format!(
        "\nSpin multiplicity is {}.\n\nNumber of alpha and beta electrons ::\n\n  \
         elec_alpha_num = {}\n  elec_beta_num  = {}\n\n",
        multiplicity, data.elec_alpha_num, data.elec_beta_num
    ))
}

/// `get` for the section: header plus body, or `""` when `read` is `None`.
pub(crate) fn get(ezfio: &EzFio, title: &str) -> Result<String> {
    match read(ezfio)? {
        Some(data) => Ok(format!("{}{}", make_header(title), to_rst(&data)?)),
        None => Ok(String::new()),
    }
}

/// `Electrons.of_rst` = `Generic_input_of_rst.of_rst t_of_sexp`.
///
/// The filter, the `(key value)` rebuilding and the wrapping are the ones of
/// `Generic_input_of_rst`: every line containing `=` becomes a pair, so the
/// prose lines of the RST are ignored.  `t_of_sexp` for a record of two
/// `int`-aliased fields accepts any integer; it does not run `of_int`, so a
/// negative count parses and is written by `write`, exactly as in OCaml.
pub(crate) fn of_rst(body: &str) -> ParseResult<Electrons> {
    let mut source = String::from("(");
    for line in body.split('\n').filter(|line| line.contains('=')) {
        source.push('(');
        source.push_str(&line.replace('=', " "));
        source.push(')');
    }
    source.push(')');

    let names = ["elec_alpha_num", "elec_beta_num"];
    let node = match sexp::parse(&source) {
        Ok(node) => node,
        Err(error) => return ParseResult::Reported(error.message),
    };
    let pairs = match sexp::record_fields(&node, &names) {
        Ok(pairs) => pairs,
        Err(error) => return ParseResult::Reported(error.message),
    };

    let mut values = [0_i64; 2];
    for (slot, name) in values.iter_mut().zip(names) {
        // `record_fields` proved completeness.
        let Some((_, _, value_node)) = pairs.iter().find(|(seen, _, _)| *seen == name) else {
            unreachable!("record_fields returned an incomplete record")
        };
        let Some(text) = value_node.atom() else {
            return ParseResult::Reported("int_of_sexp: expected an atom".to_string());
        };
        match text.parse::<i64>() {
            Ok(value) => *slot = value,
            Err(_) => {
                return ParseResult::Reported(format!("int_of_sexp: invalid integer {text:?}"));
            }
        }
    }
    ParseResult::Parsed(Electrons {
        elec_alpha_num: values[0],
        elec_beta_num: values[1],
    })
}

/// `set` for the section: parse the body and, when `write_back` is true,
/// write it.  The silent retry of `qp_edit.set` (wave S2.4) is in
/// [`super::apply_parsed`].
pub(crate) fn set(ezfio: &EzFio, title: &str, body: &str, write_back: bool) -> Result<()> {
    apply_parsed(title, body, write_back, of_rst, |data| write(ezfio, data))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::make_header;

    fn scratch() -> (tempfile::TempDir, EzFio) {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("ezfio")).expect("create EZFIO");
        (dir, ezfio)
    }

    #[test]
    fn read_returns_none_without_the_alpha_counter() {
        let (_dir, ezfio) = scratch();
        assert_eq!(read(&ezfio).expect("read"), None);
        // The guard is on alpha only: a stray beta alone is still invisible.
        ezfio
            .set::<i64>(GROUP, "elec_beta_num", &3)
            .expect("set beta");
        assert_eq!(read(&ezfio).expect("read"), None);
    }

    #[test]
    fn to_rst_matches_the_ocaml_layout() {
        let data = Electrons {
            elec_alpha_num: 7,
            elec_beta_num: 7,
        };
        assert_eq!(
            to_rst(&data).expect("render"),
            "\nSpin multiplicity is Singlet.\n\nNumber of alpha and beta electrons ::\n\n  \
             elec_alpha_num = 7\n  elec_beta_num  = 7\n\n"
        );
        let doublet = Electrons {
            elec_alpha_num: 2,
            elec_beta_num: 1,
        };
        assert!(
            to_rst(&doublet)
                .expect("render")
                .contains("Spin multiplicity is Doublet.")
        );
    }

    #[test]
    fn get_prefixes_the_section_header_and_hides_an_absent_section() {
        let (_dir, ezfio) = scratch();
        assert_eq!(get(&ezfio, "Electrons").expect("empty"), "");

        ezfio
            .set::<i64>(GROUP, "elec_alpha_num", &5)
            .expect("set alpha");
        ezfio
            .set::<i64>(GROUP, "elec_beta_num", &4)
            .expect("set beta");
        assert_eq!(
            get(&ezfio, "Electrons").expect("render"),
            make_header("Electrons")
                + "\nSpin multiplicity is Doublet.\n\nNumber of alpha and beta electrons ::\n\n  \
                  elec_alpha_num = 5\n  elec_beta_num  = 4\n\n"
        );
    }

    #[test]
    fn of_rst_ignores_prose_lines_and_accepts_the_rendered_body() {
        let data = Electrons {
            elec_alpha_num: 7,
            elec_beta_num: 3,
        };
        let body = to_rst(&data).expect("render");
        assert!(
            matches!(of_rst(&body), ParseResult::Parsed(parsed) if parsed == data),
            "expected Parsed({data:?})"
        );
    }

    #[test]
    fn of_rst_reports_a_malformed_value_instead_of_parsing_it() {
        // Wave S2.4: `of_rst` no longer prints; the driver reports the
        // `Reported` message as `qp_edit: parse error in <title>: <msg>`.
        let malformed = of_rst("  elec_alpha_num = seven\n  elec_beta_num = 7\n");
        assert!(
            matches!(&malformed, ParseResult::Reported(message) if message.contains("int_of_sexp: invalid integer")),
            "{malformed:?}"
        );
        let missing = of_rst("  elec_alpha_num = 7\n");
        assert!(
            matches!(&missing, ParseResult::Reported(message) if message.contains("missing field: elec_beta_num")),
            "{missing:?}"
        );
    }

    #[test]
    fn set_applies_only_when_asked() {
        let (_dir, ezfio) = scratch();
        ezfio
            .set::<i64>(GROUP, "elec_alpha_num", &7)
            .expect("set alpha");
        ezfio
            .set::<i64>(GROUP, "elec_beta_num", &7)
            .expect("set beta");
        let body = "  elec_alpha_num = 6\n  elec_beta_num = 5\n";

        set(&ezfio, "Electrons", body, false).expect("parse only");
        assert_eq!(ezfio.get::<i64>(GROUP, "elec_alpha_num").expect("alpha"), 7);
        assert_eq!(ezfio.get::<i64>(GROUP, "elec_beta_num").expect("beta"), 7);

        set(&ezfio, "Electrons", body, true).expect("apply");
        assert_eq!(ezfio.get::<i64>(GROUP, "elec_alpha_num").expect("alpha"), 6);
        assert_eq!(ezfio.get::<i64>(GROUP, "elec_beta_num").expect("beta"), 5);
    }

    #[test]
    fn a_section_body_without_its_header_is_a_no_op() {
        let (_dir, ezfio) = scratch();
        // `set` is called with the whole RST by the driver; the section text
        // itself does not contain its header, which is expected: the driver
        // extracts the body.  This test pins the body parser, not the header.
        set(&ezfio, "Electrons", "", true).expect("empty body");
        assert!(!ezfio.has(GROUP, "elec_alpha_num").expect("has"));
    }

    #[test]
    fn read_rejects_counts_that_the_ocaml_wrappers_assert() {
        let (_dir, ezfio) = scratch();
        ezfio.set::<i64>(GROUP, "elec_alpha_num", &0).expect("set");
        ezfio.set::<i64>(GROUP, "elec_beta_num", &0).expect("set");
        let error = read(&ezfio).expect_err("alpha must be positive");
        assert!(error.to_string().contains("x > 0"), "{error}");

        ezfio.set::<i64>(GROUP, "elec_alpha_num", &1).expect("set");
        ezfio.set::<i64>(GROUP, "elec_beta_num", &-1).expect("set");
        let error = read(&ezfio).expect_err("beta must be non-negative");
        assert!(error.to_string().contains("x >= 0"), "{error}");
    }
}
