//! Declared-type table: how an `EZFIO.cfg` `type:` string is parsed and
//! validated.
//!
//! `EZFIO.cfg` declares types by name (`integer`, `double precision`,
//! `Threshold`, `MO_class`, ...).  The Python generator deliberately stores
//! the declared string **verbatim** in the TOML config, so this table is the
//! single place that knows what each name means.  It mirrors
//! `scripts/ezfio_interface/ei_handler.py:get_type_dict`, which reads
//! `ocaml/qptypes_generator.ml`; every spelling that generator accepts has an
//! entry here, so an `EZFIO.cfg` edit never requires recompiling the binary.
//!
//! Units: none.  The table stores storage/validation rules only; physical
//! units are not declared anywhere in `EZFIO.cfg` and are not invented here.
//!
//! Constraints are copied from `qptypes_generator.ml`, which uses
//! `assert`/`raise (Invalid_argument ...)`.  Two deliberate differences are
//! documented per entry: plain `integer` is range-checked to 32 bits (its
//! Fortran storage width; the OCaml side uses a 64-bit `int`), and the editor
//! rejects values that OCaml silently clamps (`MO_occ`).

use std::ops::Bound;

use anyhow::{Result, bail};

/// Primitive representation of one EZFIO element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValueKind {
    /// Signed 64-bit integer.  Covers `integer` (32-bit storage, range
    /// checked), `integer*8`, `int`, `int64` and the integer wrapper types.
    Int,
    /// IEEE-754 double precision (`double precision`, `float`, and the float
    /// wrapper types).
    Float,
    /// Fortran `logical`: the literals `true` and `false`, case-insensitive.
    /// EZFIO's `T`/`F` rendering belongs to the P1 I/O layer, not here.
    Bool,
    /// Character string (`character*(N)` and the string wrapper types).
    String,
}

impl ValueKind {
    /// Short human name used by `qp_edit --describe`.
    pub const fn name(self) -> &'static str {
        match self {
            ValueKind::Int => "int",
            ValueKind::Float => "float",
            ValueKind::Bool => "bool",
            ValueKind::String => "string",
        }
    }
}

/// One parsed value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// Integer element.
    Int(i64),
    /// Floating-point element.
    Float(f64),
    /// Logical element.
    Bool(bool),
    /// String element (fixed-width character types included).
    Str(String),
}

/// Accepted domain of a type, checked after the primitive parse.
#[derive(Debug, Clone, Copy)]
pub enum Domain {
    /// No constraint beyond the primitive parse.
    Any,
    /// Inclusive/exclusive integer bounds.
    IntRange {
        /// Lower bound.
        min: Bound<i64>,
        /// Upper bound.
        max: Bound<i64>,
    },
    /// Inclusive/exclusive floating-point bounds.
    FloatRange {
        /// Lower bound.
        min: Bound<f64>,
        /// Upper bound.
        max: Bound<f64>,
    },
    /// Case-insensitive choice among fixed spellings.
    Choices(&'static [&'static str]),
    /// Case-insensitive choice among fixed integers.
    IntChoices(&'static [i64]),
    /// `MO_class`: `Core|Inactive|Active|Virtual|Deleted`, optionally followed
    /// by the `: (i, j-k, ...)` list printed by `MO_class.to_string`.
    MoClass,
    /// 32-character lowercase-hex-ish key, exactly as `qptypes_generator.ml`
    /// checks it (`length = 32` and every byte `< 'g'`).
    Md5,
    /// String that must not be empty.
    NonEmptyString,
}

/// Parsing and validation rules for one declared type.
#[derive(Debug, Clone, Copy)]
pub struct TypeSpec {
    declared: &'static str,
    value: ValueKind,
    domain: Domain,
}

impl TypeSpec {
    /// Primitive kind used to parse and print one element.
    pub const fn value(&self) -> ValueKind {
        self.value
    }

    /// Validation domain applied after the primitive parse.
    pub const fn domain(&self) -> Domain {
        self.domain
    }

    /// Human description of the accepted domain (`"int in [1, 10000]"`,
    /// `"one of Read|Write|None"`, ...).  Used by `--describe`; allocation is
    /// acceptable because this is not a hot path.
    pub fn describe_domain(&self) -> String {
        match self.domain {
            Domain::Any => self.value.name().to_string(),
            Domain::IntRange { min, max } => {
                format!("int in {}", int_range(min, max))
            }
            Domain::FloatRange { min, max } => {
                format!("float in {}", float_range(min, max))
            }
            Domain::Choices(choices) => format!("one of {}", choices.join("|")),
            Domain::IntChoices(choices) => {
                let rendered: Vec<String> = choices.iter().map(i64::to_string).collect();
                format!("one of {}", rendered.join("|"))
            }
            Domain::MoClass => "MO_class keyword".to_string(),
            Domain::Md5 => "32-character MD5".to_string(),
            Domain::NonEmptyString => "non-empty string".to_string(),
        }
    }

    /// Parse and validate one scalar element.
    ///
    /// `raw` is the textual value as it appears in EZFIO or in an `EZFIO.cfg`
    /// `default:`.  Whitespace around the value is ignored.  Errors carry the
    /// offending value and the declared type so the editor can point at the
    /// field.
    pub fn parse_scalar(&self, raw: &str) -> Result<Value> {
        let text = raw.trim();
        let value = match self.value {
            ValueKind::Int => Value::Int(text.parse::<i64>().map_err(|_| {
                anyhow::anyhow!(
                    "{text:?} is not a valid integer for type {:?}",
                    self.declared
                )
            })?),
            ValueKind::Float => Value::Float(text.parse::<f64>().map_err(|_| {
                anyhow::anyhow!("{text:?} is not a valid float for type {:?}", self.declared)
            })?),
            ValueKind::Bool => match text.to_ascii_lowercase().as_str() {
                "true" => Value::Bool(true),
                "false" => Value::Bool(false),
                _ => bail!(
                    "{text:?} is not a valid logical for type {:?}",
                    self.declared
                ),
            },
            ValueKind::String => Value::Str(text.to_string()),
        };
        self.check_domain(&value)?;
        Ok(value)
    }

    /// Check an already parsed value against the domain.
    fn check_domain(&self, value: &Value) -> Result<()> {
        let accepted = match (self.domain, value) {
            (Domain::Any, _) => true,
            (Domain::IntRange { min, max }, Value::Int(v)) => bound_check(*v, min, max, i64::cmp),
            (Domain::FloatRange { min, max }, Value::Float(v)) => {
                // NaN compares false against every bound: reject it there.
                !v.is_nan() && bound_check(*v, min, max, f64::total_cmp)
            }
            (Domain::Choices(choices), Value::Str(s)) => {
                choices.iter().any(|c| c.eq_ignore_ascii_case(s))
            }
            (Domain::IntChoices(choices), Value::Int(v)) => choices.contains(v),
            (Domain::MoClass, Value::Str(s)) => is_mo_class(s),
            (Domain::Md5, Value::Str(s)) => s.len() == 32 && s.bytes().all(|b| b < b'g'),
            (Domain::NonEmptyString, Value::Str(s)) => !s.is_empty(),
            (domain, value) => bail!(
                "internal type-table error: domain {domain:?} does not apply to value {value:?}"
            ),
        };
        if accepted {
            Ok(())
        } else {
            bail!(
                "value {:?} is outside the domain of type {:?} ({})",
                match value {
                    Value::Int(v) => v.to_string(),
                    Value::Float(v) => v.to_string(),
                    Value::Bool(v) => v.to_string(),
                    Value::Str(s) => s.clone(),
                },
                self.declared,
                self.describe_domain()
            )
        }
    }
}

/// Inclusive/exclusive bound test.  `cmp` is passed explicitly so the float
/// path can use a total order (`f64::total_cmp`) and stay free of NaN
/// surprises.
fn bound_check<T: Copy>(
    value: T,
    min: Bound<T>,
    max: Bound<T>,
    cmp: fn(&T, &T) -> std::cmp::Ordering,
) -> bool {
    use std::cmp::Ordering;

    let above_min = match min {
        Bound::Included(b) => cmp(&value, &b) != Ordering::Less,
        Bound::Excluded(b) => cmp(&value, &b) == Ordering::Greater,
        Bound::Unbounded => true,
    };
    let below_max = match max {
        Bound::Included(b) => cmp(&value, &b) != Ordering::Greater,
        Bound::Excluded(b) => cmp(&value, &b) == Ordering::Less,
        Bound::Unbounded => true,
    };
    above_min && below_max
}

/// Render an integer interval in mathematical notation.
fn int_range(min: Bound<i64>, max: Bound<i64>) -> String {
    let left = match min {
        Bound::Included(v) => format!("[{v}"),
        Bound::Excluded(v) => format!("({v}"),
        Bound::Unbounded => "(-inf".to_string(),
    };
    let right = match max {
        Bound::Included(v) => format!("{v}]"),
        Bound::Excluded(v) => format!("{v})"),
        Bound::Unbounded => "+inf)".to_string(),
    };
    format!("{left}, {right}")
}

/// Render a floating-point interval in mathematical notation.
fn float_range(min: Bound<f64>, max: Bound<f64>) -> String {
    let left = match min {
        Bound::Included(v) => format!("[{v}"),
        Bound::Excluded(v) => format!("({v}"),
        Bound::Unbounded => "(-inf".to_string(),
    };
    let right = match max {
        Bound::Included(v) => format!("{v}]"),
        Bound::Excluded(v) => format!("{v})"),
        Bound::Unbounded => "+inf)".to_string(),
    };
    format!("{left}, {right}")
}

/// `MO_class` grammar: a case-insensitive keyword, then optionally the list
/// printed by `MO_class.to_string` (`"Active   : (1, 2)"`).  Only the keyword
/// is validated: the list is the `mo_basis` handler's business (P4).
fn is_mo_class(raw: &str) -> bool {
    let keyword = raw.split(':').next().unwrap_or("").trim();
    matches!(
        keyword.to_ascii_lowercase().as_str(),
        "core" | "inactive" | "active" | "virtual" | "deleted"
    )
}

use Domain::{Any, Choices, FloatRange, IntChoices, IntRange, Md5, MoClass, NonEmptyString};
use ValueKind::{Bool, Float, Int, String as StrKind};

/// Every type spelling accepted by `ei_handler.get_type_dict`, sorted by
/// [`TypeSpec::declared`] (bytewise, which is also `str` order).  Keeping the
/// table sorted lets [`lookup_type`] binary-search it; a unit test enforces
/// the order.
///
/// `%s` is present in the Python dictionary as a parsing artifact of
/// `qptypes_generator.ml` and is deliberately omitted: it is not a usable
/// type.
pub const TYPES: &[TypeSpec] = &[
    // `AO_basis_name`: non-empty basis-set name.
    TypeSpec {
        declared: "AO_basis_name",
        value: StrKind,
        domain: NonEmptyString,
    },
    TypeSpec {
        declared: "AO_coef",
        value: Float,
        domain: Any,
    },
    // `AO_expo`: x >= 0.
    TypeSpec {
        declared: "AO_expo",
        value: Float,
        domain: FloatRange {
            min: Bound::Included(0.0),
            max: Bound::Unbounded,
        },
    },
    // `AO_number`: 1 <= x <= 10_000 (input_ezfio).
    TypeSpec {
        declared: "AO_number",
        value: Int,
        domain: IntRange {
            min: Bound::Included(1),
            max: Bound::Included(10_000),
        },
    },
    // `AO_prim_number`: x > 0.
    TypeSpec {
        declared: "AO_prim_number",
        value: Int,
        domain: IntRange {
            min: Bound::Included(1),
            max: Bound::Unbounded,
        },
    },
    // `Bit_kind`: 1 | 2 | 4 | 8.
    TypeSpec {
        declared: "Bit_kind",
        value: Int,
        domain: IntChoices(&[1, 2, 4, 8]),
    },
    // `Bit_kind_size`: 8 | 16 | 32 | 64.
    TypeSpec {
        declared: "Bit_kind_size",
        value: Int,
        domain: IntChoices(&[8, 16, 32, 64]),
    },
    // `Det_coef`: -1 <= x <= 1.
    TypeSpec {
        declared: "Det_coef",
        value: Float,
        domain: FloatRange {
            min: Bound::Included(-1.0),
            max: Bound::Included(1.0),
        },
    },
    // `Det_number`: 1 <= x <= 50_000_000_000 (input_ezfio).
    TypeSpec {
        declared: "Det_number",
        value: Int,
        domain: IntRange {
            min: Bound::Included(1),
            max: Bound::Included(50_000_000_000),
        },
    },
    // `Det_number_max`: x > 0; OCaml only warns above 50 billion.
    TypeSpec {
        declared: "Det_number_max",
        value: Int,
        domain: IntRange {
            min: Bound::Included(1),
            max: Bound::Unbounded,
        },
    },
    // `Disk_access`: Read | Write | None, case-insensitive.
    TypeSpec {
        declared: "Disk_access",
        value: StrKind,
        domain: Choices(&["Read", "Write", "None"]),
    },
    // `Elec_alpha_number`: x > 0.
    TypeSpec {
        declared: "Elec_alpha_number",
        value: Int,
        domain: IntRange {
            min: Bound::Included(1),
            max: Bound::Unbounded,
        },
    },
    // `Elec_beta_number`: x >= 0.
    TypeSpec {
        declared: "Elec_beta_number",
        value: Int,
        domain: IntRange {
            min: Bound::Included(0),
            max: Bound::Unbounded,
        },
    },
    // `Elec_number`: x > 0.
    TypeSpec {
        declared: "Elec_number",
        value: Int,
        domain: IntRange {
            min: Bound::Included(1),
            max: Bound::Unbounded,
        },
    },
    // `Energy`: x <= 0.
    TypeSpec {
        declared: "Energy",
        value: Float,
        domain: FloatRange {
            min: Bound::Unbounded,
            max: Bound::Included(0.0),
        },
    },
    // `MD5`: 32 characters, every byte < 'g'.
    TypeSpec {
        declared: "MD5",
        value: StrKind,
        domain: Md5,
    },
    // `MO_class`: keyword grammar; the list suffix is parsed by the handler.
    TypeSpec {
        declared: "MO_class",
        value: StrKind,
        domain: MoClass,
    },
    TypeSpec {
        declared: "MO_coef",
        value: Float,
        domain: Any,
    },
    // `MO_guess`: Huckel | HCore, case-insensitive.
    TypeSpec {
        declared: "MO_guess",
        value: StrKind,
        domain: Choices(&["Huckel", "HCore"]),
    },
    // `MO_number`: 1 <= x <= 10_000 (input_ezfio).
    TypeSpec {
        declared: "MO_number",
        value: Int,
        domain: IntRange {
            min: Bound::Included(1),
            max: Bound::Included(10_000),
        },
    },
    // `MO_occ`: OCaml clamps to [0, 2]; the editor rejects out-of-range
    // values instead of silently changing them (D4: no read-side re-stamp).
    TypeSpec {
        declared: "MO_occ",
        value: Float,
        domain: FloatRange {
            min: Bound::Included(0.0),
            max: Bound::Included(2.0),
        },
    },
    // `N_int_number`: 1 <= x <= 128 (input_ezfio).
    TypeSpec {
        declared: "N_int_number",
        value: Int,
        domain: IntRange {
            min: Bound::Included(1),
            max: Bound::Included(128),
        },
    },
    // `Negative_float`: x <= 0.
    TypeSpec {
        declared: "Negative_float",
        value: Float,
        domain: FloatRange {
            min: Bound::Unbounded,
            max: Bound::Included(0.0),
        },
    },
    // `Negative_int`: x <= 0.
    TypeSpec {
        declared: "Negative_int",
        value: Int,
        domain: IntRange {
            min: Bound::Unbounded,
            max: Bound::Included(0),
        },
    },
    // `Non_empty_string`: x <> "".
    TypeSpec {
        declared: "Non_empty_string",
        value: StrKind,
        domain: NonEmptyString,
    },
    // `Normalized_float`: 0 <= x <= 1.
    TypeSpec {
        declared: "Normalized_float",
        value: Float,
        domain: FloatRange {
            min: Bound::Included(0.0),
            max: Bound::Included(1.0),
        },
    },
    // `Nucl_number`: 1 <= x <= 10_000 (input_ezfio).
    TypeSpec {
        declared: "Nucl_number",
        value: Int,
        domain: IntRange {
            min: Bound::Included(1),
            max: Bound::Included(10_000),
        },
    },
    // `PT2_energy`: x >= 0.
    TypeSpec {
        declared: "PT2_energy",
        value: Float,
        domain: FloatRange {
            min: Bound::Included(0.0),
            max: Bound::Unbounded,
        },
    },
    // `Perturbation`: EN | HF | CFG, case-insensitive.
    TypeSpec {
        declared: "Perturbation",
        value: StrKind,
        domain: Choices(&["EN", "HF", "CFG"]),
    },
    // `Positive_float`: x >= 0.
    TypeSpec {
        declared: "Positive_float",
        value: Float,
        domain: FloatRange {
            min: Bound::Included(0.0),
            max: Bound::Unbounded,
        },
    },
    // `Positive_int`: x >= 0.
    TypeSpec {
        declared: "Positive_int",
        value: Int,
        domain: IntRange {
            min: Bound::Included(0),
            max: Bound::Unbounded,
        },
    },
    // `Positive_int64`: x >= 0, no upper bound.
    TypeSpec {
        declared: "Positive_int64",
        value: Int,
        domain: IntRange {
            min: Bound::Included(0),
            max: Bound::Unbounded,
        },
    },
    // `R_power`: -2 <= x <= 8.
    TypeSpec {
        declared: "R_power",
        value: Int,
        domain: IntRange {
            min: Bound::Included(-2),
            max: Bound::Included(8),
        },
    },
    TypeSpec {
        declared: "Rst_string",
        value: StrKind,
        domain: Any,
    },
    // `S2`: x >= 0.
    TypeSpec {
        declared: "S2",
        value: Float,
        domain: FloatRange {
            min: Bound::Included(0.0),
            max: Bound::Unbounded,
        },
    },
    // `States_number`: x > 0; OCaml only warns above 1000.
    TypeSpec {
        declared: "States_number",
        value: Int,
        domain: IntRange {
            min: Bound::Included(1),
            max: Bound::Unbounded,
        },
    },
    // `Strictly_negative_float`: x < 0.
    TypeSpec {
        declared: "Strictly_negative_float",
        value: Float,
        domain: FloatRange {
            min: Bound::Unbounded,
            max: Bound::Excluded(0.0),
        },
    },
    // `Strictly_negative_int`: x < 0.
    TypeSpec {
        declared: "Strictly_negative_int",
        value: Int,
        domain: IntRange {
            min: Bound::Unbounded,
            max: Bound::Excluded(0),
        },
    },
    // `Strictly_positive_float`: x > 0.
    TypeSpec {
        declared: "Strictly_positive_float",
        value: Float,
        domain: FloatRange {
            min: Bound::Excluded(0.0),
            max: Bound::Unbounded,
        },
    },
    // `Strictly_positive_int`: x > 0.
    TypeSpec {
        declared: "Strictly_positive_int",
        value: Int,
        domain: IntRange {
            min: Bound::Excluded(0),
            max: Bound::Unbounded,
        },
    },
    // `Threshold`: 0 <= x <= 1.
    TypeSpec {
        declared: "Threshold",
        value: Float,
        domain: FloatRange {
            min: Bound::Included(0.0),
            max: Bound::Included(1.0),
        },
    },
    // `bool`: alias of `logical`.
    TypeSpec {
        declared: "bool",
        value: Bool,
        domain: Any,
    },
    // `character*(256)`: free text; width is not enforced (OCaml does not
    // check it either, and EZFIO truncates on write).
    TypeSpec {
        declared: "character*(256)",
        value: StrKind,
        domain: Any,
    },
    TypeSpec {
        declared: "character*(32)",
        value: StrKind,
        domain: Any,
    },
    TypeSpec {
        declared: "character*(64)",
        value: StrKind,
        domain: Any,
    },
    TypeSpec {
        declared: "double precision",
        value: Float,
        domain: Any,
    },
    TypeSpec {
        declared: "float",
        value: Float,
        domain: Any,
    },
    // `int`: alias of `integer`; 32-bit Fortran storage.
    TypeSpec {
        declared: "int",
        value: Int,
        domain: IntRange {
            min: Bound::Included(i32::MIN as i64),
            max: Bound::Included(i32::MAX as i64),
        },
    },
    // `int64`: 64-bit storage.
    TypeSpec {
        declared: "int64",
        value: Int,
        domain: Any,
    },
    // `integer`: 32-bit Fortran storage.  The OCaml side uses a 64-bit `int`;
    // the editor is stricter on purpose (the Fortran writer would overflow).
    TypeSpec {
        declared: "integer",
        value: Int,
        domain: IntRange {
            min: Bound::Included(i32::MIN as i64),
            max: Bound::Included(i32::MAX as i64),
        },
    },
    // `integer*8`: 64-bit storage.
    TypeSpec {
        declared: "integer*8",
        value: Int,
        domain: Any,
    },
    // `logical`: alias of `bool`.
    TypeSpec {
        declared: "logical",
        value: Bool,
        domain: Any,
    },
];

/// Look up a declared type by its exact `EZFIO.cfg` spelling.
///
/// Returns `None` for a spelling that `ei_handler.get_type_dict` does not
/// know; the loader turns that into an error naming the module and field.
pub fn lookup_type(declared: &str) -> Option<&'static TypeSpec> {
    TYPES
        .binary_search_by(|spec| spec.declared.cmp(declared))
        .ok()
        .map(|index| &TYPES[index])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_sorted_and_unique() {
        for pair in TYPES.windows(2) {
            assert!(
                pair[0].declared < pair[1].declared,
                "type table not sorted: {:?} then {:?}",
                pair[0].declared,
                pair[1].declared
            );
        }
    }

    #[test]
    fn lookup_finds_every_entry() {
        for spec in TYPES {
            assert_eq!(
                lookup_type(spec.declared).map(|found| found.declared),
                Some(spec.declared)
            );
        }
        assert!(lookup_type("no_such_type").is_none());
        assert!(lookup_type("%s").is_none());
    }

    #[test]
    fn threshold_domain_is_closed() {
        let spec = lookup_type("Threshold").expect("Threshold in table");
        assert_eq!(
            spec.parse_scalar("0.").expect("lower bound"),
            Value::Float(0.0)
        );
        assert_eq!(
            spec.parse_scalar("1.e-14").expect("inside"),
            Value::Float(1e-14)
        );
        assert_eq!(
            spec.parse_scalar("1.0").expect("upper bound"),
            Value::Float(1.0)
        );
        assert!(spec.parse_scalar("1.5").is_err());
        assert!(spec.parse_scalar("-1.e-3").is_err());
    }

    #[test]
    fn integer_storage_is_32_bit() {
        let spec = lookup_type("integer").expect("integer in table");
        assert!(spec.parse_scalar("2147483647").is_ok());
        assert!(spec.parse_scalar("2147483648").is_err());
        assert!(spec.parse_scalar("-1").is_ok());
        let wide = lookup_type("integer*8").expect("integer*8 in table");
        assert!(wide.parse_scalar("2147483648").is_ok());
    }

    #[test]
    fn strict_bounds_exclude_zero() {
        let strict = lookup_type("Strictly_positive_int").expect("in table");
        assert!(strict.parse_scalar("0").is_err());
        assert!(strict.parse_scalar("1").is_ok());
        let strict_float = lookup_type("Strictly_negative_float").expect("in table");
        assert!(strict_float.parse_scalar("0.0").is_err());
        assert!(strict_float.parse_scalar("-1.e-30").is_ok());
    }

    #[test]
    fn choices_are_case_insensitive() {
        let spec = lookup_type("Disk_access").expect("in table");
        assert_eq!(
            spec.parse_scalar("None").expect("canonical"),
            Value::Str("None".into())
        );
        assert_eq!(
            spec.parse_scalar("read").expect("lowercase"),
            Value::Str("read".into())
        );
        assert!(spec.parse_scalar("append").is_err());
    }

    #[test]
    fn mo_class_accepts_keyword_and_printed_form() {
        let spec = lookup_type("MO_class").expect("in table");
        assert!(spec.parse_scalar("Active").is_ok());
        assert!(spec.parse_scalar("virtual").is_ok());
        assert!(spec.parse_scalar("Active   : (1, 3-5)").is_ok());
        assert!(spec.parse_scalar("Frozen").is_err());
    }

    #[test]
    fn nan_is_rejected() {
        let spec = lookup_type("Positive_float").expect("in table");
        assert!(spec.parse_scalar("NaN").is_err());
        assert!(spec.parse_scalar("inf").is_ok());
    }
}
