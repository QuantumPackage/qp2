//! Static registry of the EZFIO attributes declared by the QP2 modules.
//!
//! Each `src/<module>/EZFIO.cfg` file declares sections such as
//!
//! ```text
//! [psi_det]
//! type: integer*8
//! size: (determinants.n_int*determinants.bit_kind/8,2,determinants.n_det)
//! ```
//!
//! [`crate::keys::KEYS`] is the generated table of the 332 attributes of the
//! 43 tracked modules, sorted by `(group, key)` so that lookup is a binary
//! search. It is **advisory**: EZFIO directories written by plugins contain
//! keys this table does not know, and both the Fortran and the OCaml clients
//! read and write such keys dynamically. The typed API in [`crate::EzFio`]
//! therefore only checks the declared type and shape when the key is known,
//! and trusts the call site otherwise.
//!
//! Dimensions are stored as the raw `size:` expressions and parsed on demand
//! ([`DimExpr`]), because they routinely mix integer literals, runtime
//! references (`determinants.n_det`), arithmetic and Fortran ranges
//! (`0:pseudo.pseudo_lmax`). Parsing one expression allocates a handful of
//! `Box`es and no strings; a typed access touches at most a few expressions.

use crate::error::{Error, Result};
use crate::layout::ElementType;

/// One attribute declared in an `EZFIO.cfg` file.
///
/// The fields mirror the generated `<module>.ezfio_interface_config` files:
/// `group` is the `ezfio_dir` (module directory name, lower-cased), `key` is
/// the `ezfio_name` (section name unless renamed, lower-cased), and `dims`
/// holds the comma-separated `size:` expressions with the outer parentheses
/// and surrounding blanks removed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KeySpec {
    /// EZFIO group: the directory name, e.g. `"determinants"`.
    pub group: &'static str,
    /// Attribute name inside the group, e.g. `"psi_det"`.
    pub key: &'static str,
    /// Element type declared by the `type:` field.
    pub ty: ElementType,
    /// Declared length of a `character*(n)` attribute, `0` for other types.
    pub char_len: usize,
    /// Raw dimension expressions of an array, `None` for a scalar.
    ///
    /// Empty is not possible: a `size:` of `1` means scalar, and the generator
    /// rejects an empty dimension list.
    pub dims: Option<&'static [&'static str]>,
}

impl KeySpec {
    /// Build a spec. `const fn`, so that the generated table in
    /// [`crate::keys::KEYS`] can be a `static`.
    ///
    /// - `group`, `key`: EZFIO directory and attribute names.
    /// - `ty`: declared element type.
    /// - `char_len`: declared `character*(n)` length, `0` for other types.
    /// - `dims`: raw `size:` expressions, `None` for a scalar.
    pub const fn new(
        group: &'static str,
        key: &'static str,
        ty: ElementType,
        char_len: usize,
        dims: Option<&'static [&'static str]>,
    ) -> Self {
        Self {
            group,
            key,
            ty,
            char_len,
            dims,
        }
    }

    /// True when the attribute is a gzip array rather than a plain-text
    /// scalar.
    pub fn is_array(&self) -> bool {
        self.dims.is_some()
    }
}

/// A source of [`KeySpec`]s.
///
/// Implemented by [`StaticRegistry`] (the generated table). The trait exists
/// so that a future code-generated registry, or a test double, can be plugged
/// into [`crate::EzFio`] without touching the I/O layer.
pub trait Registry: std::fmt::Debug {
    /// The spec of `group/key`, or `None` if the registry does not know it.
    ///
    /// `group` and `key` are matched case-sensitively; EZFIO names are
    /// lower-case by construction (the module configs are lower-cased when the
    /// interface files are generated).
    fn lookup(&self, group: &str, key: &str) -> Option<KeySpec>;

    /// Every spec of the registry, sorted by `(group, key)`.
    fn specs(&self) -> &'static [KeySpec];
}

/// The registry compiled into this crate: [`crate::keys::KEYS`].
#[derive(Debug, Clone, Copy, Default)]
pub struct StaticRegistry;

impl StaticRegistry {
    /// Build the registry. It is stateless, so this is only sugar for
    /// `StaticRegistry` in expression position.
    pub const fn new() -> Self {
        Self
    }
}

impl Registry for StaticRegistry {
    fn lookup(&self, group: &str, key: &str) -> Option<KeySpec> {
        let keys = crate::keys::KEYS;
        keys.binary_search_by(|spec| (spec.group, spec.key).cmp(&(group, key)))
            .ok()
            .map(|index| keys[index])
    }

    fn specs(&self) -> &'static [KeySpec] {
        crate::keys::KEYS
    }
}

/// Resolves the runtime symbols of a dimension expression.
///
/// EZFIO dimension expressions may refer to other attributes, always scalar
/// integers in the QP2 module configs (`determinants.n_int`,
/// `pseudo.pseudo_lmax`, ...). A reference to an attribute that is not set is
/// *unresolved*, not an error: the shape of the array then simply cannot be
/// checked (e.g. `pseudo/pseudo_dz_k` is written with dimensions `(3,0)`
/// before `pseudo_klocmax` exists).
pub trait DimResolver {
    /// Value of the scalar integer attribute `group/key`.
    ///
    /// Returns `Ok(None)` when the attribute is not set, and an error when it
    /// is set but is not a plain-text integer.
    fn resolve(&self, group: &str, key: &str) -> Result<Option<i64>>;
}

/// A dimension expression, as written in the `size:` field of `EZFIO.cfg`.
///
/// Grammar (all tokens ASCII, blanks allowed between them):
///
/// ```text
/// dim    := expr (':' expr)?          Fortran range lo:hi, extent hi-lo+1
/// expr   := term (('+'|'-') term)*
/// term   := factor (('*'|'/') factor)*
/// factor := '-' factor | '+' factor | '(' expr ')' | integer | ref
/// ref    := identifier '.' identifier
/// ```
///
/// Division is Fortran integer division (truncation toward zero) and every
/// operation is checked for overflow. The corpus exercises all of it:
/// `determinants.n_int*determinants.bit_kind/8`,
/// `-pseudo.pseudo_lmax:pseudo.pseudo_lmax` and `0:pseudo.pseudo_lmax`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DimExpr<'a> {
    /// Integer literal, e.g. `3`.
    Int(i64),
    /// Reference to another attribute, e.g. `determinants.n_det`.
    Ref {
        /// EZFIO group of the referenced attribute.
        group: &'a str,
        /// Attribute name inside that group.
        key: &'a str,
    },
    /// `a + b`.
    Add(Box<DimExpr<'a>>, Box<DimExpr<'a>>),
    /// `a - b`.
    Sub(Box<DimExpr<'a>>, Box<DimExpr<'a>>),
    /// `a * b`.
    Mul(Box<DimExpr<'a>>, Box<DimExpr<'a>>),
    /// `a / b`, Fortran integer division.
    Div(Box<DimExpr<'a>>, Box<DimExpr<'a>>),
    /// Unary `-a`.
    Neg(Box<DimExpr<'a>>),
    /// Fortran range `lo:hi`, whose extent is `hi - lo + 1`.
    Range {
        /// Lower bound (inclusive).
        lo: Box<DimExpr<'a>>,
        /// Upper bound (inclusive).
        hi: Box<DimExpr<'a>>,
    },
}

impl<'a> DimExpr<'a> {
    /// Parse one dimension expression.
    ///
    /// The returned tree borrows from `expr`; only the interior nodes are
    /// boxed. Fails with [`Error::Dim`] on a syntax error or an out-of-range
    /// integer literal.
    pub fn parse(expr: &'a str) -> Result<Self> {
        let mut parser = Parser::new(expr);
        let lo = parser.parse_sum()?;
        parser.skip_whitespace();
        if parser.peek() == Some(b':') {
            parser.bump();
            let hi = parser.parse_sum()?;
            parser.skip_whitespace();
            if parser.peek().is_some() {
                return Err(parser.error("unexpected trailing characters after the range"));
            }
            Ok(DimExpr::Range {
                lo: Box::new(lo),
                hi: Box::new(hi),
            })
        } else {
            if parser.peek().is_some() {
                return Err(parser.error("unexpected trailing characters"));
            }
            Ok(lo)
        }
    }

    /// Evaluate the expression.
    ///
    /// `expr` is the source text, used only in error messages. Returns
    /// `Ok(None)` when a referenced attribute is not set; that is not an error
    /// (see [`DimResolver`]). Fails with [`Error::Dim`] on division by zero or
    /// on `i64` overflow.
    pub fn eval(&self, expr: &str, resolver: &dyn DimResolver) -> Result<Option<i64>> {
        match self {
            DimExpr::Int(value) => Ok(Some(*value)),
            DimExpr::Ref { group, key } => resolver.resolve(group, key),
            DimExpr::Neg(inner) => match inner.eval(expr, resolver)? {
                None => Ok(None),
                Some(value) => value.checked_neg().map(Some).ok_or_else(|| Error::Dim {
                    expr: expr.to_string(),
                    detail: format!("-{value} overflows i64"),
                }),
            },
            DimExpr::Add(a, b) => eval_binary(expr, resolver, a, b, i64::checked_add, "+"),
            DimExpr::Sub(a, b) => eval_binary(expr, resolver, a, b, i64::checked_sub, "-"),
            DimExpr::Mul(a, b) => eval_binary(expr, resolver, a, b, i64::checked_mul, "*"),
            DimExpr::Div(a, b) => {
                let (Some(x), Some(y)) = (a.eval(expr, resolver)?, b.eval(expr, resolver)?) else {
                    return Ok(None);
                };
                if y == 0 {
                    return Err(Error::Dim {
                        expr: expr.to_string(),
                        detail: format!("division of {x} by zero"),
                    });
                }
                Ok(Some(x / y))
            }
            DimExpr::Range { lo, hi } => {
                let (Some(lo), Some(hi)) = (lo.eval(expr, resolver)?, hi.eval(expr, resolver)?)
                else {
                    return Ok(None);
                };
                hi.checked_sub(lo)
                    .and_then(|span| span.checked_add(1))
                    .map(Some)
                    .ok_or_else(|| Error::Dim {
                        expr: expr.to_string(),
                        detail: format!("range {lo}:{hi} overflows i64"),
                    })
            }
        }
    }
}

/// Apply a checked binary operator to two evaluated operands.
fn eval_binary(
    expr: &str,
    resolver: &dyn DimResolver,
    a: &DimExpr<'_>,
    b: &DimExpr<'_>,
    op: fn(i64, i64) -> Option<i64>,
    symbol: &str,
) -> Result<Option<i64>> {
    let (Some(x), Some(y)) = (a.eval(expr, resolver)?, b.eval(expr, resolver)?) else {
        return Ok(None);
    };
    op(x, y).map(Some).ok_or_else(|| Error::Dim {
        expr: expr.to_string(),
        detail: format!("{x} {symbol} {y} overflows i64"),
    })
}

/// Evaluate every dimension expression of an array declaration.
///
/// Returns the extents in Fortran order (first dimension fastest), or
/// `Ok(None)` as soon as one runtime reference is unset — in which case the
/// caller must skip the shape check. Fails with [`Error::Dim`] on a malformed
/// expression and with [`Error::Malformed`] if it evaluates to a negative
/// extent.
pub fn resolve_dims(dims: &[&str], resolver: &dyn DimResolver) -> Result<Option<Vec<usize>>> {
    let mut extents = Vec::with_capacity(dims.len());
    for &raw in dims {
        let Some(value) = DimExpr::parse(raw)?.eval(raw, resolver)? else {
            return Ok(None);
        };
        let extent = usize::try_from(value).map_err(|_| Error::Dim {
            expr: raw.to_string(),
            detail: format!("extent {value} is negative or does not fit in usize"),
        })?;
        extents.push(extent);
    }
    Ok(Some(extents))
}

/// Recursive-descent parser over the bytes of one expression.
struct Parser<'a> {
    /// Full source text, used for error messages and for slicing identifiers.
    src: &'a str,
    /// `src.as_bytes()`, for byte-wise scanning (the grammar is ASCII-only).
    bytes: &'a [u8],
    /// Byte offset of the next token.
    pos: usize,
}

impl<'a> Parser<'a> {
    fn new(src: &'a str) -> Self {
        Self {
            src,
            bytes: src.as_bytes(),
            pos: 0,
        }
    }

    fn error(&self, detail: &str) -> Error {
        Error::Dim {
            expr: self.src.to_string(),
            detail: detail.to_string(),
        }
    }

    fn skip_whitespace(&mut self) {
        while self.peek().is_some_and(|b| b.is_ascii_whitespace()) {
            self.pos += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn bump(&mut self) {
        self.pos += 1;
    }

    fn parse_sum(&mut self) -> Result<DimExpr<'a>> {
        let mut lhs = self.parse_product()?;
        loop {
            self.skip_whitespace();
            let node = match self.peek() {
                Some(b'+') => {
                    self.bump();
                    DimExpr::Add
                }
                Some(b'-') => {
                    self.bump();
                    DimExpr::Sub
                }
                _ => return Ok(lhs),
            };
            let rhs = self.parse_product()?;
            lhs = node(Box::new(lhs), Box::new(rhs));
        }
    }

    fn parse_product(&mut self) -> Result<DimExpr<'a>> {
        let mut lhs = self.parse_factor()?;
        loop {
            self.skip_whitespace();
            let node = match self.peek() {
                Some(b'*') => {
                    self.bump();
                    DimExpr::Mul
                }
                Some(b'/') => {
                    self.bump();
                    DimExpr::Div
                }
                _ => return Ok(lhs),
            };
            let rhs = self.parse_factor()?;
            lhs = node(Box::new(lhs), Box::new(rhs));
        }
    }

    fn parse_factor(&mut self) -> Result<DimExpr<'a>> {
        self.skip_whitespace();
        match self.peek() {
            Some(b'-') => {
                self.bump();
                Ok(DimExpr::Neg(Box::new(self.parse_factor()?)))
            }
            Some(b'+') => {
                self.bump();
                self.parse_factor()
            }
            Some(b'(') => {
                self.bump();
                let inner = self.parse_sum()?;
                self.skip_whitespace();
                if self.peek() != Some(b')') {
                    return Err(self.error("expected `)`"));
                }
                self.bump();
                Ok(inner)
            }
            Some(b) if b.is_ascii_digit() => self.parse_integer(),
            Some(b) if b.is_ascii_alphabetic() || b == b'_' => self.parse_reference(),
            _ => Err(self.error("expected an integer, a `group.key` reference, `-` or `(`")),
        }
    }

    /// Parse `[0-9]+` as an `i64`.
    fn parse_integer(&mut self) -> Result<DimExpr<'a>> {
        let start = self.pos;
        while self.peek().is_some_and(|b| b.is_ascii_digit()) {
            self.bump();
        }
        let text = &self.src[start..self.pos];
        text.parse::<i64>()
            .map(DimExpr::Int)
            .map_err(|_| self.error("integer literal does not fit in i64"))
    }

    /// Parse `identifier '.' identifier`.
    fn parse_reference(&mut self) -> Result<DimExpr<'a>> {
        let group = self.parse_identifier()?;
        self.skip_whitespace();
        if self.peek() != Some(b'.') {
            return Err(self.error("expected `group.key` (a reference needs a dot)"));
        }
        self.bump();
        self.skip_whitespace();
        let key = self.parse_identifier()?;
        Ok(DimExpr::Ref { group, key })
    }

    /// Consume one identifier; the caller has checked the first byte.
    fn parse_identifier(&mut self) -> Result<&'a str> {
        self.skip_whitespace();
        match self.peek() {
            Some(b) if b.is_ascii_alphabetic() || b == b'_' => {}
            _ => return Err(self.error("expected an identifier")),
        }
        let start = self.pos;
        while self
            .peek()
            .is_some_and(|b| b.is_ascii_alphanumeric() || b == b'_')
        {
            self.bump();
        }
        Ok(&self.src[start..self.pos])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Resolver backed by a fixed map, for the parser tests.
    struct Fixed(HashMap<(&'static str, &'static str), i64>);

    impl Fixed {
        fn new() -> Self {
            let mut values = HashMap::new();
            values.insert(("det", "n_int"), 3);
            values.insert(("det", "bit_kind"), 8);
            values.insert(("det", "n_det"), 1000);
            values.insert(("pseudo", "lmax"), 4);
            Self(values)
        }
    }

    impl DimResolver for Fixed {
        fn resolve(&self, group: &str, key: &str) -> Result<Option<i64>> {
            // The parser borrows the identifiers, so compare through the
            // string contents rather than the &'static str keys.
            Ok(self
                .0
                .iter()
                .find(|((g, k), _)| *g == group && *k == key)
                .map(|(_, v)| *v))
        }
    }

    fn eval(src: &str) -> Result<Option<i64>> {
        let expr = DimExpr::parse(src)?;
        expr.eval(src, &Fixed::new())
    }

    #[test]
    fn parses_literals_and_references() {
        assert_eq!(eval("3").unwrap(), Some(3));
        assert_eq!(eval(" 0 ").unwrap(), Some(0));
        assert_eq!(eval("det.n_det").unwrap(), Some(1000));
    }

    #[test]
    fn parses_arithmetic_with_fortran_precedence() {
        // n_int * bit_kind / 8 = 3 * 8 / 8
        assert_eq!(eval("det.n_int*det.bit_kind/8").unwrap(), Some(3));
        assert_eq!(eval("2 + 3 * 4").unwrap(), Some(14));
        assert_eq!(eval("(2 + 3) * 4").unwrap(), Some(20));
        assert_eq!(eval("-pseudo.lmax").unwrap(), Some(-4));
        assert_eq!(eval("10 - 4 - 3").unwrap(), Some(3));
        assert_eq!(eval("7 / 2").unwrap(), Some(3));
        assert_eq!(eval("-7 / 2").unwrap(), Some(-3));
    }

    #[test]
    fn parses_fortran_ranges() {
        // Extent of 0:lmax is lmax + 1, of -lmax:lmax is 2*lmax + 1.
        assert_eq!(eval("0:pseudo.lmax").unwrap(), Some(5));
        assert_eq!(eval("-pseudo.lmax:pseudo.lmax").unwrap(), Some(9));
        assert_eq!(eval("-3:-1").unwrap(), Some(3));
    }

    #[test]
    fn unresolved_references_yield_none() {
        assert_eq!(eval("pseudo.klocmax").unwrap(), None);
        assert_eq!(eval("2 * pseudo.klocmax").unwrap(), None);
        assert_eq!(eval("0:pseudo.klocmax").unwrap(), None);
    }

    #[test]
    fn rejects_malformed_expressions() {
        for src in [
            "",
            "3 +",
            "3 4",
            "det",
            "det.",
            ".n_det",
            "3:",
            "(3",
            "3)",
            "3 : 4 : 5",
            "3 $ 4",
            "99999999999999999999",
        ] {
            assert!(eval(src).is_err(), "{src:?} should not parse");
        }
    }

    #[test]
    fn division_by_zero_and_overflow_are_errors() {
        assert!(eval("1/0").is_err());
        assert!(eval("9223372036854775807 + 1").is_err());
        assert!(eval("9223372036854775807 * 2").is_err());
    }

    #[test]
    fn resolve_dims_evaluates_in_order() {
        let dims = ["det.n_int*det.bit_kind/8", "2", "det.n_det"];
        assert_eq!(
            resolve_dims(&dims, &Fixed::new()).unwrap(),
            Some(vec![3, 2, 1000])
        );
        assert_eq!(
            resolve_dims(&["2", "pseudo.klocmax"], &Fixed::new()).unwrap(),
            None
        );
        assert!(resolve_dims(&["-1"], &Fixed::new()).is_err());
    }

    #[test]
    fn static_registry_is_sorted_and_complete() {
        let registry = StaticRegistry::new();
        let specs = registry.specs();
        assert_eq!(specs.len(), 332, "one spec per EZFIO.cfg section");
        assert!(
            specs
                .windows(2)
                .all(|pair| (pair[0].group, pair[0].key) < (pair[1].group, pair[1].key)),
            "KEYS must be sorted by (group, key) for binary search"
        );
        // Every dimension reference names a scalar integer attribute.
        for spec in specs {
            let Some(dims) = spec.dims else { continue };
            assert!(
                !dims.is_empty(),
                "{}/{} has empty dims",
                spec.group,
                spec.key
            );
            for raw in dims {
                for name in references(raw) {
                    let (group, key) = name.split_once('.').expect("checked by references()");
                    let target = registry.lookup(group, key).unwrap_or_else(|| {
                        panic!("{}/{}: unknown dim ref {name}", spec.group, spec.key)
                    });
                    assert!(
                        !target.is_array(),
                        "{}/{}: dim ref {name} must be a scalar",
                        spec.group,
                        spec.key
                    );
                    assert!(
                        matches!(target.ty, ElementType::Integer | ElementType::Integer64),
                        "{}/{}: dim ref {name} must be an integer",
                        spec.group,
                        spec.key
                    );
                }
            }
        }
    }

    /// Every `group.key` reference in a dimension expression, in source order.
    fn references(raw: &str) -> Vec<String> {
        let expr = DimExpr::parse(raw).expect("generated expressions must parse");
        let mut refs = Vec::new();
        collect_references(&expr, &mut refs);
        refs
    }

    fn collect_references(expr: &DimExpr<'_>, out: &mut Vec<String>) {
        match expr {
            DimExpr::Ref { group, key } => out.push(format!("{group}.{key}")),
            DimExpr::Int(_) => {}
            DimExpr::Neg(inner) => collect_references(inner, out),
            DimExpr::Add(a, b) | DimExpr::Sub(a, b) | DimExpr::Mul(a, b) | DimExpr::Div(a, b) => {
                collect_references(a, out);
                collect_references(b, out);
            }
            DimExpr::Range { lo, hi } => {
                collect_references(lo, out);
                collect_references(hi, out);
            }
        }
    }

    #[test]
    fn static_registry_knows_the_documented_keys() {
        let registry = StaticRegistry::new();
        let typ = registry.lookup("basis", "typ").expect("basis/typ");
        assert_eq!(typ.ty, ElementType::Character);
        assert_eq!(typ.char_len, 32);
        assert!(!typ.is_array());

        let n_int = registry.lookup("determinants", "n_int").expect("n_int");
        assert_eq!(n_int.ty, ElementType::Integer);
        assert!(!n_int.is_array());

        let psi_det = registry.lookup("determinants", "psi_det").expect("psi_det");
        assert_eq!(psi_det.ty, ElementType::Integer64);
        assert_eq!(
            psi_det.dims,
            Some(
                &[
                    "determinants.n_int*determinants.bit_kind/8",
                    "2",
                    "determinants.n_det"
                ][..]
            )
        );

        let label = registry.lookup("nuclei", "nucl_label").expect("nucl_label");
        assert_eq!(label.ty, ElementType::Character);
        assert_eq!(label.char_len, 32);
        assert_eq!(label.dims, Some(&["nuclei.nucl_num"][..]));

        assert_eq!(registry.lookup("nope", "nope"), None);
        // `basis_type` was renamed by `ezfio_name: typ`.
        assert_eq!(registry.lookup("basis", "basis_type"), None);
    }
}
