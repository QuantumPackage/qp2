//! Minimal S-expression reader for the RST parse-back path.
//!
//! `qp_edit` never sends S-expressions to EZFIO.  It rebuilds one from the
//! RST text it has just rendered (or that the user edited), because the
//! OCaml tool parses every section with
//! `Sexp.of_string_conv` on the sexp derived from the `key = value` lines.
//! This module is that reader plus the record-shape checks of
//! `ppx_sexp_conv`, no more: it does not print S-expressions.
//!
//! # What is reproduced
//!
//! * Atoms, quoted atoms (`"..."` with the OCaml escapes), lists, and `;`
//!   comments to end of line.
//! * Byte spans for every node.  The OCaml `fail_msg` reports the
//!   offending sub-expression by slicing the source string, so the spans
//!   must refer to the string that was parsed (the `"(" ^ s ^ ")"` wrapper
//!   in the caller), not to an unescaped copy.
//! * `ppx_sexp_conv`'s record rules: a list of two-element lists, known
//!   names only, no duplicates, every field present, `()` rejected.
//!
//! # Deliberate deviations
//!
//! * Error *messages* are not sexplib's.  Sexplib's converter exceptions
//!   are phrased by the generated code (`t_of_sexp`, `int_of_sexp`, ...) and
//!   are not part of the byte contract (the C1–C7 gate never triggers
//!   them; only `-r` on malformed text does).  The *snippet* that follows
//!   the message is reproduced by [`snippet`] exactly.
//! * `\x` escapes are accepted on input; sexplib only produces them.
//!
//! Units: none.

/// A parse or conversion failure: a message and the byte span it refers to
/// in the parsed string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SexpError {
    /// Human-readable description, printed after `Error:  (`.
    pub message: String,
    /// Start byte offset of the offending sub-expression (inclusive).
    pub start: usize,
    /// End byte offset (exclusive).
    pub end: usize,
}

impl SexpError {
    /// Build an error spanning `node`.
    fn at(message: impl Into<String>, node: &Sexp) -> SexpError {
        let (start, end) = node.span();
        SexpError {
            message: message.into(),
            start,
            end,
        }
    }
}

/// One parsed node, with its byte span in the source string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Sexp {
    /// A bare or quoted atom; `value` is the unescaped text.
    Atom {
        /// Unescaped atom text.
        value: String,
        /// Start offset of the token, including its opening quote if any.
        start: usize,
        /// End offset, one past the closing quote.
        end: usize,
    },
    /// A parenthesised list.
    List {
        /// Items in source order.
        items: Vec<Sexp>,
        /// Offset of the opening parenthesis.
        start: usize,
        /// Offset one past the closing parenthesis.
        end: usize,
    },
}

impl Sexp {
    /// Byte span of the node in the parsed string.
    pub(crate) fn span(&self) -> (usize, usize) {
        match self {
            Sexp::Atom { start, end, .. } | Sexp::List { start, end, .. } => (*start, *end),
        }
    }

    /// Atom text, or `None` for a list.
    pub(crate) fn atom(&self) -> Option<&str> {
        match self {
            Sexp::Atom { value, .. } => Some(value),
            Sexp::List { .. } => None,
        }
    }
}

/// Parse one S-expression; trailing whitespace and comments are allowed,
/// anything else after the first node is an error.
///
/// # Errors
///
/// Returns a [`SexpError`] for an unbalanced parenthesis, an unterminated
/// quoted atom, or a second top-level node.
pub(crate) fn parse(input: &str) -> Result<Sexp, SexpError> {
    let mut parser = Parser {
        bytes: input.as_bytes(),
        input,
        pos: 0,
    };
    parser.skip_trivia();
    let node = parser.node()?;
    parser.skip_trivia();
    if parser.pos != parser.bytes.len() {
        return Err(SexpError {
            message: "Sexp.of_string: trailing content after the expression".to_string(),
            start: parser.pos,
            end: parser.bytes.len(),
        });
    }
    Ok(node)
}

struct Parser<'a> {
    bytes: &'a [u8],
    input: &'a str,
    pos: usize,
}

impl Parser<'_> {
    /// Skip whitespace and `;` comments.
    fn skip_trivia(&mut self) {
        loop {
            match self.bytes.get(self.pos) {
                Some(b' ' | b'\t' | b'\n' | b'\r' | 0x0c) => self.pos += 1,
                Some(b';') => {
                    while let Some(byte) = self.bytes.get(self.pos) {
                        if *byte == b'\n' {
                            break;
                        }
                        self.pos += 1;
                    }
                }
                _ => return,
            }
        }
    }

    /// Parse one node at the current position (trivia already skipped).
    fn node(&mut self) -> Result<Sexp, SexpError> {
        match self.bytes.get(self.pos) {
            Some(b'(') => self.list(),
            Some(b')') => Err(SexpError {
                message: "Sexp.of_string: unexpected ')'".to_string(),
                start: self.pos,
                end: self.pos + 1,
            }),
            Some(b'"') => self.quoted_atom(),
            Some(_) => self.bare_atom(),
            None => Err(SexpError {
                message: "Sexp.of_string: unexpected end of input".to_string(),
                start: self.pos,
                end: self.pos,
            }),
        }
    }

    fn list(&mut self) -> Result<Sexp, SexpError> {
        let start = self.pos;
        self.pos += 1; // '('
        let mut items = Vec::new();
        loop {
            self.skip_trivia();
            match self.bytes.get(self.pos) {
                Some(b')') if items.is_empty() => {
                    // sexplib lexes `()` — and `(   )` — as the unit atom,
                    // not as an empty list; the record reader then reports
                    // "expected a list", exactly like OCaml.
                    self.pos += 1;
                    return Ok(Sexp::Atom {
                        value: "()".to_string(),
                        start,
                        end: self.pos,
                    });
                }
                Some(b')') => {
                    self.pos += 1;
                    return Ok(Sexp::List {
                        items,
                        start,
                        end: self.pos,
                    });
                }
                None => {
                    return Err(SexpError {
                        message: "Sexp.of_string: unterminated list".to_string(),
                        start,
                        end: self.bytes.len(),
                    });
                }
                Some(_) => items.push(self.node()?),
            }
        }
    }

    /// A quoted atom.  Escapes follow OCaml's lexical conventions because
    /// sexplib uses `Scanf`-style escapes: `\\`, `\"`, `\n`, `\t`, `\r`,
    /// `\b`, `\xHH`, and `\DDD` (decimal character code).
    fn quoted_atom(&mut self) -> Result<Sexp, SexpError> {
        let start = self.pos;
        self.pos += 1; // opening quote
        let mut value = String::new();
        loop {
            let Some(&byte) = self.bytes.get(self.pos) else {
                return Err(SexpError {
                    message: "Sexp.of_string: unterminated string".to_string(),
                    start,
                    end: self.bytes.len(),
                });
            };
            match byte {
                b'"' => {
                    self.pos += 1;
                    return Ok(Sexp::Atom {
                        value,
                        start,
                        end: self.pos,
                    });
                }
                b'\\' => {
                    if let Some(escaped) = self.escape(start)? {
                        value.push(escaped);
                    }
                }
                _ => {
                    // Copy one UTF-8 scalar, not one byte.
                    let ch = self.input[self.pos..]
                        .chars()
                        .next()
                        .expect("pos < len is a char boundary here");
                    value.push(ch);
                    self.pos += ch.len_utf8();
                }
            }
        }
    }

    /// Decode one backslash escape; `self.pos` is on the backslash.
    ///
    /// Returns `None` for a backslash-newline line continuation, which
    /// contributes no character.
    fn escape(&mut self, atom_start: usize) -> Result<Option<char>, SexpError> {
        self.pos += 1;
        let Some(&byte) = self.bytes.get(self.pos) else {
            return Err(SexpError {
                message: "Sexp.of_string: unterminated escape".to_string(),
                start: atom_start,
                end: self.bytes.len(),
            });
        };
        let simple = match byte {
            b'\\' => Some('\\'),
            b'"' => Some('"'),
            b'\'' => Some('\''),
            b'n' => Some('\n'),
            b't' => Some('\t'),
            b'r' => Some('\r'),
            b'b' => Some('\u{8}'),
            b' ' | b'\n' => {
                // Backslash-newline: skip the following whitespace.
                self.pos += 1;
                while matches!(self.bytes.get(self.pos), Some(b' ' | b'\t' | b'\n' | b'\r')) {
                    self.pos += 1;
                }
                return Ok(None);
            }
            _ => None,
        };
        if let Some(ch) = simple {
            self.pos += 1;
            return Ok(Some(ch));
        }
        if byte == b'x' {
            // \xHH: one or two hex digits.
            self.pos += 1;
            let mut code = 0u32;
            let mut digits = 0;
            while digits < 2 {
                match self
                    .bytes
                    .get(self.pos)
                    .and_then(|b| (*b as char).to_digit(16))
                {
                    Some(digit) => {
                        code = code * 16 + digit;
                        self.pos += 1;
                        digits += 1;
                    }
                    None => break,
                }
            }
            if digits == 0 {
                return Err(SexpError {
                    message: "Sexp.of_string: bad \\x escape".to_string(),
                    start: atom_start,
                    end: self.pos,
                });
            }
            return Ok(Some(char::from_u32(code & 0xff).unwrap_or('\u{0}')));
        }
        if byte.is_ascii_digit() {
            // \DDD: up to three decimal digits.
            let mut code = 0u32;
            let mut digits = 0;
            while digits < 3 {
                match self.bytes.get(self.pos) {
                    Some(d) if d.is_ascii_digit() => {
                        code = code * 10 + u32::from(d - b'0');
                        self.pos += 1;
                        digits += 1;
                    }
                    _ => break,
                }
            }
            return Ok(Some(char::from_u32(code).unwrap_or('\u{0}')));
        }
        // Unknown escape: keep the escaped byte, as OCaml's lexer does for
        // a character that has no special meaning.
        let ch = self.input[self.pos..]
            .chars()
            .next()
            .expect("byte exists, so a char exists");
        self.pos += ch.len_utf8();
        Ok(Some(ch))
    }

    /// A bare atom: everything up to whitespace, `(`, `)`, `"` or `;`.
    fn bare_atom(&mut self) -> Result<Sexp, SexpError> {
        let start = self.pos;
        while let Some(&byte) = self.bytes.get(self.pos) {
            if matches!(
                byte,
                b' ' | b'\t' | b'\n' | b'\r' | 0x0c | b'(' | b')' | b'"' | b';'
            ) {
                break;
            }
            self.pos += 1;
        }
        Ok(Sexp::Atom {
            value: self.input[start..self.pos].to_string(),
            start,
            end: self.pos,
        })
    }
}

/// One accepted record field: `(name, pair, value)`.
///
/// `pair` is the whole `(name value)` node, whose span is what a value
/// conversion error reports in `fail_msg`; `value` is the value node the
/// converter reads.
pub(crate) type RecordField<'a> = (&'a str, &'a Sexp, &'a Sexp);

/// Check a parsed node against a record schema and return its fields.
///
/// `expected` lists the accepted field names.  The returned fields are in
/// **input order**; callers that must apply values in record order index
/// them by name.
///
/// This reproduces `ppx_sexp_conv`'s record converter, which is what
/// `[@@deriving sexp]` puts behind `t_of_sexp` in the generated
/// `Input_*.ml` files: the node must be a list of `(name value)` pairs,
/// every name must be a declared field, no name may repeat, and every
/// declared field must appear.
///
/// # Errors
///
/// The error span is the offending pair for shape, unknown-name and
/// duplicate errors, and the whole record for a missing field, matching
/// where `Sexp.Annotated.conv` reports the failure.
pub(crate) fn record_fields<'a>(
    node: &'a Sexp,
    expected: &[&str],
) -> Result<Vec<RecordField<'a>>, SexpError> {
    let Sexp::List { items, .. } = node else {
        return Err(SexpError::at(
            "record_of_sexp: expected a list of (field value) pairs",
            node,
        ));
    };
    let mut pairs: Vec<RecordField> = Vec::with_capacity(items.len());
    for item in items {
        let Sexp::List { items: pair, .. } = item else {
            return Err(SexpError::at(
                "record_of_sexp: expected a (field value) pair",
                item,
            ));
        };
        if pair.len() != 2 {
            return Err(SexpError::at(
                "record_of_sexp: expected a (field value) pair",
                item,
            ));
        }
        let Some(name) = pair[0].atom() else {
            return Err(SexpError::at(
                "record_of_sexp: the field name must be an atom",
                item,
            ));
        };
        if !expected.contains(&name) {
            return Err(SexpError::at(
                format!("record_of_sexp: unknown field name: {name}"),
                item,
            ));
        }
        if pairs.iter().any(|(seen, _, _)| *seen == name) {
            return Err(SexpError::at(
                format!("record_of_sexp: duplicate field name: {name}"),
                item,
            ));
        }
        pairs.push((name, item, &pair[1]));
    }
    for name in expected {
        if !pairs.iter().any(|(seen, _, _)| seen == name) {
            return Err(SexpError::at(
                format!("record_of_sexp: missing field: {name}"),
                node,
            ));
        }
    }
    Ok(pairs)
}

/// Reproduce `Generic_input_of_rst.fail_msg`: print the conversion error and
/// its `##`-marked excerpt to stderr.
///
/// The format string is the OCaml one, `"Error:  (%s)\n\n  %s\n\n"`;
/// `source` is the string handed to [`parse`], so the snippet's byte spans
/// are valid.
pub(crate) fn fail_msg(source: &str, error: &SexpError) {
    eprint!(
        "Error:  ({})\n\n  {}\n\n",
        error.message,
        snippet(source, error.start, error.end)
    );
}

/// Reproduce `Generic_input_of_rst.fail_msg`'s snippet.
///
/// `source` is the exact string handed to [`parse`] (including the
/// enclosing parentheses); `start..end` is the failing sub-expression.
/// OCaml marks it with `##`, replaces every `(` with a space, splits on
/// `)`, trims each piece, keeps the pieces that contain `##`, and joins
/// them with newlines.  The result is the text printed after
/// `Error:  (<message>)` and two spaces of indentation.
pub(crate) fn snippet(source: &str, start: usize, end: usize) -> String {
    let marked = format!(
        "{} ## {} ## {}",
        &source[..start],
        &source[start..end],
        &source[end..]
    );
    let replaced = marked.replace('(', " ");
    replaced
        .split(')')
        .map(str::trim)
        .filter(|piece| piece.contains("##"))
        .collect::<Vec<_>>()
        .join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atoms_and_lists_keep_spans() {
        let source = "((     key   value)(other 1))";
        let node = parse(source).expect("valid sexp");
        let Sexp::List { items, start, end } = &node else {
            panic!("top level must be a list");
        };
        assert_eq!((*start, *end), (0, source.len()));
        assert_eq!(items.len(), 2);
        assert_eq!(items[0].atom(), None);
        let Sexp::List { items: pair, .. } = &items[0] else {
            panic!("item must be a list");
        };
        assert_eq!(pair[0].atom(), Some("key"));
        assert_eq!(pair[1].atom(), Some("value"));
    }

    #[test]
    fn unit_atom_is_rejected_by_the_record_reader() {
        // sexplib lexes `()` as the atom "()", not as an empty list.
        let node = parse("()").expect("valid sexp");
        assert_eq!(node.atom(), Some("()"));
        let error = record_fields(&node, &["a"]).expect_err("not a record");
        assert!(error.message.contains("list of (field value) pairs"));
        assert_eq!(node.span(), (0, 2));
        assert_eq!(parse("(   )").expect("valid sexp").atom(), Some("()"));
    }

    #[test]
    fn record_accepts_all_fields_in_any_order() {
        let node = parse("((b 2)(a 1))").expect("valid sexp");
        let pairs = record_fields(&node, &["a", "b"]).expect("complete record");
        assert_eq!(
            pairs.iter().map(|(name, _, _)| *name).collect::<Vec<_>>(),
            ["b", "a"]
        );
        // The pair node spans `(a 1)`, the value node only `1`.
        assert_eq!(pairs[1].1.span(), (6, 11));
        assert_eq!(pairs[1].2.span(), (9, 10));
    }

    #[test]
    fn record_rejects_missing_extra_duplicate_and_non_pair() {
        let node = parse("((a 1))").expect("valid sexp");
        let err = record_fields(&node, &["a", "b"]).expect_err("missing b");
        assert!(err.message.contains("missing field: b"));
        let node = parse("((a 1)(c 2))").expect("valid sexp");
        let err = record_fields(&node, &["a"]).expect_err("unknown c");
        assert!(err.message.contains("unknown field name: c"));
        let node = parse("((a 1)(a 2))").expect("valid sexp");
        let err = record_fields(&node, &["a"]).expect_err("duplicate a");
        assert!(err.message.contains("duplicate field name: a"));
        let node = parse("((a 1 2))").expect("valid sexp");
        assert!(record_fields(&node, &["a"]).is_err());
        let node = parse("(a)").expect("valid sexp");
        assert!(record_fields(&node, &["a"]).is_err());
    }

    #[test]
    fn quoted_atoms_unescape_ocaml_escapes() {
        let node = parse(r#""a\"b\\c\n""#).expect("valid sexp");
        assert_eq!(node.atom(), Some("a\"b\\c\n"));
    }

    #[test]
    fn quoted_atom_span_includes_the_quotes() {
        let source = "(x \"y\")";
        let node = parse(source).expect("valid sexp");
        let Sexp::List { items, .. } = &node else {
            panic!("list");
        };
        assert_eq!(items[1].span(), (3, 6));
    }

    #[test]
    fn comments_and_whitespace_are_trivia() {
        let node = parse(" ; hello\n (a 1) ; bye\n").expect("valid sexp");
        assert!(matches!(node, Sexp::List { .. }));
    }

    #[test]
    fn unbalanced_input_is_an_error() {
        assert!(parse("((a 1)").is_err());
        assert!(parse("(a 1))").is_err());
        assert!(parse("").is_err());
        assert!(parse("(a) (b)").is_err());
    }

    /// Worked by hand from `Generic_input_of_rst.fail_msg` for the source
    /// `((a 1)(b 2))` and the range of `(b 2)`.
    #[test]
    fn snippet_marks_the_offending_pair() {
        let source = "((a 1)(b 2))";
        assert_eq!(snippet(source, 6, 11), "##  b 2\n##");
    }

    /// Whole-record range: `"("`→`" "` turns the two opening parentheses
    /// into two spaces, which follow the single space of the leading
    /// `" ## "` marker.
    #[test]
    fn snippet_of_the_whole_record_is_one_marked_line() {
        let source = "((a 1)(b 2))";
        assert_eq!(snippet(source, 0, source.len()), "##   a 1\n##");
    }
}
