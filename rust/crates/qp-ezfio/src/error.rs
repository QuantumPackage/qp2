//! Error type shared by the whole crate.

use std::path::PathBuf;

/// Convenience alias for results returned by this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Every failure mode of the EZFIO reader/writer.
///
/// The messages deliberately echo the reference implementations where they
/// exist (`Attribute <group>/<key> is not set` is the Fortran string), so that
/// the tool layer's user-visible diagnostics stay stable.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// An I/O operation failed.
    #[error("{path}: {source}")]
    Io {
        /// File or directory the operation was attempted on.
        path: PathBuf,
        /// Underlying OS error.
        #[source]
        source: std::io::Error,
    },

    /// A delegated `ezfio_core` operation failed.
    ///
    /// The string is the core's own message (`Error reading data in
    /// <path>.`, `Unable to create directory <path>`, `Unable to write
    /// <path>`). The core folds the underlying `io::Error` and the line
    /// number of a malformed payload into that text, so no structured source
    /// is available here; a caller that needs the OS error gets
    /// [`Error::Io`] from one of the crate's own filesystem probes.
    #[error("{0}")]
    Core(String),

    /// The directory is not an EZFIO directory.
    #[error("{path} is not an EZFIO directory (no .version file)")]
    NotEzFio {
        /// Directory that was opened.
        path: PathBuf,
    },

    /// The requested attribute is not set (no scalar file, no `.gz` array).
    #[error("Attribute {group}/{key} is not set")]
    NotSet {
        /// EZFIO group (directory name).
        group: String,
        /// Attribute name inside the group.
        key: String,
    },

    /// The attribute exists but its on-disk kind is not the requested one.
    #[error("EZFIO attribute {group}/{key}: {detail}")]
    Kind {
        /// EZFIO group (directory name).
        group: String,
        /// Attribute name inside the group.
        key: String,
        /// What went wrong, e.g. "is a gzip array, not a scalar".
        detail: String,
    },

    /// A payload has a shape this crate cannot represent.
    ///
    /// `ezfio_core` reports its own parse failures through [`Error::Core`];
    /// this variant is the crate-level guard on a shape the core returned but
    /// that is outside the format (a negative array extent). The core's
    /// parser rejects negative extents too, so the branch is defensive.
    #[error("{path}: {detail}")]
    Malformed {
        /// File being parsed.
        path: PathBuf,
        /// What is wrong.
        detail: String,
    },

    /// The registry's declared type does not accept the requested Rust type.
    #[error(
        "EZFIO attribute {group}/{key}: registry declares {expected}, call site requested {found}"
    )]
    TypeMismatch {
        /// EZFIO group (directory name).
        group: String,
        /// Attribute name inside the group.
        key: String,
        /// Type declared by the [`crate::Registry`].
        expected: String,
        /// Type implied by the Rust value type.
        found: String,
    },

    /// Array dimensions disagree with the shape declared by the registry.
    ///
    /// A **write-path** error ([`crate::EzFio::set_array`]): reads never
    /// consult the declared extents, so a stale tree is readable.
    #[error(
        "EZFIO attribute {group}/{key}: dimensions {found:?} differ from the declared {expected:?}"
    )]
    Shape {
        /// EZFIO group (directory name).
        group: String,
        /// Attribute name inside the group.
        key: String,
        /// Shape computed from the registry's dimension expressions.
        expected: Vec<usize>,
        /// Shape the caller tried to write.
        found: Vec<usize>,
    },

    /// An [`crate::Array`] was built with dimensions that do not match its data.
    #[error("inconsistent array: dims {dims:?} require {expected} elements, got {found}")]
    InconsistentArray {
        /// Dimensions passed to the constructor.
        dims: Vec<usize>,
        /// Number of elements the dimensions require.
        expected: usize,
        /// Number of elements actually supplied.
        found: usize,
    },

    /// A 1-based multi-index does not address an element of the array.
    #[error("index {indices:?} is out of bounds for dims {dims:?}")]
    IndexOutOfBounds {
        /// Dimensions of the array.
        dims: Vec<usize>,
        /// The offending 1-based indices, or a prefix of them.
        indices: Vec<usize>,
    },

    /// A dimension expression is syntactically invalid or cannot be resolved.
    #[error("dimension expression {expr:?}: {detail}")]
    Dim {
        /// The expression as written in the EZFIO config (`size:` field).
        expr: String,
        /// Why it could not be parsed or evaluated.
        detail: String,
    },

    /// An array rank cannot be stored in the EZFIO header.
    ///
    /// The header writes the rank with `I3`, but [`crate::MAX_RANK`] (10) is
    /// the rank `ezfio_core` accepts: its Fortran wrapper hands the caller a
    /// fixed `dims(10)` output buffer, so a deeper array could be written but
    /// never read back through this layer.
    #[error("EZFIO array rank {rank} is not supported (the limit is 1..=10)")]
    Rank {
        /// The offending rank (number of dimensions).
        rank: usize,
    },

    /// An array extent is too large for the `I20` dimension field.
    #[error("dimension {index} = {value} does not fit in a Fortran I20 field")]
    DimensionTooLarge {
        /// 0-based position of the extent in the dimension list.
        index: usize,
        /// The extent that does not fit in an `i64`.
        value: usize,
    },

    /// A group or attribute name is not usable as a path component.
    #[error("invalid EZFIO {kind} name {name:?}: {detail}")]
    InvalidName {
        /// `"group"` or `"attribute"`.
        kind: &'static str,
        /// The offending name.
        name: String,
        /// Why it was rejected.
        detail: &'static str,
    },
}

impl Error {
    /// Wrap an [`std::io::Error`] with the path it happened on.
    ///
    /// `path` is the file or directory the failed call was made on; it is
    /// copied into the error, so the caller keeps ownership.
    pub(crate) fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        Self::Io {
            path: path.into(),
            source,
        }
    }

    /// Build a [`Error::Malformed`] with a formatted detail string.
    pub(crate) fn malformed(path: impl Into<PathBuf>, detail: impl Into<String>) -> Self {
        Self::Malformed {
            path: path.into(),
            detail: detail.into(),
        }
    }
}

/// Attach the path an I/O call was made on to its error.
pub(crate) trait IoContext<T> {
    /// Map the error of `self` into [`Error::Io`] with `path` as context.
    fn at(self, path: &std::path::Path) -> Result<T>;
}

impl<T> IoContext<T> for std::result::Result<T, std::io::Error> {
    fn at(self, path: &std::path::Path) -> Result<T> {
        self.map_err(|e| Error::io(path, e))
    }
}
