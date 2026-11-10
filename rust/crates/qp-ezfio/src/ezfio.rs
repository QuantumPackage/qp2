//! The EZFIO directory handle: groups, scalars, gzip arrays.
//!
//! This module is the QP2 tool layer of the EZFIO format. **No byte is
//! formatted or parsed here**: every payload is read and written by
//! [`ezfio_core::ops`], the implementation of the on-disk format that the
//! Fortran, Python and OCaml layers of this build also use. What remains is
//! the handle ([`EzFio`]), the registry checks, and the directory-level
//! operations the core does not have (listing, removal).
//!
//! Two invariants drive the code:
//!
//! * **Presence is the filesystem.** An attribute exists iff its payload file
//!   exists; there is no index to keep in sync. Writes therefore create the
//!   group directory on demand, and `ezfio_core` installs them with a rename,
//!   so that a concurrent reader sees either the old or the new payload,
//!   never a partial one.
//! * **Reading is dynamic.** Neither reference client consults a schema before
//!   reading, so a key that the registry does not know (a plugin's attribute)
//!   must still be readable and writable. [`crate::Registry`] checks are
//!   applied only when the key is known, they cover the element **type** on
//!   both paths, and the declared **shape** is checked on writes only: a tree
//!   whose scalars moved after the arrays were written stays readable, as it
//!   does under the reference getters.
//!
//! # Writer dialect
//!
//! `ezfio_core` exposes two byte dialects for its writers. This crate always
//! passes [`WriterDialect::Ocaml`] (the private `WRITE_DIALECT`): it replaces
//! the reference OCaml tool layer, and the golden harness compares EZFIO state
//! reference-vs-Rust byte for byte. The dialect never leaks through [`EzFio`]:
//! there is one writer, and it is not configurable.
//!
//! # Costs
//!
//! * [`EzFio::get_array`] decompresses the payload **twice**: once in
//!   `ops::array_info` for the shape (the core exposes shape and data as two
//!   calls, and its `array_info` decompresses the whole stream to parse the
//!   two header lines) and once in `ops::read_array` for the elements. Peak
//!   memory is one array plus the decompressed text buffer of one call,
//!   because the text is dropped before the second call starts. A concurrent
//!   write between the two calls makes the shapes disagree and the read fails
//!   with the core's `BadDims` (no tool writes an attribute it is reading).
//! * [`EzFio::set_array`] copies the elements once into the owned
//!   [`ezfio_core::ArrayData`] the core takes, and the core materialises the
//!   whole decompressed payload as a `String` before compressing it, so peak
//!   memory is about twice the array plus the text image (roughly 25 bytes
//!   per real element, 21 per integer).
//! * [`EzFio::get`]/[`EzFio::get_array`]/[`DimResolver::resolve`] open the
//!   payload once to surface an OS-level failure as [`Error::Io`] (the core
//!   folds "cannot open" and "cannot parse" into one error and drops the
//!   `io::Error`); the core re-opens it.

use std::fs::{self, File};
use std::io::ErrorKind;
use std::path::{Path, PathBuf};

use ezfio_core::ops;
use ezfio_core::{EzfioError, WriterDialect};

use crate::array::Array;
use crate::error::{Error, IoContext, Result};
use crate::layout::{Element, ElementType};
use crate::registry::{DimResolver, KeySpec, Registry, StaticRegistry, resolve_dims};

/// Version written into every `.version` file.
///
/// Re-exported from `ezfio_core`, so the marker this crate writes and the one
/// the rest of the build writes cannot drift. Readers do not enforce it; the
/// file is only a marker that a directory belongs to an EZFIO tree.
pub const EZFIO_VERSION: &str = ezfio_core::VERSION;

/// Maximum number of array dimensions.
///
/// This is `ezfio_core`'s limit (the Fortran wrapper's `dims(10)`): a payload
/// with a larger rank cannot be returned by `ops::array_info`, so a writer
/// that produced one would emit an array this crate cannot read back. The
/// reference readers accept it; the core does not, and this crate matches the
/// core it links.
pub const MAX_RANK: usize = 10;

/// Name of the marker file at the root of the tree and in every group.
const VERSION_FILE: &str = ".version";

/// Byte dialect of every write.
///
/// [`WriterDialect::Ocaml`] reproduces the reference client this crate
/// replaces: lower-case `e` exponents, characters written verbatim (blanks
/// are data), and the doubled newline of a logical scalar (`T\n\n`), which
/// every OCaml-client-written tree of the golden corpus carries. The dialect
/// is an implementation detail of this module, not a property of the handle.
const WRITE_DIALECT: WriterDialect = WriterDialect::Ocaml;

/// A handle on one EZFIO directory tree.
///
/// The handle is lightweight: it holds the root path and a registry. All
/// state lives in the filesystem, so two handles opened on the same directory
/// see each other's writes immediately.
#[derive(Debug)]
pub struct EzFio {
    /// Root directory, i.e. the one holding `.version`.
    root: PathBuf,
    /// Declared attributes, consulted by the typed accessors when the key is
    /// known. Boxed so that a code-generated registry can replace the static
    /// one without making [`EzFio`] generic.
    registry: Box<dyn Registry>,
}

impl EzFio {
    /// True when `path` is the root of an EZFIO tree.
    ///
    /// That is exactly the presence of `<path>/.version`, which is also what
    /// the reference clients test (`ezfio_core::ops::group_exists`). `false`
    /// is returned for a missing path or a directory without the marker,
    /// including on permission errors (the predicate is a cheap probe, not a
    /// diagnosis).
    pub fn is_ezfio_dir(path: &Path) -> bool {
        ops::group_exists(path)
    }

    /// Open an existing EZFIO tree.
    ///
    /// - `root`: directory containing `.version`; symlinks are followed.
    ///
    /// Fails with [`Error::NotEzFio`] when the marker is absent; the directory
    /// is not scanned or validated beyond that.
    pub fn open(root: impl Into<PathBuf>) -> Result<Self> {
        let root = root.into();
        if !Self::is_ezfio_dir(&root) {
            return Err(Error::NotEzFio { path: root });
        }
        Ok(Self {
            root,
            registry: Box::new(StaticRegistry::new()),
        })
    }

    /// Create an EZFIO tree, or open it if it already exists.
    ///
    /// - `root`: directory to create (parents are created too).
    ///
    /// Delegates to `ezfio_core::ops::mkdir`, which creates the directory
    /// (`create_dir_all`) and stamps it with `.version` unless the stamp is
    /// already there. The `ezfio` group stamps (`ezfio/creation`,
    /// `ezfio/user`, `ezfio/library`) are *attributes*, not format metadata:
    /// both reference clients set them in their creation path, so a caller
    /// that wants a byte-identical directory writes them with [`EzFio::set`]
    /// (`library` is the caller's basis library, which this crate cannot
    /// know).
    ///
    /// On failure the error is [`Error::Core`]: the core reports
    /// `Unable to create directory <path>` and, unlike the previous local
    /// implementation, does not carry the underlying `io::Error`.
    pub fn create(root: impl Into<PathBuf>) -> Result<Self> {
        let root = root.into();
        ops::mkdir(&root).map_err(core_error)?;
        Ok(Self {
            root,
            registry: Box::new(StaticRegistry::new()),
        })
    }

    /// Replace the registry used for typed checks (builder style).
    ///
    /// - `registry`: table consulted by [`EzFio::get`], [`EzFio::set`],
    ///   [`EzFio::get_array`] and [`EzFio::set_array`]; [`StaticRegistry`] is
    ///   the default.
    pub fn with_registry(mut self, registry: Box<dyn Registry>) -> Self {
        self.registry = registry;
        self
    }

    /// The root directory this handle was opened on.
    pub fn root(&self) -> &Path {
        &self.root
    }

    /// The version string stored in the root `.version` file, trimmed.
    ///
    /// Reference trees contain `2.0.8`. This is informational only; no
    /// compatibility check is performed when opening.
    pub fn version(&self) -> Result<String> {
        let path = self.root.join(VERSION_FILE);
        let text = fs::read_to_string(&path).at(&path)?;
        Ok(text.lines().next().unwrap_or("").trim().to_string())
    }

    /// List the groups of the tree, sorted by name.
    ///
    /// A group is a subdirectory carrying its own `.version` marker, which is
    /// what both reference `mkdir` paths write. Non-directories and entries
    /// whose names are not valid UTF-8 are skipped.
    pub fn groups(&self) -> Result<Vec<String>> {
        let entries = match fs::read_dir(&self.root) {
            Ok(entries) => entries,
            Err(e) if e.kind() == ErrorKind::NotFound => return Ok(Vec::new()),
            Err(e) => return Err(Error::io(&self.root, e)),
        };
        let mut groups = Vec::new();
        for entry in entries {
            let entry = entry.map_err(|e| Error::io(&self.root, e))?;
            let path = entry.path();
            if path.is_dir()
                && let Some(name) = entry.file_name().to_str()
                && Self::is_ezfio_dir(&path)
            {
                groups.push(name.to_string());
            }
        }
        groups.sort();
        Ok(groups)
    }

    /// List the attributes of `group`, sorted by name and without duplicates.
    ///
    /// Both scalars (`<key>`) and arrays (`<key>.gz`) are listed under `key`;
    /// an attribute present in both forms appears once. Files whose name starts
    /// with `.` (the `.version` marker, the core's leftover temporary files)
    /// are skipped. A missing group directory yields an empty list, like the
    /// reference `Sys.file_exists`-based probes.
    pub fn keys(&self, group: &str) -> Result<Vec<String>> {
        Self::check_name("group", group)?;
        let dir = self.root.join(group);
        let entries = match fs::read_dir(&dir) {
            Ok(entries) => entries,
            Err(e) if e.kind() == ErrorKind::NotFound => return Ok(Vec::new()),
            Err(e) => return Err(Error::io(&dir, e)),
        };
        let mut keys = Vec::new();
        for entry in entries {
            let entry = entry.map_err(|e| Error::io(&dir, e))?;
            if !entry.path().is_file() {
                continue;
            }
            let Some(name) = entry.file_name().to_str().map(str::to_string) else {
                continue;
            };
            if name.starts_with('.') {
                continue;
            }
            let key = name.strip_suffix(".gz").unwrap_or(&name).to_string();
            keys.push(key);
        }
        keys.sort();
        keys.dedup();
        Ok(keys)
    }

    /// True when `group/key` is set.
    ///
    /// Mirrors the reference client: the group directory must carry its
    /// `.version` marker, and the payload file must exist. When the registry
    /// declares the key, only the declared form is probed (`<key>` for a
    /// scalar, `<key>.gz` for an array, `ezfio_core::ops::attr_exists`); a
    /// registry-unknown key, like a plugin attribute, is present in either
    /// form. Errors are reserved for invalid names; filesystem probes treat
    /// "not accessible" as "absent", exactly like `Sys.file_exists`.
    pub fn has(&self, group: &str, key: &str) -> Result<bool> {
        Self::check_name("group", group)?;
        Self::check_name("attribute", key)?;
        let dir = self.root.join(group);
        if !Self::is_ezfio_dir(&dir) {
            return Ok(false);
        }
        match self.registry.lookup(group, key) {
            Some(spec) => Ok(ops::attr_exists(&dir, key, spec.is_array())),
            None => Ok(dir.join(key).is_file() || dir.join(format!("{key}.gz")).is_file()),
        }
    }

    /// Remove `group/key`, in scalar and/or array form.
    ///
    /// Returns `true` if at least one file was removed. A missing attribute is
    /// not an error, matching the reference clients' `rm -f` / no-op delete;
    /// the group directory itself is kept. `ezfio_core` has no removal
    /// operation (the C ABI does not expose one), so the two `remove_file`
    /// calls stay here. The `.version` markers cannot be removed through this
    /// method because names starting with `.` are rejected.
    pub fn remove(&self, group: &str, key: &str) -> Result<bool> {
        Self::check_name("group", group)?;
        Self::check_name("attribute", key)?;
        let dir = self.root.join(group);
        let mut removed = false;
        for path in [dir.join(key), dir.join(format!("{key}.gz"))] {
            match fs::remove_file(&path) {
                Ok(()) => removed = true,
                Err(e) if e.kind() == ErrorKind::NotFound => {}
                Err(e) => return Err(Error::io(&path, e)),
            }
        }
        Ok(removed)
    }

    /// Create `group` if needed, with its `.version` marker.
    ///
    /// Called automatically by [`EzFio::set`] and [`EzFio::set_array`], like
    /// the reference OCaml `create_group`; a group directory that exists
    /// without a marker gets one, like the Fortran `ezfio_mkdir`. Delegates to
    /// `ezfio_core::ops::mkdir`: idempotent, and a no-op (`stat` only) when
    /// the marker is already present.
    pub fn create_group(&self, group: &str) -> Result<()> {
        Self::check_name("group", group)?;
        ops::mkdir(&self.root.join(group)).map_err(core_error)
    }

    /// Read the scalar `group/key` as `T`.
    ///
    /// - `group`, `key`: EZFIO names; validated as single path components.
    ///
    /// The file is read dynamically: its on-disk kind wins over the registry,
    /// so a plugin key or a scalar written by either reference client is read
    /// the same way. `ezfio_core` reads the first line and trims it with the
    /// reference clients' whitespace set, so padded Fortran payloads and OCaml
    /// payloads are equivalent.
    ///
    /// Fails with [`Error::NotSet`] when neither a scalar nor an array exists,
    /// with [`Error::Kind`] when only the array form exists, with
    /// [`Error::TypeMismatch`] when the registry declares another element type,
    /// with [`Error::Io`] when the file exists but cannot be opened, and with
    /// [`Error::Core`] when the value does not parse (the core reports that as
    /// `Error reading data in <path>.`, without a line number).
    pub fn get<T: Element>(&self, group: &str, key: &str) -> Result<T> {
        Self::check_name("group", group)?;
        Self::check_name("attribute", key)?;
        self.check_scalar_spec::<T>(group, key)?;
        let dir = self.root.join(group);
        let path = dir.join(key);
        if !path.is_file() {
            return Err(self.not_set_or_kind(group, key, false));
        }
        probe_readable(&path)?;
        let value = ops::read_scalar(&dir, key, T::TYPE.to_core_ty()).map_err(core_error)?;
        T::from_scalar(value).ok_or_else(|| core_value_mismatch::<T>(&path))
    }

    /// Write the scalar `group/key` with value `value`.
    ///
    /// - `group`, `key`: EZFIO names; the group is created on demand.
    /// - `value`: handed to `ezfio_core`, which formats it with the OCaml
    ///   dialect (see the module note on the writer dialect) and installs
    ///   `<field>\n` atomically.
    ///
    /// Character scalars are written verbatim: the reference OCaml client
    /// does not pad a scalar, and padding carries no information because every
    /// reader trims.
    ///
    /// Fails with [`Error::Kind`] when the registry declares an array, with
    /// [`Error::TypeMismatch`] when it declares another element type, and with
    /// [`Error::Core`] when the core cannot create the group or write the
    /// file (`Unable to create directory <path>` / `Unable to write <path>`).
    pub fn set<T: Element>(&self, group: &str, key: &str, value: &T) -> Result<()> {
        Self::check_name("group", group)?;
        Self::check_name("attribute", key)?;
        self.check_scalar_spec::<T>(group, key)?;
        self.create_group(group)?;
        let dir = self.root.join(group);
        let scalar = value.to_scalar();
        ops::write_scalar(&dir, key, T::TYPE.to_core_ty(), &scalar, WRITE_DIALECT)
            .map_err(core_error)
    }

    /// Read the gzip array `group/key` as `T`.
    ///
    /// The dimensions are read from the file, **never from the registry**: the
    /// reference getters (`read_float_array`, `read_int64_array`, …) do not
    /// consult the declared `size:` either, so a stale or truncated tree stays
    /// readable, and a section reader that needs a shape carries its own check
    /// or repair (e.g. the `determinants` record reader and its
    /// `state_average_weight` re-materialisation).
    ///
    /// The shape comes from `ops::array_info` and the elements from
    /// `ops::read_array`, so the payload is decompressed twice (see the module
    /// docs). The reader accepts every spelling the historical writers
    /// produced, ignores records past the declared element count, and treats a
    /// zero extent as legal, yielding an empty [`Array`].
    ///
    /// Fails with [`Error::NotSet`] when no payload exists, [`Error::Kind`]
    /// when the attribute exists only as a scalar, [`Error::TypeMismatch`]
    /// when the registry declares another element type, [`Error::Io`] when the
    /// file exists but cannot be opened, [`Error::Malformed`] when the header
    /// declares rank 0 (the core accepts it, the format does not), and
    /// [`Error::Core`] on a corrupt header or an unparseable element (the core
    /// reports both as `Error reading data in <path>.`, without a line
    /// number).
    pub fn get_array<T: Element>(&self, group: &str, key: &str) -> Result<Array<T>> {
        Self::check_name("group", group)?;
        Self::check_name("attribute", key)?;
        // Kind and element type only: the declared extents are not resolved
        // on this path (see [`EzFio::set_array`] for the shape guard).
        let _ = self.check_array_spec::<T>(group, key)?;
        let dir = self.root.join(group);
        let path = dir.join(format!("{key}.gz"));
        if !path.is_file() {
            return Err(self.not_set_or_kind(group, key, true));
        }
        probe_readable(&path)?;

        let stored = ops::array_info(&dir, key).map_err(core_error)?;
        // A `%3d` rank of 0 parses in the core (`prod([]) == 1`), but no
        // writer produces one and the reference getters reject it.
        if stored.is_empty() {
            return Err(Error::malformed(&path, "rank 0 is outside 1..=10"));
        }
        let mut dims = Vec::with_capacity(stored.len());
        for &extent in &stored {
            // `array_info` only returns shapes whose product exists, and
            // `prod` rejects negative extents, so this cannot fire; the
            // checked conversion is kept rather than an unchecked cast.
            dims.push(usize::try_from(extent).map_err(|_| {
                Error::malformed(&path, format!("negative extent {extent} in the header"))
            })?);
        }
        let data = ops::read_array(&dir, key, T::TYPE.to_core_ty(), &stored).map_err(core_error)?;
        let values = T::from_array_data(data).ok_or_else(|| core_value_mismatch::<T>(&path))?;
        Array::new(dims, values)
    }

    /// Write the gzip array `group/key`.
    ///
    /// - `group`, `key`: EZFIO names; the group is created on demand.
    /// - `value`: dimensions (Fortran order) and flattened elements.
    ///
    /// The dimensions are checked against the registry when they resolve (the
    /// read path never checks them: see [`EzFio::get_array`]), the rank must be
    /// `1..=10` ([`MAX_RANK`], the core's limit), and every extent must fit in
    /// an `I20` field. The declared element type must accept `T`. `ezfio_core`
    /// writes the header and the records with the OCaml dialect (the private
    /// `WRITE_DIALECT`), which stores characters verbatim: the padding the
    /// Fortran writer applies to a declared `character*(n)` element is not
    /// reproduced.
    ///
    /// The gzip stream is deterministic (`mtime` 0, no embedded file name) and
    /// the file is installed with a rename, so two writes of the same array
    /// produce identical bytes.
    pub fn set_array<T: Element>(&self, group: &str, key: &str, value: &Array<T>) -> Result<()> {
        Self::check_name("group", group)?;
        Self::check_name("attribute", key)?;
        let spec = self.check_array_spec::<T>(group, key)?;
        self.check_shape(group, key, spec, value.dims())?;

        let rank = value.dims().len();
        if rank == 0 || rank > MAX_RANK {
            return Err(Error::Rank { rank });
        }
        let mut dims = Vec::with_capacity(rank);
        for (index, &extent) in value.dims().iter().enumerate() {
            let extent = i64::try_from(extent).map_err(|_| Error::DimensionTooLarge {
                index,
                value: extent,
            })?;
            dims.push(extent);
        }

        self.create_group(group)?;
        // `ArrayData` owns its elements, so this copies the array once; the
        // core takes no borrowed payload (see the module docs).
        let data = T::to_array_data(value.data());
        ops::write_array(
            &self.root.join(group),
            key,
            T::TYPE.to_core_ty(),
            &dims,
            &data,
            WRITE_DIALECT,
        )
        .map_err(core_error)
    }

    /// Registry check for a scalar access.
    fn check_scalar_spec<T: Element>(&self, group: &str, key: &str) -> Result<()> {
        let Some(spec) = self.registry.lookup(group, key) else {
            return Ok(());
        };
        if spec.is_array() {
            return Err(Error::Kind {
                group: group.to_string(),
                key: key.to_string(),
                detail: "the registry declares an array, but a scalar was requested".to_string(),
            });
        }
        check_element_type::<T>(group, key, spec)
    }

    /// Registry check for an array access; returns the spec, whose declared
    /// extents the write path needs.
    fn check_array_spec<T: Element>(&self, group: &str, key: &str) -> Result<Option<KeySpec>> {
        let Some(spec) = self.registry.lookup(group, key) else {
            return Ok(None);
        };
        if !spec.is_array() {
            return Err(Error::Kind {
                group: group.to_string(),
                key: key.to_string(),
                detail: "the registry declares a scalar, but an array was requested".to_string(),
            });
        }
        check_element_type::<T>(group, key, spec)?;
        Ok(Some(spec))
    }

    /// Compare `found` with the registry shape, when it can be resolved.
    ///
    /// Write path only ([`EzFio::set_array`]): a reader must accept a tree
    /// whose declared extents have since changed, exactly like the reference
    /// getters, which never consult the declaration.
    ///
    /// An unresolved dimension expression (a referenced attribute not set yet)
    /// means the tree is mid-construction, e.g. `pseudo/pseudo_dz_k` is written
    /// with extents `(3,0)` before `pseudo/pseudo_klocmax` exists; the check is
    /// skipped rather than failing.
    fn check_shape(
        &self,
        group: &str,
        key: &str,
        spec: Option<KeySpec>,
        found: &[usize],
    ) -> Result<()> {
        let Some(spec) = spec else {
            return Ok(());
        };
        let Some(raw_dims) = spec.dims else {
            return Ok(());
        };
        let Some(expected) = resolve_dims(raw_dims, self)? else {
            return Ok(());
        };
        if expected.as_slice() != found {
            return Err(Error::Shape {
                group: group.to_string(),
                key: key.to_string(),
                expected,
                found: found.to_vec(),
            });
        }
        Ok(())
    }

    /// The set/not-set error for an attribute, given which form was requested.
    ///
    /// `array_requested` selects the wording: when the *other* form exists, the
    /// attribute is set but of the wrong kind, which is a much better
    /// diagnostic than "not set".
    fn not_set_or_kind(&self, group: &str, key: &str, array_requested: bool) -> Error {
        let dir = self.root.join(group);
        let (own, other, detail) = if array_requested {
            (
                dir.join(format!("{key}.gz")),
                dir.join(key),
                "on disk this attribute is a plain-text scalar, not a gzip array",
            )
        } else {
            (
                dir.join(key),
                dir.join(format!("{key}.gz")),
                "on disk this attribute is a gzip array, not a plain-text scalar",
            )
        };
        if own.exists() {
            // Present but not a regular file (e.g. a directory): expose the
            // kind instead of claiming the attribute is missing.
            Error::Kind {
                group: group.to_string(),
                key: key.to_string(),
                detail: "the payload exists but is not a regular file".to_string(),
            }
        } else if other.is_file() {
            Error::Kind {
                group: group.to_string(),
                key: key.to_string(),
                detail: detail.to_string(),
            }
        } else {
            Error::NotSet {
                group: group.to_string(),
                key: key.to_string(),
            }
        }
    }

    /// Validate a group or attribute name as a single safe path component.
    ///
    /// EZFIO names are lower-case ASCII identifiers in practice, but the
    /// reference clients perform no validation at all and plugin attributes
    /// exist, so only the characters that would escape the tree or collide with
    /// the markers are rejected:
    ///
    /// * empty, `.` and `..`;
    /// * a leading `.` (reserved for `.version` and temporary files);
    /// * `/`, `\` and NUL anywhere.
    fn check_name(kind: &'static str, name: &str) -> Result<()> {
        let detail = if name.is_empty() {
            "a name cannot be empty"
        } else if name == "." || name == ".." {
            "`.` and `..` are not names"
        } else if name.starts_with('.') {
            "a name cannot start with `.` (reserved for the `.version` marker)"
        } else if name.contains('/') || name.contains('\\') || name.contains('\0') {
            "a name cannot contain a path separator or NUL"
        } else {
            return Ok(());
        };
        Err(Error::InvalidName {
            kind,
            name: name.to_string(),
            detail,
        })
    }
}

/// True when the registry's declared type accepts the Rust element type.
///
/// The EZFIO `integer` and `integer*8` declarations share one on-disk
/// representation (`I20`), and [`i64`] is the single Rust integer used for
/// both, so the two are interchangeable. Everything else must match exactly.
fn types_compatible(declared: ElementType, requested: ElementType) -> bool {
    match (declared, requested) {
        (
            ElementType::Integer | ElementType::Integer64,
            ElementType::Integer | ElementType::Integer64,
        ) => true,
        (declared, requested) => declared == requested,
    }
}

/// Apply [`types_compatible`] and build the [`Error::TypeMismatch`] otherwise.
fn check_element_type<T: Element>(group: &str, key: &str, spec: KeySpec) -> Result<()> {
    let requested = T::TYPE;
    if types_compatible(spec.ty, requested) {
        return Ok(());
    }
    Err(Error::TypeMismatch {
        group: group.to_string(),
        key: key.to_string(),
        expected: spec.ty.as_str().to_string(),
        found: requested.as_str().to_string(),
    })
}

/// Turn an `ezfio_core` failure into the crate's error type.
///
/// The core reports the failing ABI entry point (`where_`) and a message that
/// already names the file (`Error reading data in <path>.`); only the message
/// is user-facing here, so only it is kept. The underlying `io::Error` and the
/// line number of a malformed payload are not available from the core.
fn core_error(error: EzfioError) -> Error {
    Error::Core(error.txt)
}

/// The error for a value the core returned with an unexpected variant.
///
/// [`Element::from_scalar`] and [`Element::from_array_data`] return `None` when
/// the core answers with another value type. Unreachable: the crate asks the
/// core for [`Element::TYPE`]'s core type, and the core parses exactly that
/// type. Returned rather than panicked so that a future core change surfaces
/// as an error.
fn core_value_mismatch<T: Element>(path: &Path) -> Error {
    Error::Core(format!(
        "ezfio-core returned another value type for the {} request on {}",
        T::TYPE,
        path.display()
    ))
}

/// Surface the OS error a delegated read would otherwise hide.
///
/// `ezfio_core::ops` folds "cannot open" and "cannot parse" into one
/// `Rc::IoError` and discards the `io::Error`, so a file that exists but
/// cannot be opened (permissions) would look like a malformed payload. One
/// `open` keeps [`Error::Io`] for the case the caller can act on; the handle
/// is closed at once and the core re-opens the file.
fn probe_readable(path: &Path) -> Result<()> {
    File::open(path).at(path)?;
    Ok(())
}

impl DimResolver for EzFio {
    /// Read `group/key` as a scalar integer, or `Ok(None)` if it is not set.
    ///
    /// Used while evaluating the dimension expressions of the registry; the
    /// integer is parsed by `ezfio_core`, with the same tolerance as any other
    /// scalar read. A reference to an array (no plain file) is therefore
    /// unresolved, which is the conservative answer: the shape check is
    /// skipped rather than guessed.
    fn resolve(&self, group: &str, key: &str) -> Result<Option<i64>> {
        let dir = self.root.join(group);
        let path = dir.join(key);
        if !path.is_file() {
            return Ok(None);
        }
        probe_readable(&path)?;
        let value =
            ops::read_scalar(&dir, key, <i64 as Element>::TYPE.to_core_ty()).map_err(|error| {
                Error::Core(format!("dimension reference {group}/{key}: {}", error.txt))
            })?;
        <i64 as Element>::from_scalar(value)
            .map(Some)
            .ok_or_else(|| core_value_mismatch::<i64>(&path))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::{Read, Write};

    use flate2::read::GzDecoder;
    use flate2::{Compression, GzBuilder};

    /// Create a fresh EZFIO tree in a temporary directory.
    fn tree() -> (tempfile::TempDir, EzFio) {
        let dir = tempfile::tempdir().expect("temp dir");
        let ezfio = EzFio::create(dir.path().join("tree")).expect("create");
        (dir, ezfio)
    }

    /// Decompress an array file and return the payload as text.
    fn payload(ezfio: &EzFio, group: &str, key: &str) -> String {
        let path = ezfio.root().join(group).join(format!("{key}.gz"));
        let file = File::open(path).expect("open payload");
        let mut text = String::new();
        GzDecoder::new(file)
            .read_to_string(&mut text)
            .expect("gunzip");
        text
    }

    /// gzip a byte slice with the deterministic settings the core uses.
    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = GzBuilder::new()
            .mtime(0)
            .write(Vec::new(), Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn create_writes_the_version_marker() {
        let (dir, ezfio) = tree();
        assert_eq!(
            fs::read_to_string(dir.path().join("tree/.version")).unwrap(),
            "2.0.8\n"
        );
        assert_eq!(ezfio.version().unwrap(), "2.0.8");
        assert!(EzFio::is_ezfio_dir(ezfio.root()));
        assert_eq!(ezfio.groups().unwrap(), Vec::<String>::new());
    }

    #[test]
    fn open_requires_the_marker() {
        let dir = tempfile::tempdir().unwrap();
        assert!(EzFio::open(dir.path().join("missing")).is_err());
        fs::create_dir(dir.path().join("plain")).unwrap();
        let err = EzFio::open(dir.path().join("plain")).unwrap_err();
        assert!(matches!(err, Error::NotEzFio { .. }));
        // `create` on the plain directory adds the marker and succeeds.
        let ezfio = EzFio::create(dir.path().join("plain")).unwrap();
        assert!(EzFio::is_ezfio_dir(ezfio.root()));
    }

    #[test]
    fn scalar_bytes_match_the_ocaml_layout() {
        let (dir, ezfio) = tree();
        ezfio.set("t", "int", &7i64).unwrap();
        ezfio.set("t", "real", &2.5f64).unwrap();
        ezfio.set("t", "bool", &true).unwrap();
        ezfio.set("t", "string", &"hello".to_string()).unwrap();
        let read = |name: &str| fs::read_to_string(dir.path().join("tree/t").join(name)).unwrap();
        assert_eq!(read("int"), "                   7\n");
        assert_eq!(read("real"), "   2.500000000000000e+00\n");
        // The OCaml client's `print_bool` emits `T\n` and the write format
        // appends a second `\n`: the logical scalar file is 3 bytes.
        assert_eq!(read("bool"), "T\n\n");
        assert_eq!(read("string"), "hello\n");
        // Writing created the group and its marker.
        assert_eq!(
            fs::read_to_string(dir.path().join("tree/t/.version")).unwrap(),
            "2.0.8\n"
        );
    }

    #[test]
    fn scalar_roundtrip_all_types() {
        let (_dir, ezfio) = tree();
        ezfio.set("g", "i", &-42i64).unwrap();
        ezfio.set("g", "f", &-1.25e-3f64).unwrap();
        ezfio.set("g", "b", &false).unwrap();
        ezfio.set("g", "s", &"a b ".to_string()).unwrap();
        assert_eq!(ezfio.get::<i64>("g", "i").unwrap(), -42);
        assert_eq!(ezfio.get::<f64>("g", "f").unwrap(), -1.25e-3);
        assert!(!ezfio.get::<bool>("g", "b").unwrap());
        // The core trims both ends on read (`String.trim`), so the trailing
        // blank is dropped (OCaml parity) — and so would a leading one be.
        assert_eq!(ezfio.get::<String>("g", "s").unwrap(), "a b");
    }

    #[test]
    fn scalar_get_reads_fortran_and_ocaml_payloads() {
        let (_dir, ezfio) = tree();
        ezfio.set("g", "x", &0i64).unwrap();
        let path = ezfio.root().join("g/x");
        // Fortran implied exponent (`E24.15` without the `E`).
        fs::write(&path, "   0.100000000000000+101\n").unwrap();
        assert_eq!(ezfio.get::<f64>("g", "x").unwrap(), 1.0e100);
        fs::write(&path, "  .TRUE.\n").unwrap();
        assert!(ezfio.get::<bool>("g", "x").unwrap());
        fs::write(&path, "  1.307093214000000e+02   \n").unwrap();
        assert_eq!(ezfio.get::<f64>("g", "x").unwrap(), 130.7093214);
    }

    #[test]
    fn empty_scalar_file_is_a_data_error() {
        let (_dir, ezfio) = tree();
        ezfio.set("g", "x", &0i64).unwrap();
        fs::write(ezfio.root().join("g/x"), "").unwrap();
        // The previous local parser called this `Malformed("empty scalar
        // file")`; the core folds it into its single data-error message.
        let err = ezfio.get::<i64>("g", "x").unwrap_err();
        assert!(matches!(err, Error::Core(_)), "{err}");
        assert!(format!("{err}").contains("Error reading data in"), "{err}");
    }

    #[test]
    fn array_payload_matches_the_ocaml_layout() {
        let (_dir, ezfio) = tree();
        let value = Array::new(vec![2, 3], vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).unwrap();
        ezfio.set_array("g", "a", &value).unwrap();
        let text = payload(&ezfio, "g", "a");
        assert_eq!(
            text,
            concat!(
                "  2\n",
                "                   2                    3 \n",
                "   1.000000000000000e+00\n",
                "   2.000000000000000e+00\n",
                "   3.000000000000000e+00\n",
                "   4.000000000000000e+00\n",
                "   5.000000000000000e+00\n",
                "   6.000000000000000e+00\n",
            )
        );
        let back = ezfio.get_array::<f64>("g", "a").unwrap();
        assert_eq!(back, value);
    }

    /// The one array dialect this crate writes: `%20d ` dimensions (one
    /// trailing blank), `%20d` integers, `%24.15e` reals, `T`/`F` logicals,
    /// characters verbatim.
    #[test]
    fn ocaml_array_bytes_for_every_element_type() {
        let (_dir, ezfio) = tree();
        ezfio
            .set_array("g", "i", &Array::from_vec(vec![1i64, -2]))
            .unwrap();
        ezfio
            .set_array("g", "r", &Array::from_vec(vec![1.0f64, -0.5]))
            .unwrap();
        ezfio
            .set_array("g", "l", &Array::from_vec(vec![true, false]))
            .unwrap();
        ezfio
            .set_array("g", "c", &Array::from_vec(vec!["abc".to_string()]))
            .unwrap();
        assert_eq!(
            payload(&ezfio, "g", "i"),
            concat!(
                "  1\n",
                "                   2 \n",
                "                   1\n",
                "                  -2\n",
            )
        );
        assert_eq!(
            payload(&ezfio, "g", "r"),
            concat!(
                "  1\n",
                "                   2 \n",
                "   1.000000000000000e+00\n",
                // `{:.15E}` normalises to one digit before the point, so
                // -0.5 is `-5.000…e-01`, right-aligned in 24 columns.
                "  -5.000000000000000e-01\n",
            )
        );
        assert_eq!(
            payload(&ezfio, "g", "l"),
            concat!("  1\n", "                   2 \n", "T\n", "F\n")
        );
        assert_eq!(
            payload(&ezfio, "g", "c"),
            concat!("  1\n", "                   1 \n", "abc\n")
        );
    }

    /// `nuclei/nucl_label` is declared `character*(32)`; the delegated writer
    /// stores characters verbatim, so the declared width does not pad them
    /// (deviation E5: the Fortran padder is gone).
    #[test]
    fn declared_characters_are_written_verbatim() {
        let (_dir, ezfio) = tree();
        ezfio
            .set_array(
                "nuclei",
                "nucl_label",
                &Array::from_vec(vec!["H".to_string()]),
            )
            .unwrap();
        assert_eq!(
            payload(&ezfio, "nuclei", "nucl_label"),
            concat!("  1\n", "                   1 \n", "H\n")
        );
    }

    #[test]
    fn gzip_stream_is_deterministic() {
        let (_dir, ezfio) = tree();
        let value = Array::from_vec(vec![0.25f64; 1000]);
        ezfio.set_array("g", "a", &value).unwrap();
        let path = ezfio.root().join("g/a.gz");
        let first = fs::read(&path).unwrap();
        ezfio.set_array("g", "a", &value).unwrap();
        let second = fs::read(&path).unwrap();
        assert_eq!(first, second);
        // mtime field of the gzip header is zero.
        assert_eq!(&first[4..8], &[0, 0, 0, 0]);
    }

    #[test]
    fn zero_extent_roundtrips() {
        let (_dir, ezfio) = tree();
        let value = Array::<f64>::new(vec![3, 0], Vec::new()).unwrap();
        ezfio.set_array("g", "empty", &value).unwrap();
        let text = payload(&ezfio, "g", "empty");
        assert_eq!(text, "  2\n                   3                    0 \n");
        let back = ezfio.get_array::<f64>("g", "empty").unwrap();
        assert!(back.is_empty());
        assert_eq!(back.dims(), &[3, 0]);
    }

    #[test]
    fn array_reader_accepts_a_fortran_payload() {
        let (_dir, ezfio) = tree();
        fs::create_dir_all(ezfio.root().join("g")).unwrap();
        let bytes = gzip(b"  2\n                   2                    2\n   0.100000000000000E+01\n   0.200000000000000E+01\n   0.300000000000000E+01\n   0.400000000000000E+01\n");
        fs::write(ezfio.root().join("g/a.gz"), bytes).unwrap();
        let array = ezfio.get_array::<f64>("g", "a").unwrap();
        assert_eq!(array.dims(), &[2, 2]);
        assert_eq!(array.data(), &[1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn array_reader_rejects_malformed_payloads() {
        let (_dir, ezfio) = tree();
        fs::create_dir_all(ezfio.root().join("g")).unwrap();
        let path = ezfio.root().join("g/a.gz");
        let check = |payload: &str| {
            fs::write(&path, gzip(payload.as_bytes())).unwrap();
            ezfio.get_array::<i64>("g", "a").unwrap_err()
        };
        // Rank 0 is outside the format (the core parses it, this layer
        // rejects it).
        assert!(matches!(check("  0\n\n"), Error::Malformed { .. }));
        // The core accepts at most 10 dimensions; the header parses, the
        // limit check fires.
        let rank11 = format!(" 11\n{}\n", "1 ".repeat(11));
        let err = check(&rank11);
        assert!(
            format!("{err}").contains("exceeds the maximum of 10"),
            "{err}"
        );
        // Fewer dimensions than the rank announces.
        let err = check("  2\n                   1 \n");
        assert!(matches!(err, Error::Core(_)), "{err}");
        assert!(format!("{err}").contains("Error reading data in"), "{err}");
        // A non-numeric element: the core reports the file but not the line.
        let err = check("  1\n                   1 \nx\n");
        assert!(matches!(err, Error::Core(_)), "{err}");
        assert!(format!("{err}").contains("Error reading data in"), "{err}");
        // Blank trailing lines are tolerated.
        fs::write(
            &path,
            gzip(b"  1\n                   1 \n                   7\n\n\n"),
        )
        .unwrap();
        assert_eq!(ezfio.get_array::<i64>("g", "a").unwrap().data(), &[7]);
        // Records past the declared count are ignored (the core stops after
        // `prod(dims)` elements; the previous local reader rejected them).
        fs::write(
            &path,
            gzip(b"  1\n                   1 \n                   7\n                   8\n"),
        )
        .unwrap();
        assert_eq!(ezfio.get_array::<i64>("g", "a").unwrap().data(), &[7]);
    }

    #[test]
    fn get_distinguishes_not_set_from_wrong_kind() {
        let (_dir, ezfio) = tree();
        let err = ezfio.get::<i64>("g", "nope").unwrap_err();
        assert!(matches!(err, Error::NotSet { .. }));
        ezfio
            .set_array("g", "a", &Array::from_vec(vec![1i64]))
            .unwrap();
        let err = ezfio.get::<i64>("g", "a").unwrap_err();
        assert!(matches!(err, Error::Kind { .. }));
        ezfio.set("g", "s", &1i64).unwrap();
        let err = ezfio.get_array::<i64>("g", "s").unwrap_err();
        assert!(matches!(err, Error::Kind { .. }));
    }

    #[test]
    fn registry_rejects_type_contradictions() {
        let (_dir, ezfio) = tree();
        // `nuclei/nucl_num` is declared `integer`.
        assert!(matches!(
            ezfio
                .set("nuclei", "nucl_num", &"2".to_string())
                .unwrap_err(),
            Error::TypeMismatch { .. }
        ));
        assert!(matches!(
            ezfio.get::<f64>("nuclei", "nucl_num").unwrap_err(),
            Error::TypeMismatch { .. }
        ));
        // `i64` is accepted for both integer widths.
        ezfio.set("nuclei", "nucl_num", &2i64).unwrap();
        assert_eq!(ezfio.get::<i64>("nuclei", "nucl_num").unwrap(), 2);
        // `nuclei/nucl_coord` is declared as an array.
        assert!(matches!(
            ezfio.set("nuclei", "nucl_coord", &0.0f64).unwrap_err(),
            Error::Kind { .. }
        ));
        // `nuclei/nucl_label` is declared as an array of characters.
        assert!(matches!(
            ezfio.get_array::<f64>("nuclei", "nucl_label").unwrap_err(),
            Error::TypeMismatch { .. }
        ));
        // Unknown keys bypass the registry entirely.
        ezfio.set("plugin_group", "whatever", &1.5f64).unwrap();
        assert_eq!(ezfio.get::<f64>("plugin_group", "whatever").unwrap(), 1.5);
    }

    #[test]
    fn registry_checks_array_shape_on_writes_only() {
        let (_dir, ezfio) = tree();
        // `nuclei/nucl_coord` has dims `(nuclei.nucl_num, 3)`.
        ezfio.set("nuclei", "nucl_num", &2i64).unwrap();
        ezfio
            .set_array(
                "nuclei",
                "nucl_coord",
                &Array::new(vec![2, 3], vec![0.0; 6]).unwrap(),
            )
            .unwrap();
        let err = ezfio
            .set_array(
                "nuclei",
                "nucl_coord",
                &Array::new(vec![3, 3], vec![0.0; 9]).unwrap(),
            )
            .unwrap_err();
        assert!(matches!(err, Error::Shape { .. }));

        // The read path never consults the declaration: growing `nucl_num`
        // makes the stored shape stale, and the array still comes back with
        // its on-disk dims (the reference getters do not look at `size:`).
        ezfio.set("nuclei", "nucl_num", &3i64).unwrap();
        let stale = ezfio
            .get_array::<f64>("nuclei", "nucl_coord")
            .expect("a read must not fail on a stale declared shape");
        assert_eq!(stale.dims(), &[2, 3]);
        assert_eq!(stale.data().len(), 6);

        // ... while the write guard still rejects the stale shape: the
        // leniency is on reads only.
        let err = ezfio
            .set_array(
                "nuclei",
                "nucl_coord",
                &Array::new(vec![2, 3], vec![0.0; 6]).unwrap(),
            )
            .unwrap_err();
        assert!(
            matches!(err, Error::Shape { .. }),
            "the write guard must not be weakened: {err}"
        );

        // Without `nucl_num`, the expression is unresolved: no shape check.
        let (dir, ezfio) = tree();
        fs::create_dir_all(dir.path().join("tree/nuclei")).unwrap();
        ezfio
            .set_array(
                "nuclei",
                "nucl_coord",
                &Array::new(vec![3, 3], vec![0.0; 9]).unwrap(),
            )
            .unwrap();
    }

    #[test]
    fn groups_keys_and_removal() {
        let (_dir, ezfio) = tree();
        ezfio.set("b", "x", &1i64).unwrap();
        ezfio.set("a", "y", &1i64).unwrap();
        ezfio
            .set_array("a", "z", &Array::from_vec(vec![1i64]))
            .unwrap();
        assert_eq!(ezfio.groups().unwrap(), vec!["a", "b"]);
        assert_eq!(ezfio.keys("a").unwrap(), vec!["y", "z"]);
        assert!(ezfio.keys("missing").unwrap().is_empty());
        assert!(ezfio.has("a", "y").unwrap());
        assert!(ezfio.has("a", "z").unwrap());
        assert!(!ezfio.has("a", "nope").unwrap());
        assert!(!ezfio.has("missing", "y").unwrap());
        assert!(ezfio.remove("a", "z").unwrap());
        assert!(!ezfio.remove("a", "z").unwrap());
        assert!(ezfio.has("a", "y").unwrap());
        // The group marker is not listed and cannot be removed as a key.
        assert_eq!(ezfio.keys("a").unwrap(), vec!["y"]);
        assert!(ezfio.remove("a", ".version").is_err());
    }

    /// A declared key is probed in its declared form only; a registry-unknown
    /// key probes both forms (plugin attributes).
    #[test]
    fn has_follows_the_declared_kind() {
        let (_dir, ezfio) = tree();
        ezfio.set("nuclei", "nucl_num", &1i64).unwrap();
        assert!(ezfio.has("nuclei", "nucl_num").unwrap());
        ezfio
            .set_array(
                "nuclei",
                "nucl_coord",
                &Array::new(vec![1, 3], vec![0.0; 3]).unwrap(),
            )
            .unwrap();
        assert!(ezfio.has("nuclei", "nucl_coord").unwrap());
        ezfio.set("plugin", "x", &1i64).unwrap();
        assert!(ezfio.has("plugin", "x").unwrap());
    }

    #[test]
    fn names_cannot_escape_the_tree() {
        let (_dir, ezfio) = tree();
        for name in ["", ".", "..", ".version", "a/b", "a\\b"] {
            assert!(
                ezfio.set(name, "k", &1i64).is_err(),
                "group {name:?} must be rejected"
            );
            assert!(
                ezfio.set("g", name, &1i64).is_err(),
                "key {name:?} must be rejected"
            );
        }
        assert!(ezfio.get::<i64>("../x", "k").is_err());
    }

    #[test]
    fn atomic_writes_leave_no_temporary_files() {
        let (_dir, ezfio) = tree();
        ezfio.set("g", "x", &1i64).unwrap();
        ezfio
            .set_array("g", "a", &Array::from_vec(vec![1i64]))
            .unwrap();
        let mut entries: Vec<String> = fs::read_dir(ezfio.root().join("g"))
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        entries.sort();
        assert_eq!(entries, vec![".version", "a.gz", "x"]);
    }

    #[test]
    fn dimension_resolution_reads_scalars_only() {
        let (_dir, ezfio) = tree();
        ezfio.set("determinants", "n_int", &3i64).unwrap();
        assert_eq!(ezfio.resolve("determinants", "n_int").unwrap(), Some(3));
        assert_eq!(ezfio.resolve("determinants", "nope").unwrap(), None);
        // A reference that is set but does not hold an integer is an error
        // (the registry declares `bit_kind` as an integer, so write it behind
        // the typed API's back to exercise the resolver alone).
        fs::write(ezfio.root().join("determinants/bit_kind"), "x\n").unwrap();
        assert!(ezfio.resolve("determinants", "bit_kind").is_err());
    }

    #[test]
    fn invalid_rank_is_rejected_on_write() {
        // `Array` cannot hold rank 0, so only the > 10 side is reachable.
        let (_dir, ezfio) = tree();
        let dims = vec![1usize; MAX_RANK + 1];
        let value = Array::new(dims, vec![0i64; 1]).unwrap();
        assert!(matches!(
            ezfio.set_array("g", "a", &value).unwrap_err(),
            Error::Rank { rank } if rank == MAX_RANK + 1
        ));
    }
}
