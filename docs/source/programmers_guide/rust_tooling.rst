.. _rust_tooling:

=====================
The Rust tool layer
=====================

The command-line tools of |qp| are written in Rust.  The Fortran and
|IRPF90| modules are untouched, and the |EZFIO| on-disk format is
unchanged: a directory written by :command:`qp create_ezfio` is read by
the Fortran programs, and the tools read what the Fortran programs
wrote.  The port is documented in :file:`rust/README.md` and
:file:`rust/CONTRACTS.md`; every user-visible difference from
the former |OCaml| implementation, with its reason, is listed in
:file:`rust/DEVIATIONS.md`.

The seven tools
===============

* ``qp_create_ezfio`` — build an |EZFIO| directory from an xyz file or a
  z-matrix.
* ``qp_edit`` — read the parameters as an RST document, let the user
  modify them, write them back; ``-c`` checks the input instead.
* ``qp_run`` — run a |qp| binary on an |EZFIO| directory; ``-s`` is the
  slave mode of a distributed calculation, and the same binary hosts the
  task server.
* ``qp_set_mo_class`` — assign the core, active, inactive, virtual and
  deleted classes of the MOs.
* ``qp_tunnel`` — open an ssh tunnel for slave mode.
* ``qp_print_basis`` — print the AO and MO bases of an |EZFIO| directory
  as RST.
* ``element_create_db`` — generate the element database used by
  ``qp_create_ezfio``.

Crates
======

The workspace lives in :file:`rust/` and is split so that a change to
one tool does not rebuild the others:

* ``qp-ezfio`` — native |EZFIO| on-disk format I/O (gzipped scalar and
  array files) and the typed registry.
* ``qp-types`` — physics-adjacent data types (atoms, molecules, basis
  functions, point group, ...).
* ``qp-input`` — the declarative configuration, the generic ``qp edit``
  engine and the seven hand-written section handlers.
* ``qp-zmq`` — the plain-text message protocol and the task server
  shared with the Fortran workers.
* ``qp-cli`` — program lookup and ``$QP_ROOT`` resolution.

Each binary has its own crate under :file:`rust/bins/`; the binary name
is the crate name (``qp_run``, ``qp_edit``, ...).

.. _rust_tooling_build:

Building
========

The toolchain version is pinned by :file:`rust/rust-toolchain.toml`
(Rust 1.92.0, edition 2024), with the ``rustfmt`` and ``clippy``
components.  ``cargo`` downloads the pinned toolchain on first use.

The normal way is the |Ninja| build: at the root of the repository,

.. code:: bash

   ninja

builds the modules and then the tools.  The ``rust_target`` rule runs
``make -C $QP_ROOT/ocaml``, which compiles the workspace in release mode
into :file:`rust/target/` and publishes the results as symbolic links in
:file:`ocaml/`, under both the plain names and the historical
``qp_edit.native`` / ``qp_run.native`` / ``element_create_db.byte``
spellings that :file:`bin/qp_export_as_tgz` still collects.  The
directory :file:`ocaml/` is on the :envvar:`PATH` after sourcing
:file:`quantum_package.rc`.

To build and test the tools alone, without the Fortran modules:

.. code:: bash

   source quantum_package.rc
   cd rust
   cargo build --release          # binaries in rust/target/release/
   cargo test --workspace -j4     # unit, golden and integration tests
   cargo clippy --workspace --all-targets -j4 -- -D warnings
   cargo fmt --all --check

The golden corpus (:file:`rust/tests/golden/`) holds reference |EZFIO|
fixtures, the reference ``qp edit -c`` transcripts and the basis-MD5
corpus; :file:`rust/tests/golden/README.md` gives the replay and diff
commands.  The byte-level contracts are in :file:`rust/CONTRACTS.md`.

The configuration read by :command:`qp edit`
============================================

:command:`qp edit` is a *fixed* binary.  The set of editable sections,
their fields, types, defaults and help texts are read at run time from
:file:`$QP_ROOT/data/qp_edit_config/`, which the build generates from
the same :file:`EZFIO.cfg` files that drive the Fortran side:

* :file:`registry.toml` lists the modules in visit order (the order of
  the sections of the document);
* :file:`<module>.toml` describes one module: its title, whether it is
  visited as a section, the handler that owns it, and one entry per
  field with the provider name, the |EZFIO| key, the declared type, the
  interface list, the default value, the size and the documentation.

Both files are build artifacts: they are gitignored, and editing an
:file:`EZFIO.cfg` followed by a rebuild is enough to change what
:command:`qp edit` shows.  The binary is not recompiled, so the
configuration cannot drift from the sources.

The hand-written behaviours that are not a field table (basis ordering,
the MD5 stamps, the RST and sexp parsers of the hand-edited sections,
the derived bit masks) remain imperative, in Rust, behind a named
handler.  The correspondence between the generator's handler names and
the Rust handlers is checked at build time.

When a key is absent from an |EZFIO| directory, the value comes from
:file:`$QP_ROOT/data/ezfio_defaults/<section>.ezfio_interface_default`,
exactly as in the |OCaml| implementation; the default carried by the
TOML file is only a fallback for a tree without the generated data.
The format and the lookup rules are documented in
:file:`rust/DESIGN-qp_edit-config.md`.

The configuration directory is found, in this order, from
:envvar:`QP_EDIT_CONFIG` and then from the sibling :file:`data/` directory
of the executable.  :file:`quantum_package.rc` exports
:envvar:`QP_EDIT_CONFIG` (in :file:`etc/paths.rc`), so the editor finds
the generated directory from any working directory.  The maintenance
modes (``qp_edit --list``, ``qp_edit --describe --module <NAME>``,
``qp_edit --check-config``) accept an explicit ``--config <DIR>``; they
are additions of the port and are not part of the reference command
line.

Testing
=======

``cargo test`` covers the crates.  The end-to-end tests of the tools are
the |Bats| tests of the modules, which :command:`qp_test` runs
(:ref:`qp_test`); they exercise the tools on real |EZFIO| directories.

Documentation
=============

The HTML pages and the man pages are built by Sphinx from
:file:`docs/source/`; the pages in :file:`man/` are generated output and
must not be edited by hand.  From the root of the repository, with the
documentation requirements installed (:file:`docs/INSTALL.md`):

.. code:: bash

    source quantum_package.rc
    cd docs
    make auto     # module pages, from the built Fortran sources
    make html
    make man      # the man pages, in docs/build/man/
