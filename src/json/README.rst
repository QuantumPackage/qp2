====
json
====

JSON result document of a calculation.

Each run of a QP3 program writes one document holding the EZFIO inputs the
calculation inspected and the numbers it produced.  Consumers (TREX2 and
friends) read that document instead of parsing the standard output.  The
file is valid JSON at every instant, including after ``kill -9``; validity
never depends on the program reaching an orderly exit.

Logging a value
===============

One line at the call site::

  call json_add('/results/scf/energy', SCF_energy)

and, once per program, before the first ``json_add``::

  call json_init('scf')

``json_init`` (``src/json/json.irp.f``) chooses the file name, opens the
document and writes its envelope.  ``json_add`` is the generic of
``src/json/qp_json.f90``: integers, reals, logicals, character strings and
rank-1 integer or real arrays.  An unsupported actual argument (single
precision, ``integer(8)``, a rank-2 array) is a compile error by design.

Paths are absolute, ``snake_case`` and contain no spaces.  Intermediate
objects are created on the way; logging the same path twice replaces the
value, so a call site may re-log a key every iteration.  A malformed path
prints one warning and is dropped; the run continues.

The file
========

``<EZFIO>/json/NNNNN.json``, ``NNNNN`` being the first unused number.  A
restart gets the next number and never touches the previous file.  Every
write replaces the whole file, so the newest document is a complete,
self-contained record of that run: a ``program fci`` run contains the SCF
results computed in process as well as the CIPSI results.

Root keys, all written by ``qp_json_begin``:

``schema``
  ``"qp3/1"``.  Anything else is not a document of this format.

``creator``
  ``{"name": "Quantum Package", "version": ...}``.  ``version`` reads
  ``"unknown"`` until a build-time ``git describe`` string is generated.

``provenance``
  ``ezfio_file``, ``host``, ``pid``, ``started`` (UTC, seconds).

``units``
  ``energy`` in hartree, ``length`` in bohr, ``time`` in second.

``run``
  ``kind`` is the program name passed to ``json_init``, ``timestamp`` and
  ``elapsed`` advance once per actual write and describe the run so far.

``input``
  One object per EZFIO group, one key per scalar input actually read, for
  example ``input.determinants.n_states``.  Generated automatically: the
  EZFIO provider generator emits the ``json_add`` call next to the
  ``ezfio_get_`` call, so only inputs the calculation really read appear.

``results``
  One object per method (``scf``, ``cipsi``, ``casscf``, ``mo_opt``,
  ``ks``), written by the call sites.

``null`` means "not finite or not available": NaN and Inf are logged as
``null``, and a consumer must treat a missing key as unavailable rather
than assume it exists.

Atomic writes
=============

Each write builds the text in memory and writes
``<path>.tmp.<pid>`` in the same directory, ``fflush``, ``fsync``,
``fclose``, then ``rename(2)`` onto ``<path>``.  ``rename`` within a
directory is atomic, so a reader sees either the previous version or the
new one, never a truncated document.  Readers accept only
``^[0-9]{5}\.json$``; a ``.tmp.<pid>`` left by ``SIGKILL`` is garbage.

``QP_JSON_MIN_INTERVAL`` (seconds, floating point) sets the maximum
staleness:

===============  =====================================================
value            behaviour
===============  =====================================================
``> 0`` (1.0)    writer thread, at most one write per interval
``0``            no thread, a synchronous write after every ``json_add``
invalid, < 0     default 1.0
===============  =====================================================

A ``SIGKILL`` therefore leaves the last throttled write: valid and at most
one interval old.  A normal exit, including ``stop 1``, forces a final
write from an ``atexit`` handler.  Non-master MPI ranks build no document
and do no I/O.  A program that calls ``json_add`` but never ``json_init``
writes nothing and warns once at exit.

Memory
======

The document lives in yyjson's per-document memory pool.  Overwriting a path
abandons the previous value there, so memory grows with the number of
``json_add`` calls, not with the number of distinct keys; it is released when
the process exits.  Measured with ``QP_JSON_MIN_INTERVAL=0`` over 10^4
iterations of the hot site (80 values and 10 keys per CIPSI iteration):
4.2 kB per iteration, so 42 MB after 10^4 iterations.  Real runs are far
below that.  A 22 570-determinant ``fci`` run on water/cc-pVDZ (13
iterations, 63 s, this machine) peaks at 82 MiB resident, and the document
itself is 4.5 kB.  Only a run of 10^5 iterations would need the path-to-value
cache of ``DESIGN.md`` §7.9; nothing in this tree runs that long.

Reading several files
=====================

When a calculation is restarted, several documents exist.  The rule is:

1. enumerate ``<EZFIO>/json/[0-9][0-9][0-9][0-9][0-9].json``, sorted by
   number, and ignore anything that is not ``"schema": "qp3/1"``;
2. take the highest-numbered file as the state of the calculation; it is
   complete, including in-process SCF and the inputs the restart read;
3. never merge fields from two files: older files are history only, and
   the sequence of ``results.cipsi.states`` ordered by number is the CIPSI
   trajectory across restarts;
4. ``run.kind`` says which program wrote a file,
   ``provenance.ezfio_file`` which calculation it belongs to and
   ``run.elapsed`` how long that run took.

``scripts/utility/qp_json.py`` implements this rule (``files``, ``load``,
``load_all``, ``load_last``); ``scripts/qp_extract_cipsi_data.py`` and
``scripts/qp_exc_energy.py`` use it.

Iteration counters
==================

Consumers must not assume one numbering convention: the SCF ``iteration``
and the CASSCF ``iteration`` are 1-based, while ``results.mo_opt.iteration``
is the 0-based CIPSI cycle number printed by its banner.  They come from
different loops.

Implementation and licensing
============================

Layers: ``src/json/json.irp.f`` (IRPF90 glue), ``src/json/qp_json.f90``
(plain Fortran generic), ``src/json/qp_json_shim.c`` (document, path
resolution, envelope, throttle, writer thread) and ``src/json/yyjson.c``,
``src/json/yyjson.h``, a flat copy of yyjson 0.13.0 under the MIT licence.

yyjson's ``dtoa`` must print doubles that read back bit-exactly, and
``CFLAGS = -Ofast`` implies ``-ffast-math``.  The phase-0 gate compiled
those two files with the build's own flags and logged ``0.1d0``,
``1.d0/3.d0``, ``1.d-300``, ``1.d-308``, ``huge(1.d0)`` and ``-0.d0``;
``jq`` read every value back equal to its ``%.17g`` image.  The
measurement is why the build keeps ``-Ofast``.

The former text writers (``json_unit``, ``json_close``, the
``json_formats.irp.f`` format providers) are gone.  A calculation adds a
number with one ``json_add`` call and no open/close bookkeeping.
