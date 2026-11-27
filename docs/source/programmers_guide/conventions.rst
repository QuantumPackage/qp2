==================
Coding conventions
==================


General conventions
===================

All program names and executable file names are in lowercase.

Tabs are forbidden everywhere.

Try to keep the maximum line length to 80 characters. Long lines can be
automatically rewrapped in vim with :kbd:`gqj`.

Use blank lines between blocks to improve readability.

For existing files, stay faithful to the existing indentation.


Shell scripts
=============

The command-line tools live in :file:`$QP_ROOT/bin/` and have no
extension. The helper scripts live under :file:`$QP_ROOT/scripts/`.
The :command:`file` command tells whether a file is binary and in which
language it is written.

Bash and Python3 are the only scripting languages used by the
repository.

The exit code of a script should be 0 upon success only.


Bash
----

* Bash scripts should start with ``#!/bin/bash``.

* All error messages should go to standard error, and should be
  prefixed with the name of the command, as in
  :file:`$QP_ROOT/bin/qp_export_as_tgz`::

      function error()
      {
          >&2 echo "$(basename $0): $@"
      }

* The script should check that the command-line arguments are
  consistent.

* Long options should be preferred to short options.

* Always quote strings containing variables, command substitutions,
  spaces or shell meta characters, unless careful unquoted expansion is
  required.

* Use ``"$@"`` unless you have a specific reason to use ``$*``.

* Use ``$(command)`` instead of backticks, because they can be easily
  nested.

* ``[[ ... ]]`` is preferred over ``[``, ``test`` and ``/usr/bin/[``.

* Declare function-specific variables with :command:`local`; the
  declaration and the assignment should be on different lines.

* Pipelines should be split one per line if they don't all fit on one
  line.

* Put ``; do`` and ``; then`` on the same line as the ``while``,
  ``for`` or ``if``.


Python
------

Only Python3 is supported.

Python scripts should start with ``#!/usr/bin/env python3`` to mention
explicitly that Python3 has to be used, possibly from a conda
installation.

The Python scripts of the repository parse their command line with
``docopt``: the usage message is written in the docstring of the
script, and the options are derived from it.

:command:`pylint` should be used to increase the quality of the source
code.


IRPF90
======

The code can be automatically indented with :command:`irpf90_indent`.

Lines should not be longer than 80 characters.

Mathematical formulas in the ``BEGIN_DOC``...``END_DOC`` sections
should be written in the LaTeX format, between ``$`` symbols.

All the providers, subroutines and functions should have a
``BEGIN_DOC``...``END_DOC`` block; the preprocessor prints a warning
when one is missing.

Providers should never be present in the same file as a main program.

Strings should use single quotes (``'``) rather than double quotes
(``"``).

A ``read`` unit list is not followed by a comma: write
``read (iunit) a, b``, not ``read (iunit), a, b``.

Only standard Fortran is allowed; Intel or GNU extensions are
forbidden.

The name of a program should be the same as the name of the file. For
example, the :ref:`fci` program is declared as

.. code-block:: fortran

    program fci

in the file :file:`fci.irp.f`. This matters because the executable
takes its name from the file, not from the ``program`` statement: a
file named :file:`foo.irp.f` produces the executable :file:`foo`,
whatever the name of the program it contains.

Integer kinds
 -----------

Integers that index or size system-sized data (AO/MO counts,
integral-map keys, determinant-space indices, grid-point counts)
use the 64-bit kind ``integer(size_t)``, where ``size_t = 8`` is
defined by the hand-written module :file:`src/utils/sizet.f90`.
Small chemistry-bounded or machine-bounded values (state counts,
CAS orbital counts, irreps, flags, file descriptors, thread
counts) keep the default 32-bit ``integer``.

Placement of ``use sizet``:

* ``.irp.f`` providers: at file level (IRPF90 hoists it into the
  generated module scope).
* hand-written procedures in ``.irp.f``: inside the procedure
  body (IRPF90 does not hoist file-level uses into them).
* program files: inside the ``program`` body, never at file level
  — a file-level use in a program file lands in the generated
  wrapper module, whose compile edge has no dependency on
  :file:`sizet.o` and races on :file:`sizet.mod`.
* hand-written ``.f90``: at module scope. The five
  directly-compiled ``.f90`` files built in an unordered ninja
  pool (``selection_types.f90``, ``qp_json.f90``,
  ``pair_block_stats.f90``, ``map_module.f90``, ``spin_exc.F90``)
  spell the kind ``integer(8)`` or ``integer(c_int64_t)`` instead,
  which is type-identical.

Kinds never mix across a call boundary: gfortran rejects an
actual/dummy integer kind mismatch on explicit-interface calls
and only warns (while silently corrupting the value) on
same-file calls. Cast at the call site — ``int(x, size_t)`` up,
``int(x, 4)`` down — and name 32-bit boundary temporaries
``*_4``. The same casts guard the 32-bit external interfaces
(BLAS/LAPACK, MPI, libcint, json-C, the EZFIO file format, all
``bind(C)`` code); those interfaces stay 32-bit.

IRPF90 cannot parse a typed function header such as
``integer(size_t) function f(a)`` (it exits 0 and the damage is
in the generated file). Declare the result in the body instead::

    function f(a)
      implicit none
      integer(size_t) :: f
      ...
    end function f
