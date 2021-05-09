==================
Coding conventions
==================


General conventions
===================

All executable files should have a name with lowercase.

Tabs are forbidden everywhere.

Try to set the maximum line length to 80 characters.  Long lines can be
automatically reformatted in vim by pressing :kbd:`gqj`.

Use blank lines between blocks to improve readability.

For existing files, stay faithful to the existing indentation.



Shell scripts
=============

Executables should have no extension.  To know if the file is binary, or in
what shell scripting language it was written, the :command:`file` command can
be used. In addition, all the shell scripts should be under
:file:`${QP_ROOT}/scripts/`.

The exit code of the script should be 0 upon success only.

Bash and Python3 are the only shell scripting language permitted for
executables.


Bash
----

* Bash scripts should start with ``#!/bin/bash``

* All error messages should go to standard error, and should be prefixed with
  the name of the command. For example, in Bash use

  .. code:: bash

      function echo_err() {
         2>& echo $(basename $0)": error"
      }

* The command-line options should be handled with ``getopt``.

* The script should check that the command-line arguments are consistent.

* Long options should be preferred to short options.

* Always quote strings containing variables, command substitutions, spaces or
  shell meta characters, unless careful unquoted expansion is required.

* Use ``"$@"`` unless you have a specific reason to use ``$*``.

* Use ``$(command)`` instead of backticks, because they can be easily nested.

* ``[[ ... ]]`` is preferred over ``[``, ``test`` and ``/usr/bin/[``.

* Declare function-specific variables with local. Declaration and assignment
  should be on different lines.

* Pipelines should be split one per line if they don't all fit on one line.

* Put ``; do`` and ``; then`` on the same line as the ``while``, ``for`` or ``if``.


Python
------

Only Python3 is supported.

Python scripts should start with ``#!/usr/bin/env python3`` to mention
explicitly that Python3 has to be used, possibly from a conda installation.

:command:`pylint` should be used to increase the quality of the source code.



IRPF90
======

The code can be automatically indented with :command:`irp_indent`.

Lines sould not be longer than 80 characters.

Mathematical formulas in the `BEGIN_DOC...END_DOC` sections sould be written in
LaTeX format, between `$` symbols.

All the providers, subroutines and functions should have a
`BEGIN_DOC...END_DOC` block.

Providers should never be present in the same file as a main program.

String must not use double quotes (`"`) but single quotes (`'`).

After a `read` statement there should be no comma.

Only standard Fortran is allowed : Intel or GNU extensions are forbidden.

The name of a program should be the same as the name of the file. For example,
for the :ref:`fci` program, we have

.. code-block:: fortran

   program fci

and the file is named :file:`fci.irp.f`.

