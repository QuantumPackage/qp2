=======================
Programming in the |qp|
=======================

To program in the |qp|, it is required that you are familiar with the
|IRPF90| code generator, for which a manual is available at
`http://irpf90.ups-tlse.fr <http://irpf90.ups-tlse.fr>`_. Programmers
are encouraged to read it before writing a new module.

|IRPF90| makes programming very simple. The only information a
programmer needs in order to write a new program is the list of the
|IRPF90| entities (providers, subroutines and functions) it uses, which
may already exist in other modules. For example, writing a program which
prints the Hartree-Fock energy is as simple as:

.. code:: fortran

    program print_hf_energy
      implicit none
      BEGIN_DOC
    ! Program which prints the Hartree-Fock energy
    ! to the standard output
      END_DOC
      print *, 'HF energy = ', HF_energy
    end

The only required information is the existence of a provider for
``HF_energy``, provided by the :ref:`module_hartree_fock` module. The
|IRPF90| preprocessor inserts the call that ensures the provider is
computed before it is read. A detailed list of all the providers,
subroutines and functions of the |qp| can be found in the index of this
manual.


Architecture
============

As |IRPF90| is used, the programmer doesn't have full control of the
sequence of instructions in the produced Fortran code. This explains
why the input data is stored in a database rather than in sequential
text files. Consequently, the programmer can't know in advance the order
in which the files will be read, so a simple random access to the
persistent data is needed. The |EZFIO| library generator is a
practical answer to this problem.

The |qp| uses a collection of programs inter-operating together. Each
of these programs reads and/or modifies information in the |EZFIO|
database, most of the time from the command line or from scripts.

.. important::

    Each command modifies the state of the |EZFIO| database, so running
    the same program twice on the same database can have different
    behavior. For reproducibility, a script should build a fresh
    |EZFIO| directory at the beginning of the script rather than reusing
    a modified one.

Modules
-------

The computational part of the |qp| is organized in **modules**. A module
is a directory that contains the |IRPF90| sources of the module
(``*.irp.f`` files, and optionally ``*.template.f`` and
``*.include.f`` codelet templates and ``*.ezfio_config`` group files),
its :file:`NEED` file, the :file:`README.rst` that documents the module,
and optionally a :file:`EZFIO.cfg` file that declares its control
parameters (see :ref:`ezfio_cfg`). The documentation is built with
`Sphinx <http://www.sphinx-doc.org>`_, and the :file:`README.rst` files
are written in the |rst| format.

The :file:`NEED` file contains one module name per line. When the module
is built, each module listed there is symlinked into the module
directory, and all its |IRPF90| files are compiled in, so the needed
entities are directly available. As a consequence, a text search inside
a module directory also scans the symlinked modules. To search the real
sources of the repository, use :command:`git grep` rather than a plain
:command:`grep`.

Build
-----

The build is driven by the |Ninja| build system and is fully automated,
such that the programmer never has to touch a build file by hand.
Running :command:`ninja` at the root of the repository builds all the
modules and then the |OCaml| command-line tools, while running
:command:`ninja` inside a module directory builds only that module (and
the |OCaml| tools). Every build first re-runs the build-graph generator,
so new files, new modules and modified :file:`NEED` files are picked up
automatically.

Files that the build creates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The build creates a number of files that must not be edited by hand and
simply disappear and reappear when the sources change:

* in each module directory that is built as a target: the symbolic links
  to the needed modules, a :file:`.gitignore` file, the
  :file:`ezfio_interface.irp.f` file generated from :file:`EZFIO.cfg`,
  the :file:`IRPF90_temp/` directory that holds the generated Fortran
  code, and a :file:`build.ninja` file;
* at the root of the repository: the :file:`build.ninja` and
  :file:`config/build.ninja` files, :file:`data/executables` (the list
  of the built programs), the generated |OCaml| files
  (``ocaml/qp_edit.ml``, ``ocaml/Input_*.ml``, ``ocaml/Qptypes.ml``,
  ``ocaml/ezfio.ml`` and their compiled ``.native`` versions), and the
  content of the :file:`lib/`, :file:`include/` and :file:`share/`
  directories.

The rule is: edit the sources, then build; never patch the generated
output.

Algorithms
==========

The `PhD thesis of Yann Garniron
<https://doi.org/10.5281/zenodo.2558127>`_ gives all the details about
the implementation of:

* The data structure for the two-electron integrals
  (:file:`src/utils/map_module.f90`)
* The Davidson diagonalization (module :ref:`module_davidson`)
* The CIPSI selection (module :ref:`module_cipsi`)
* The hybrid stochastic/deterministic PT2 correction (module
  :ref:`module_cipsi`)
* The hybrid stochastic/deterministic matrix dressing (module
  :ref:`module_dressing`)


Extracting results for use with other codes
===========================================

The |AOs| and |MOs| can be inspected with :ref:`qp_edit`. We also
provide a utility to create a file which can be read by Molden for
visualizing the |MOs| (see :ref:`molden`). For using external |CI|
solvers, we provide a utility that generates a file containing the
two-electron integrals in the |MO| basis set in the ``FCIDUMP`` format
(see :ref:`fcidump`).

All the results are stored in the |EZFIO| directory, so users willing to
fetch data such as the |MOs| or the |CI| coefficients should use the
|EZFIO| API. There are three major ways to do this:

* Write a script in Python or |OCaml| and use the |EZFIO| API. The
  script :file:`$QP_ROOT/bin/qp_convert_output_to_ezfio` is a good
  example to understand how to use the |EZFIO| API in Python,
* Write an independent program in Fortran or C, link it with the
  |EZFIO| library located at
  :file:`$QP_ROOT/external/ezfio/lib/libezfio.a` and call directly the
  |EZFIO| routines,
* Write a new module for the |qp| printing the desired quantities in a
  suitable text format. The program :ref:`fcidump` is an example of
  such a program.
