=======================
Programming in the |qp|
=======================

To program in the |qp|, it is required that you are familiar with the |IRPF90|
code generator. A GitBook can be found `here <http://scemama.gitbooks.io/irpf90>`_,
and programmers are encouraged to visit this manual.

|IRPF90| makes programming very simple. The only information a programmer needs
in order to write a new program is the name of the required |IRPF90| entities
which may already exist in other modules.  For example, writing a program which
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


The only required information was the existence of a provider for
:command:`HF_energy`. A detailed list of all the providers, subroutines
and functions of the |qp| can be found in the appendix of this manual.



Architecture
============

As |IRPF90| is used, the programmer doesn't have full control of the sequence
of instructions in the produced Fortran code. This explains why the input data
is stored in a database rather than in sequential text files. Consequently, the
programmer can't know in advance the order in which the files will be read, so a
simple random access to persistent data is needed. The |EZFIO| library generator
is a practical answer to this problem. 

The |qp| uses a collection of programs inter-operating together. Each of these
programs is reading and/or modifying information in the |EZFIO| database.
This is done mostly using the command line or scripting.

.. important::

    Each command modifies the state of the |EZFIO| database, so running the
    same program twice on the same database may have different behavior because of the
    state of the database. For reproducibility, users are encouraged to run scripts
    where a fresh new |EZFIO| database is created at the beginning of the
    script. This way of running the |qp| makes calculations reproducible.


The computational part |qp| is organized in **modules**. A module is a
directory which contains multiple |IRPF90| files, a |README| and a |NEED| file.

The |README| file contains documentation about the module, that is
automatically included in the documentation of the |qp|. The documentation is
generated by the `Sphinx documentation builder <http://www.sphinx-doc.org>`_,
and it should be written using the |rst| format.

The |NEED| file contains the list of the modules which are needed for the
current module. When a module is needed, it means that all the |IRPF90| files
it contains should be included in the current module. This is done
automatically during the building process, by creating symbolic links in the
current directory.

To compile the program, the |Ninja| build system is used, and all the building
process is fully automated such that the programmer will never have to modify a
file by hand. Running :command:`ninja` inside a module will compile only the
module, and running :command:`ninja` at the root of the |qp| will build all the
modules, as well as the tools.


Algorithms
==========

The `PhD thesis of Yann Garniron <https://doi.org/10.5281/zenodo.2558127>`_
gives all the details about the implementation of:

* The data structure for the two-electron integrals (:file:`utils/map_module.f`)
* The Davdison diagonalization (module :ref:`module_davidson`)
* The CIPSI selection (module :ref:`module_cipsi`)
* The hybrid stochastic/deterministic PT2 correction (module :ref:`module_cipsi`)
* The hybrid stochastic/deterministic matrix dressing (module :ref:`module_dressing`)


Extracting results for use with other codes
===========================================

The |AOs| and |MOs| can be seen with :ref:`qp_edit`. We also provide a utility
to create a file which can be read by `molden` for visualizing the |MOs| (see
:ref:`molden`). For using external |CI| solvers, we provide a utility that
generates a file containing the two-electron integrals in the |MO| basis set
in the `FCIDUMP` format (see :ref:`fcidump`).

All the results are stored in the |EZFIO| directory, so users willing to fetch
data such as the |MOs| or the |CI| coefficients should use the |EZFIO| API.
There are multiple major ways to do this:

* Write a script in Python or OCaml and use the Python |EZFIO| API. The script
  :file:`$QP_ROOT/bin/qp_convert_output_to_ezfio` is a good example to understand
  how to use the |EZFIO| API in Python,
* Write an independent program in Fortran or C, link it with the |EZFIO| library
  located at :file:`$QP_ROOT/external/ezfio/lib/libezfio.a` and call directly
  the |EZFIO| routines,
* Write a new module for the |qp| printing the desired quantities in a suitable
  text format. The program :ref:`fcidump` is an example of such a program.


