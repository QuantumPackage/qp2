.. _ezfio_cfg:

=====
EZFIO
=====

An |EZFIO| "file" is a directory: one subdirectory per |EZFIO| group
(normally the name of a module), one plain-text file per parameter, a
:file:`.version` file and a :file:`work/` directory that holds large
scratch files and can be deleted before archiving. All the values are
stored in atomic units (Bohr, Hartree).

The current |EZFIO| directory is set with::

   qp set_file my_molecule.ezfio

and stored in the ``$EZFIO_FILE`` environment variable. The parameters
are then read and written with::

   qp get                          # list the groups
   qp get determinants             # list the parameters of a group
   qp get determinants n_det_max   # read one value
   qp set determinants n_det_max 1000000
   qp has determinants n_det_max   # success (0) if the parameter is set

:command:`qp edit <directory>` opens a text editor to modify the
parameters by hand; it can also be used in batch mode, see
:ref:`qp_edit`.

.. note::

    As each program modifies the state of the |EZFIO| directory, a
    reproducible script should build a fresh |EZFIO| directory with
    :command:`qp create_ezfio` at the beginning of the script, rather
    than reusing a modified one.

EZFIO.cfg
=========

The simplest way to declare a new control parameter of the |qp| is to
create an :file:`EZFIO.cfg` file in the module. The example below is
extracted from the :file:`EZFIO.cfg` of the :ref:`module_scf_utils`
module::

   [max_dim_diis]
   type: integer
   doc: Maximum size of the |DIIS| extrapolation procedure
   interface: ezfio,provider,ocaml
   default: 15

   [threshold_diis]
   type: Threshold
   doc: Threshold on the convergence of the DIIS error vector during a
        Hartree-Fock calculation. If 0. is chosen, the square root of
        thresh_scf will be used.
   interface: ezfio,provider,ocaml
   default: 0.

The syntax obeys the following rules. The name of the parameter in the
|EZFIO| directory is the section header, lowered by the generator.

Required
--------

* ``[<name>]``: the name of the provider in the |IRPF90| files and of
  the parameter in the |EZFIO| directory.
* ``doc:<str>``: the plain text documentation of the parameter.
* ``type:<str>``: a type supported by the generator. The complete list
  of the supported types is obtained by::

      python3 $QP_ROOT/scripts/ezfio_interface/ei_handler.py list_supported_types

  Besides the plain Fortran types (``integer``, ``double precision``,
  ``logical``, ``character*(32)``, ...), the list contains structured
  types such as ``Energy``, ``Det_number``, ``MO_guess`` or
  ``AO_basis_name``.
* ``interface:<str>``: a comma separated list of the interfaces the
  generator builds for the parameter. The recognized values are:

  * ``ezfio``: the parameter is stored in the |EZFIO| directory;
  * ``provider``: an |IRPF90| provider is generated;
  * ``ocaml``: the parameter is bound to the |OCaml| tools
    (in particular :command:`qp edit`).

Optional
--------

* ``default:<str>``: the default value. It is required when
  ``ocaml`` is in the interface list. No default can be set for
  arrays.
* ``size:<str>``: the size of the variable, which is one (scalar) by
  default. The size can be an expression in the |IRPF90| entities, for
  example ``=sum(ao_num)`` or ``(ao_basis.ao_num,3)``. When the
  expression refers to an entity of another module, the module and the
  entity are separated by ``.`` rather than ``_``; for example
  ``(determinants.n_det)``.
* ``ezfio_name:<str>``: the name of the parameter in the |EZFIO|
  directory; by default it is the name of the section, lowered.
* ``ezfio_dir:<str>``: the |EZFIO| group (the subdirectory of the
  |EZFIO| file) in which the parameter is stored; by default it is the
  name of the module, lowered.

Regenerating the interface
--------------------------

When a module is compiled, the build system runs the generator
(:file:`$QP_ROOT/scripts/ezfio_interface/ei_handler.py`) on the
:file:`EZFIO.cfg` file and creates a number of files from it:

* the :file:`ezfio_interface.irp.f` file of the module, which declares
  the generated providers (the file starts with a
  ``DO NOT MODIFY BY HAND`` comment);
* the :file:`<module>.ezfio_interface_config` file in the |EZFIO|
  configuration directory (:file:`$QP_ROOT/external/ezfio/config/`),
  which tells the |EZFIO| library which parameters the group holds;
* the :file:`ocaml/Input_<module>.ml` file and the default value file
  (:file:`$QP_ROOT/data/ezfio_defaults/<module>.ezfio_interface_default`)
  that feed :command:`qp edit`.

These files are never edited by hand: they are regenerated at each
build, so any manual edit is simply lost. The way to change the
declared parameters of a module is to modify its :file:`EZFIO.cfg`
file and to rebuild the module.

The ``.ezfio_config`` files
---------------------------

A module can also declare its |EZFIO| group with a file named
``<group>.ezfio_config`` in the module directory, instead of the
:file:`EZFIO.cfg` syntax. The file lists the parameters of the group
one per line, as in the :file:`bitmasks.ezfio_config` file of the
:ref:`module_bitmask` module::

   bitmasks
     N_int         integer
     bit_kind      integer

