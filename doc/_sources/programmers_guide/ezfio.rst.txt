=====
EZFIO
=====


EZFIO.cfg
=========

The simplest way to add control parameters in the |EZFIO| directory is to create a
:file:`EZFIO.cfg` file in the module. An example can be found in existing modules
such as :ref:`hartree_fock`::

   [max_dim_diis]
   type: integer
   doc: Maximum size of the |DIIS| extrapolation procedure
   interface: ezfio,provider,ocaml
   default: 15

   [threshold_diis]
   type: Threshold
   doc: Threshold on the convergence of the |DIIS| error vector during a Hartree-Fock calculation. If 0. is chosen, the square root of thresh_scf will be used.
   interface: ezfio,provider,ocaml
   default: 0.

   [thresh_scf]
   type: Threshold
   doc: Threshold on the convergence of the Hartree Fock energy.
   interface: ezfio,provider,ocaml
   default: 1.e-10



The syntax obeys the following rules:

Required
--------

.. option::  [<provider_name>]

   The name of the provider in irp.f90 and in the EZFIO lib

.. option:: doc:<str>

   The plain text documentation

.. option:: type:<str>

   A type supported by the |OCaml| modules. The complete list of supported
   types can be obtained by::

       ei_handler.py list_supported_types


.. option:: interface:<str>

   The interface is a list of strings sepeared by ","  which can contain :

   - ``ezfio`` : to build the |EZFIO| API
   - ``provider`` : to build the corresponding providers
   - ``ocaml`` : to build the corresponding bindings in |OCaml|

If an ``EZFIO.cfg`` file is used, the compilation of the module will generate
the ``ezfio_interface.irp.f`` file which contains the generated providers.
This file should not be added to the repository

Optional
--------

.. option:: default:<str>

   The default value needed if ``ocaml`` is in interface list.
   No default can be set for arrays.

.. option:: size:<str>

   The size of the variable, which is one by default (scalar).

   Examples : ``1``; ``=sum(ao_num)``; ``(ao_basis.ao_num,3)``

   .. warning::

      The module and the value are separed by a ``.`` not a ``_``.
      For example ``(determinants.n_det)``

.. option:: ezfio_name:<str>

   The name in the |EZFIO| API (by default is ``<provider_name>``)


\*.ezfio_config
===============

It is possible to directly add to the current module |EZFIO| configuration
files, named with the ``.ezfio_config`` suffix. An example is in the
:ref:`bitmask` module. 

.. code:: text

   bitmasks
      N_int         integer
      bit_kind      integer
      N_mask_gen    integer
      generators    integer*8 (bitmasks_N_int*bitmasks_bit_kind/8,2,6,bitmasks_N_mask_gen)
      N_mask_cas    integer
      cas           integer*8 (bitmasks_N_int*bitmasks_bit_kind/8,2,bitmasks_N_mask_cas)







