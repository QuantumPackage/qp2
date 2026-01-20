.. _module_generators_full: 
 
.. program:: generators_full 
 
.. default-role:: option 
 
===============
generators_full
===============

Module defining the generator determinants as all the determinants of the
variational space.

This module is intended to be included in the :file:`NEED` file to define
a full set of generators.
 
 
 
Providers 
--------- 
 
.. c:var:: degree_max_generators


    File : :file:`generators_full/generators.irp.f`

    .. code:: fortran

        integer	:: degree_max_generators	


    Max degree of excitation (respect to HF) of the generators

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_det_generators`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`


