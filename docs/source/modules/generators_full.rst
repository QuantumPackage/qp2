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


 
.. c:var:: n_det_generators


    File : :file:`generators_full/generators.irp.f`

    .. code:: fortran

        integer	:: n_det_generators	


    For Single reference wave functions, the number of generators is 1 : the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`output_wall_time_0`
       * :c:data:`psi_det_sorted`
       * :c:data:`threshold_generators`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_generators`
       * :c:data:`global_selection_buffer`
       * :c:data:`n_det_selectors`
       * :c:data:`pt2_f`
       * :c:data:`pt2_j`
       * :c:data:`pt2_n_tasks`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_u`
       * :c:data:`pt2_w`

 
.. c:var:: psi_coef_generators


    File : :file:`generators_full/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_generators`

 
.. c:var:: psi_coef_sorted_gen


    File : :file:`generators_full/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: psi_det_generators


    File : :file:`generators_full/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`degree_max_generators`

 
.. c:var:: psi_det_sorted_gen


    File : :file:`generators_full/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: psi_det_sorted_gen_order


    File : :file:`generators_full/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_sorted_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: select_max


    File : :file:`generators_full/generators.irp.f`

    .. code:: fortran

        double precision, allocatable	:: select_max	(size_select_max)


    Memo to skip useless selectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`size_select_max`


 
.. c:var:: size_select_max


    File : :file:`generators_full/generators.irp.f`

    .. code:: fortran

        integer	:: size_select_max	


    Size of the select_max array

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`select_max`

