.. _module_generators_cas: 
 
.. program:: generators_cas 
 
.. default-role:: option 
 
==============
generators_cas
==============

Module defining the generator determinants as those belonging to a |CAS|.
The |MOs| belonging to the |CAS| are those which were set as active with
the :ref:`qp_set_mo_class` command.

This module is intended to be included in the :file:`NEED` file to define
the generators as the |CAS| determinants, which can be useful to define post-CAS approaches (see cassd module for instance).


 
 
 
Providers 
--------- 
 
.. c:var:: n_det_generators


    File : :file:`generators_cas/generators.irp.f`

    .. code:: fortran

        integer	:: n_det_generators	


    Number of generator detetrminants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`output_wall_time_0`
       * :c:data:`psi_det_sorted`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`global_selection_buffer`
       * :c:data:`n_det_selectors`
       * :c:data:`psi_det_generators`
       * :c:data:`pt2_f`
       * :c:data:`pt2_j`
       * :c:data:`pt2_n_tasks`
       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_u`
       * :c:data:`pt2_w`

 
.. c:var:: psi_coef_generators


    File : :file:`generators_cas/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)
        integer(bit_kind), allocatable	:: psi_det_sorted_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: psi_coef_sorted_gen


    File : :file:`generators_cas/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)
        integer(bit_kind), allocatable	:: psi_det_sorted_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: psi_det_generators


    File : :file:`generators_cas/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)
        integer(bit_kind), allocatable	:: psi_det_sorted_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: psi_det_sorted_gen


    File : :file:`generators_cas/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)
        integer(bit_kind), allocatable	:: psi_det_sorted_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: psi_det_sorted_gen_order


    File : :file:`generators_cas/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)
        integer(bit_kind), allocatable	:: psi_det_sorted_gen	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_sorted_gen	(psi_det_size,N_states)
        integer, allocatable	:: psi_det_sorted_gen_order	(psi_det_size)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`act_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt2_n_teeth`
       * :c:data:`pt2_w`

 
.. c:var:: select_max


    File : :file:`generators_cas/generators.irp.f`

    .. code:: fortran

        double precision, allocatable	:: select_max	(size_select_max)


    Memo to skip useless selectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`size_select_max`


 
.. c:var:: size_select_max


    File : :file:`generators_cas/generators.irp.f`

    .. code:: fortran

        integer	:: size_select_max	


    Size of the select_max array

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`select_max`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: extract_cas:


    File : :file:`generators_cas/extract_cas.irp.f`

    .. code:: fortran

        subroutine extract_cas


    Replaces the total wave function by the normalized projection on the CAS.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_generators`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

