.. _module_generators_cas: 
 
.. program:: generators_cas 
 
.. default-role:: option 
 
==============
generators_cas
==============

Generators for the perturbation energy, as the determinants of a
|CAS|.

The generators are the determinants of the wave function that belong
to the |CAS|: the |MOs| of the |CAS| are the ones set as active with
the :ref:`qp_set_mo_class` command. The subroutine
:c:func:`extract_cas` replaces the total wave function by its
normalized projection on the |CAS|.

This module is intended to be included in the :file:`NEED` file to
define the generators as the |CAS| determinants, which is useful to
define post-CAS approaches (see the cassd module for instance).
 
 
 
Providers 
--------- 
 
.. c:var:: double psi_coef_sorted_gen


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

 
.. c:var:: int psi_det_sorted_gen


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

 
.. c:var:: int psi_det_sorted_gen_order


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

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void extract_cas(void)


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

