.. _module_selectors_full: 
 
.. program:: selectors_full 
 
.. default-role:: option 
 
==============
selectors_full
==============

All the determinants are possible selectors. Only the largest contributions are kept, where
a threshold is applied to the squared norm of the wave function, with the :option:`determinants
threshold_selectors` flag.
 
 
 
Providers 
--------- 
 
.. c:var:: n_det_selectors


    File : :file:`selectors_full/selectors.irp.f`

    .. code:: fortran

        integer	:: n_det_selectors	


    For Single reference wave functions, the number of selectors is 1 : the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`output_wall_time_0`
       * :c:data:`psi_det_sorted`
       * :c:data:`threshold_selectors`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_hf_selector`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_coef_transp`
       * :c:data:`psi_selectors_diag_h_mat`

 
.. c:var:: psi_selectors


    File : :file:`selectors_full/selectors.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_selectors	(N_int,2,psi_selectors_size)
        double precision, allocatable	:: psi_selectors_coef	(psi_selectors_size,N_states)


    Determinants on which we apply <i|H|psi> for perturbation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_selectors_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_hf_selector`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`psi_selectors_coef_transp`
       * :c:data:`psi_selectors_diag_h_mat`

 
.. c:var:: psi_selectors_coef


    File : :file:`selectors_full/selectors.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_selectors	(N_int,2,psi_selectors_size)
        double precision, allocatable	:: psi_selectors_coef	(psi_selectors_size,N_states)


    Determinants on which we apply <i|H|psi> for perturbation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_selectors_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_hf_selector`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`psi_selectors_coef_transp`
       * :c:data:`psi_selectors_diag_h_mat`

 
.. c:var:: threshold_selectors


    File : :file:`selectors_full/selectors.irp.f`

    .. code:: fortran

        double precision	:: threshold_selectors	


    Thresholds on selectors (fraction of the square of the norm)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`threshold_generators`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`

