.. _module_dav_general_mat: 
 
.. program:: dav_general_mat 
 
.. default-role:: option 
 
===============
dav_general_mat
===============

Davidson diagonalization of a general symmetric matrix.

Unlike the Davidson diagonalization of the Hamiltonian (see the
:ref:`module_davidson` module), this module diagonalizes any symmetric
matrix, with two ways to provide the matrix-vector product:

* :c:func:`davidson_general` takes the matrix directly
* :c:func:`davidson_general_ext_rout` takes a subroutine that
  computes :math:`v = H u`

Both have dressed variants, :c:func:`davidson_general_ext_rout_dressed`
and :c:func:`davidson_general_ext_rout_diag_dressed`, that add the
dressing operator of the :ref:`module_dressing` module, as used in
the |CIPSI| iteration, and non symmetric variants that do not assume
the symmetry of the matrix.

The program :ref:`test_dav` is a small example of the use of these
routines (see :file:`test_dav.irp.f`).
 
 
 
Programs 
-------- 
 
 * :ref:`test_dav` 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: dav_double_dressed:


    File : :file:`dav_double_dress_ext_rout.irp.f`

    .. code:: fortran

        subroutine dav_double_dressed(u_in,H_jj,Dress_jj,Dressing_vec,idx_dress,energies,sze,N_st,N_st_diag,converged,hcalc)


    Generic Davidson diagonalization with TWO DRESSING VECTORS
    
    Dress_jj : DIAGONAL DRESSING of the Hamiltonian
    
    Dressing_vec : COLUMN / LINE DRESSING VECTOR
    
    idx_dress : position of the basis function used to use the Dressing_vec (usually the largest coeff)
    
    H_jj : specific diagonal H matrix elements to diagonalize de Davidson
    
    u_in : guess coefficients on the various states. Overwritten on exit
    
    sze : leftmost dimension of u_in
    
    sze : Number of determinants
    
    N_st : Number of eigenstates
    
    N_st_diag : Number of states in which H is diagonalized. Assumed > sze
    
    Initial guess vectors are not necessarily orthonormal
    
    hcalc subroutine to compute W = H U (see routine hcalc_template for template of input/output)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`disk_based_davidson`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`qp_max_mem`
       * :c:data:`state_following`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_davidson_pt2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`dgemm`
       * :c:func:`dressing_diag_uv`
       * :c:func:`dsygv`
       * :c:func:`hcalc`
       * :c:func:`normalize`
       * :c:func:`nullify_small_elements`
       * :c:func:`random_number`
       * :c:func:`resident_memory`
       * :c:func:`write_double`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:macro:: davidson_general:


    File : :file:`dav_general.irp.f`

    .. code:: fortran

        subroutine davidson_general(u_in,H_jj,energies,dim_in,sze,N_st,N_st_diag_in,converged,h_mat)


    Davidson diagonalization with specific diagonal elements of the H matrix
    
    H_jj : specific diagonal H matrix elements to diagonalize de Davidson
    
    u_in : guess coefficients on the various states. Overwritten on exit
    
    dim_in : leftmost dimension of u_in
    
    sze : Number of determinants
    
    N_st : Number of eigenstates
    
    N_st_diag_in : Number of states in which H is diagonalized. Assumed > sze
    
    Initial guess vectors are not necessarily orthonormal

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`disk_based_davidson`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`qp_max_mem`
       * :c:data:`state_following`
       * :c:data:`threshold_davidson`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`routine`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`dgemm`
       * :c:func:`hpsi`
       * :c:func:`lapack_diag`
       * :c:func:`normalize`
       * :c:func:`ortho_qr`
       * :c:func:`random_number`
       * :c:func:`resident_memory`
       * :c:func:`write_double`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:macro:: davidson_general_diag_dressed_ext_rout_nonsym_b1space:


    File : :file:`dav_diag_dressed_ext_rout_nonsym_B1space.irp.f`

    .. code:: fortran

        subroutine davidson_general_diag_dressed_ext_rout_nonsym_b1space(u_in, H_jj, Dress_jj,energies, sze, N_st, N_st_diag_in, converged, hcalc)


    Generic modified-Davidson diagonalization
    
    H_jj : specific diagonal H matrix elements to diagonalize de Davidson
    
    u_in : guess coefficients on the various states. Overwritten on exit by right eigenvectors
    
    sze : Number of determinants
    
    N_st : Number of eigenstates
    
    N_st_diag_in : Number of states in which H is diagonalized. Assumed > N_st
    
    Initial guess vectors are not necessarily orthonormal
    
    hcalc subroutine to compute W = H U (see routine hcalc_template for template of input/output)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`disk_based_davidson`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`qp_max_mem`
       * :c:data:`threshold_nonsym_davidson`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`dgemm`
       * :c:func:`diag_nonsym_right`
       * :c:func:`dress_calc`
       * :c:func:`hcalc`
       * :c:func:`normalize`
       * :c:func:`ortho_qr`
       * :c:func:`random_number`
       * :c:func:`resident_memory`
       * :c:func:`write_double`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:macro:: davidson_general_ext_rout:


    File : :file:`dav_ext_rout.irp.f`

    .. code:: fortran

        subroutine davidson_general_ext_rout(u_in,H_jj,energies,sze,N_st,N_st_diag_in,converged,hcalc)


    Generic Davidson diagonalization
    
    H_jj : specific diagonal H matrix elements to diagonalize de Davidson
    
    u_in : guess coefficients on the various states. Overwritten on exit
    
    sze : leftmost dimension of u_in
    
    sze : Number of determinants
    
    N_st : Number of eigenstates
    
    N_st_diag_in : Number of states in which H is diagonalized. Assumed > sze
    
    Initial guess vectors are not necessarily orthonormal
    
    hcalc subroutine to compute W = H U (see routine hcalc_template for template of input/output)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`disk_based_davidson`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`qp_max_mem`
       * :c:data:`state_following`
       * :c:data:`threshold_davidson`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`routine`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`dgemm`
       * :c:func:`hcalc`
       * :c:func:`lapack_diag`
       * :c:func:`normalize`
       * :c:func:`ortho_qr`
       * :c:func:`random_number`
       * :c:func:`resident_memory`
       * :c:func:`write_double`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:macro:: davidson_general_ext_rout_diag_dressed:


    File : :file:`dav_diag_dressed_ext_rout.irp.f`

    .. code:: fortran

        subroutine davidson_general_ext_rout_diag_dressed(u_in,H_jj,Dress_jj,energies,sze,N_st,N_st_diag_in,converged,hcalc)


    Generic Davidson diagonalization with ONE DIAGONAL DRESSING OPERATOR
    
    Dress_jj : DIAGONAL DRESSING of the Hamiltonian
    
    H_jj : specific diagonal H matrix elements to diagonalize de Davidson
    
    u_in : guess coefficients on the various states. Overwritten on exit
    
    sze : leftmost dimension of u_in
    
    sze : Number of determinants
    
    N_st : Number of eigenstates
    
    N_st_diag_in : Number of states in which H is diagonalized. Assumed > sze
    
    Initial guess vectors are not necessarily orthonormal
    
    hcalc subroutine to compute W = H U (see routine hcalc_template for template of input/output)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`disk_based_davidson`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`qp_max_mem`
       * :c:data:`state_following`
       * :c:data:`threshold_davidson`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`dgemm`
       * :c:func:`dressing_diag_uv`
       * :c:func:`hcalc`
       * :c:func:`lapack_diag`
       * :c:func:`normalize`
       * :c:func:`ortho_qr`
       * :c:func:`random_number`
       * :c:func:`resident_memory`
       * :c:func:`write_double`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:macro:: davidson_general_ext_rout_dressed:


    File : :file:`dav_dressed_ext_rout.irp.f`

    .. code:: fortran

        subroutine davidson_general_ext_rout_dressed(u_in,H_jj,energies,sze,N_st,N_st_diag,dressing_state,dressing_vec,idress,converged,hcalc)


    Davidson diagonalization.
    
    u_in : guess coefficients on the various states. Overwritten
      on exit
    
    sze : leftmost dimension of u_in
    
    sze : Number of determinants
    
    N_st : Number of eigenstates
    
    Initial guess vectors are not necessarily orthonormal

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`disk_based_davidson`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`qp_max_mem`
       * :c:data:`state_following`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_davidson_pt2`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_dav_dress`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`dgemm`
       * :c:func:`dsygv`
       * :c:func:`hcalc`
       * :c:func:`normalize`
       * :c:func:`nullify_small_elements`
       * :c:func:`ortho_qr`
       * :c:func:`random_number`
       * :c:func:`resident_memory`
       * :c:func:`write_double`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:macro:: davidson_general_ext_rout_nonsym_b1space:


    File : :file:`dav_ext_rout_nonsym_B1space.irp.f`

    .. code:: fortran

        subroutine davidson_general_ext_rout_nonsym_b1space(u_in, H_jj, energies, sze, N_st, N_st_diag_in, converged, hcalc)


    Generic modified-Davidson diagonalization
    
    H_jj : specific diagonal H matrix elements to diagonalize de Davidson
    
    u_in : guess coefficients on the various states. Overwritten on exit by right eigenvectors
    
    sze : Number of determinants
    
    N_st : Number of eigenstates
    
    N_st_diag_in : Number of states in which H is diagonalized. Assumed > N_st
    
    Initial guess vectors are not necessarily orthonormal
    
    hcalc subroutine to compute W = H U (see routine hcalc_template for template of input/output)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`davidson_sze_max`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`qp_max_mem`
       * :c:data:`threshold_nonsym_davidson`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`check_mem`
       * :c:func:`dgemm`
       * :c:func:`diag_nonsym_right`
       * :c:func:`hcalc`
       * :c:func:`normalize`
       * :c:func:`ortho_qr`
       * :c:func:`random_number`
       * :c:func:`resident_memory`
       * :c:func:`write_double`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:macro:: dress_calc:


    File : :file:`dav_diag_dressed_ext_rout_nonsym_B1space.irp.f`

    .. code:: fortran

        subroutine dress_calc(v,dress,u,N_st,sze)


    Routine that computed the action of the diagonal dressing dress
    
    WARNING :: v is not initialiazed !!!

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_general_diag_dressed_ext_rout_nonsym_b1space`

 
.. c:macro:: dressing_diag_uv:


    File : :file:`dav_double_dress_ext_rout.irp.f`

    .. code:: fortran

        subroutine dressing_diag_uv(v,u,dress_diag,N_st,sze)


    Routine that computes the diagonal part of the dressing
    
    v(i) =! v(i) +( u(i) * dress_diag(i))
    
    !!!!!!!! WARNING !!!!!!!! the vector v is not initialized
    
    !!!!!!!! SO MAKE SURE THERE ARE SOME MEANINGFUL VALUES IN THERE

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`dav_double_dressed`
       * :c:func:`davidson_general_ext_rout_diag_dressed`

 
.. c:macro:: hcalc_template:


    File : :file:`dav_ext_rout.irp.f`

    .. code:: fortran

        subroutine hcalc_template(v,u,N_st,sze)


    Template of routine for the application of H
    
    Here, it is done with the Hamiltonian matrix
    
    on the set of determinants of psi_det
    
    Computes $v = H | u \rangle$
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_all_dets`
       * :c:data:`nuclear_repulsion`

 
.. c:macro:: hpsi:


    File : :file:`dav_general.irp.f`

    .. code:: fortran

        subroutine hpsi(v,u,N_st,sze,h_mat)


    Computes $v = H | u \rangle$ and

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_general`

 
.. c:macro:: routine:


    File : :file:`test_dav.irp.f`

    .. code:: fortran

        subroutine routine



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_all_dets`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nthreads_davidson`
       * :c:data:`nuclear_repulsion`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_dav`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_general`
       * :c:func:`davidson_general_ext_rout`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:macro:: test_dav_dress:


    File : :file:`test_dav.irp.f`

    .. code:: fortran

        subroutine test_dav_dress



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`h_matrix_all_dets`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nuclear_repulsion`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_dav`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_general_ext_rout_dressed`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

