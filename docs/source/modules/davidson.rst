.. _module_davidson: 
 
.. program:: davidson 
 
.. default-role:: option 
 
========
davidson
========

Abstract module for Davidson's diagonalization.
It contains everything required for the Davidson algorithm, dressed or
not. If a dressing is used, the dressing column should be defined and
the :ref:`module_davidson_dressed` module should be used. If no dressing
is required, the :ref:`module_davidson` module should be used, and it
has a default zero dressing vector.

The important providers for that module are:

#. :c:data:`psi_energy` which is the expectation value over the wave
   function (:c:data:`psi_det`, :c:data:`psi_coef`) of the Hamiltonian,
   dressed or not. It uses the general subroutine :c:func:`u_0_H_u_0`.

#. :c:data:`psi_energy_two_e` which is the expectation value over the
   wave function (:c:data:`psi_det`, :c:data:`psi_coef`) of the standard
   two-electron Coulomb operator. It uses the general routine
   :c:func:`u_0_H_u_0_two_e`.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: threshold_davidson
 
    Thresholds of Davidson's algorithm
 
    Default: 1.e-10
 
.. option:: n_states_diag
 
    Number of states to consider during the Davdison diagonalization
 
    Default: 4
 
.. option:: davidson_sze_max
 
    Number of micro-iterations before re-contracting
 
    Default: 15
 
.. option:: state_following
 
    If |true|, the states are re-ordered to match the input states
 
    Default: False
 
.. option:: disk_based_davidson
 
    If |true|, a memory-mapped file may be used to store the W and S2 vectors if not enough RAM is available
 
    Default: True
 
.. option:: distributed_davidson
 
    If |true|, use the distributed algorithm
 
    Default: True
 
.. option:: only_expected_s2
 
    If |true|, use filter out all vectors with bad |S^2| values
 
    Default: True
 
.. option:: n_det_max_full
 
    Maximum number of determinants where |H| is fully diagonalized
 
    Default: 1000
 
 
Providers 
--------- 
 
.. c:var:: ci_eigenvectors


    File : :file:`davidson/diagonalize_ci.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ci_electronic_energy	(N_states_diag)
        double precision, allocatable	:: ci_eigenvectors	(N_det,N_states_diag)
        double precision, allocatable	:: ci_s2	(N_states_diag)


    Eigenvectors/values of the |CI| matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`diag_algorithm`
       * :c:data:`dressing_column_h`
       * :c:data:`expected_s2`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nthreads_davidson`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`s2_eig`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`threshold_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy`

 
.. c:var:: ci_electronic_energy


    File : :file:`davidson/diagonalize_ci.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ci_electronic_energy	(N_states_diag)
        double precision, allocatable	:: ci_eigenvectors	(N_det,N_states_diag)
        double precision, allocatable	:: ci_s2	(N_states_diag)


    Eigenvectors/values of the |CI| matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`diag_algorithm`
       * :c:data:`dressing_column_h`
       * :c:data:`expected_s2`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nthreads_davidson`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`s2_eig`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`threshold_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy`

 
.. c:var:: ci_energy


    File : :file:`davidson/diagonalize_ci.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ci_energy	(N_states_diag)


    :c:data:`n_states` lowest eigenvalues of the |CI| matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nuclear_repulsion`
       * :c:data:`output_wall_time_0`


 
.. c:var:: ci_s2


    File : :file:`davidson/diagonalize_ci.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ci_electronic_energy	(N_states_diag)
        double precision, allocatable	:: ci_eigenvectors	(N_det,N_states_diag)
        double precision, allocatable	:: ci_s2	(N_states_diag)


    Eigenvectors/values of the |CI| matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`diag_algorithm`
       * :c:data:`dressing_column_h`
       * :c:data:`expected_s2`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`nthreads_davidson`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`s2_eig`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`threshold_davidson`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy`

 
.. c:var:: diag_algorithm


    File : :file:`davidson/diagonalization_hs2_dressed.irp.f`

    .. code:: fortran

        character*(64)	:: diag_algorithm	


    Diagonalization algorithm (Davidson or Lapack)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`

 
.. c:var:: dressed_column_idx


    File : :file:`davidson/diagonalization_hs2_dressed.irp.f`

    .. code:: fortran

        integer, allocatable	:: dressed_column_idx	(N_states)


    Index of the dressed columns

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`


 
.. c:var:: n_states_diag


    File : :file:`davidson/input.irp.f`

    .. code:: fortran

        integer	:: n_states_diag	


    Number of states to consider during the Davdison diagonalization

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`n_states`
       * :c:data:`output_wall_time_0`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
       * :c:data:`psi_energy`

 
.. c:var:: nthreads_davidson


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        integer	:: nthreads_davidson	


    Number of threads for Davidson

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`nproc`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`

 
.. c:var:: psi_energy


    File : :file:`davidson/u0_h_u0.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_energy	(N_states)
        double precision, allocatable	:: psi_s2	(N_states)


    psi_energy(i) = :math:`\langle \Psi_i | H | \Psi_i \rangle` 
    
    psi_s2(i) = :math:`\langle \Psi_i | S^2 | \Psi_i \rangle` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`distributed_davidson`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`s2_matrix_all_dets`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_energy_with_nucl_rep`

 
.. c:var:: psi_energy_two_e


    File : :file:`davidson/u0_wee_u0.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_energy_two_e	(N_states)


    Energy of the current wave function

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`


 
.. c:var:: psi_energy_with_nucl_rep


    File : :file:`davidson/u0_h_u0.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_energy_with_nucl_rep	(N_states)


    Energy of the wave function with the nuclear repulsion energy.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_energy`


 
.. c:var:: psi_s2


    File : :file:`davidson/u0_h_u0.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_energy	(N_states)
        double precision, allocatable	:: psi_s2	(N_states)


    psi_energy(i) = :math:`\langle \Psi_i | H | \Psi_i \rangle` 
    
    psi_s2(i) = :math:`\langle \Psi_i | S^2 | \Psi_i \rangle` 

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`distributed_davidson`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`n_det`
       * :c:data:`n_det_max_full`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`n_states_diag`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`s2_matrix_all_dets`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_energy_with_nucl_rep`

 
.. c:var:: threshold_davidson_pt2


    File : :file:`davidson/diagonalization_hs2_dressed.irp.f`

    .. code:: fortran

        double precision	:: threshold_davidson_pt2	


    Threshold of Davidson's algorithm, using PT2 as a guide

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`threshold_davidson`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: davidson_collector:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        subroutine davidson_collector(zmq_to_qp_run_socket, zmq_socket_pull, v0, s0, sze, N_st)


    Routine collecting the results of the workers in Davidson's algorithm.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_zmq`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_pull_results`

 
.. c:function:: davidson_diag_hjj_sjj:


    File : :file:`davidson/diagonalization_hs2_dressed.irp.f`

    .. code:: fortran

        subroutine davidson_diag_hjj_sjj(dets_in,u_in,H_jj,s2_out,energies,dim_in,sze,N_st,N_st_diag_in,Nint,dressing_state,converged)


    Davidson diagonalization with specific diagonal elements of the H matrix
    
    H_jj : specific diagonal H matrix elements to diagonalize de Davidson
    
    S2_out : Output : s^2
    
    dets_in : bitmasks corresponding to determinants
    
    u_in : guess coefficients on the various states. Overwritten
      on exit
    
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
       * :c:data:`distributed_davidson`
       * :c:data:`dressed_column_idx`
       * :c:data:`dressing_column_h`
       * :c:data:`expected_s2`
       * :c:data:`ezfio_work_dir`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nproc`
       * :c:data:`nthreads_davidson`
       * :c:data:`nuclear_repulsion`
       * :c:data:`only_expected_s2`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`qp_max_mem`
       * :c:data:`s2_eig`
       * :c:data:`state_following`
       * :c:data:`threshold_davidson`
       * :c:data:`threshold_davidson_from_pt2`
       * :c:data:`threshold_davidson_pt2`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_diag_hs2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`c_f_pointer`
       * :c:func:`check_mem`
       * :c:func:`dgemm`
       * :c:func:`dswap`
       * :c:func:`h_s2_u_0_nstates_openmp`
       * :c:func:`h_s2_u_0_nstates_zmq`
       * :c:func:`lapack_diag`
       * :c:func:`mmap`
       * :c:func:`munmap`
       * :c:func:`normalize`
       * :c:func:`nullify_small_elements`
       * :c:func:`ortho_qr`
       * :c:func:`ortho_qr_unblocked`
       * :c:func:`random_number`
       * :c:func:`resident_memory`
       * :c:func:`sgemm`
       * :c:func:`write_double`
       * :c:func:`write_int`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:function:: davidson_diag_hs2:


    File : :file:`davidson/diagonalization_hs2_dressed.irp.f`

    .. code:: fortran

        subroutine davidson_diag_hs2(dets_in,u_in,s2_out,dim_in,energies,sze,N_st,N_st_diag,Nint,dressing_state,converged)


    Davidson diagonalization.
    
    dets_in : bitmasks corresponding to determinants
    
    u_in : guess coefficients on the various states. Overwritten
      on exit
    
    dim_in : leftmost dimension of u_in
    
    sze : Number of determinants
    
    N_st : Number of eigenstates
    
    Initial guess vectors are not necessarily orthonormal

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dressing_column_h`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_diag_hjj_sjj`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:function:: davidson_pull_results:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        subroutine davidson_pull_results(zmq_socket_pull, v_t, s_t, imin, imax, task_id)


    Pull the results of $H | U \rangle$ on the master.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states_diag`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_collector`

 
.. c:function:: davidson_push_results:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        subroutine davidson_push_results(zmq_socket_push, v_t, s_t, imin, imax, task_id)


    Push the results of $H | U \rangle$ from a worker to the master.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states_diag`

 
.. c:function:: davidson_push_results_async_recv:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        subroutine davidson_push_results_async_recv(zmq_socket_push,sending)


    Push the results of $H | U \rangle$ from a worker to the master.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_slave_work`

 
.. c:function:: davidson_push_results_async_send:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        subroutine davidson_push_results_async_send(zmq_socket_push, v_t, s_t, imin, imax, task_id,sending)


    Push the results of $H | U \rangle$ from a worker to the master.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states_diag`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_slave_work`

 
.. c:function:: davidson_run_slave:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        subroutine davidson_run_slave(thread,iproc)


    Slave routine for Davidson's diagonalization.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states_diag`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_slave_inproc`
       * :c:func:`davidson_slave_tcp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_slave_work`
       * :c:func:`end_zmq_push_socket`
       * :c:func:`end_zmq_to_qp_run_socket`
       * :c:func:`sleep`

 
.. c:function:: davidson_slave_inproc:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        subroutine davidson_slave_inproc(i)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_zmq`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_run_slave`

 
.. c:function:: davidson_slave_tcp:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        subroutine davidson_slave_tcp(i)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_run_slave`

 
.. c:function:: davidson_slave_work:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        subroutine davidson_slave_work(zmq_to_qp_run_socket, zmq_socket_push, N_st, sze, worker_id)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_initialized`
       * :c:data:`mpi_rank`
       * :c:data:`n_det`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`ref_bitmask_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_run_slave`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_push_results_async_recv`
       * :c:func:`davidson_push_results_async_send`
       * :c:func:`h_s2_u_0_nstates_openmp_work`

 
.. c:function:: diagonalize_ci:


    File : :file:`davidson/diagonalize_ci.irp.f`

    .. code:: fortran

        subroutine diagonalize_CI


    Replace the coefficients of the |CI| states by the coefficients of the
    eigenstates of the |CI| matrix.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`psi_coef`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy`

 
.. c:function:: h_s2_u_0_nstates_openmp:


    File : :file:`davidson/u0_h_u0.irp.f`

    .. code:: fortran

        subroutine H_S2_u_0_nstates_openmp(v_0,s_0,u_0,N_st,sze)


    Computes $v_0 = H | u_0\rangle$ and $s_0 = S^2  | u_0\rangle$.
    
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep are used in ZMQ parallelization.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`u_0_h_u_0`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dset_order`
       * :c:func:`dtranspose`
       * :c:func:`h_s2_u_0_nstates_openmp_work`

 
.. c:function:: h_s2_u_0_nstates_openmp_work:


    File : :file:`davidson/u0_h_u0.irp.f`

    .. code:: fortran

        subroutine H_S2_u_0_nstates_openmp_work(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t\rangle$ and $s_t = S^2  | u_t\rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_slave_work`
       * :c:func:`h_s2_u_0_nstates_openmp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_nstates_openmp_work_n_int`

 
.. c:function:: h_s2_u_0_nstates_openmp_work_1:


    File : :file:`davidson/u0_h_u0.irp.f_template_708`

    .. code:: fortran

        subroutine H_S2_u_0_nstates_openmp_work_1(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t \rangle$ and $s_t = S^2 | u_t\rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_1`
       * :c:func:`get_all_spin_singles_and_doubles_1`
       * :c:func:`get_s2`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`i_h_j_single_spin`

 
.. c:function:: h_s2_u_0_nstates_openmp_work_2:


    File : :file:`davidson/u0_h_u0.irp.f_template_708`

    .. code:: fortran

        subroutine H_S2_u_0_nstates_openmp_work_2(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t \rangle$ and $s_t = S^2 | u_t\rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_2`
       * :c:func:`get_all_spin_singles_and_doubles_2`
       * :c:func:`get_s2`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`i_h_j_single_spin`

 
.. c:function:: h_s2_u_0_nstates_openmp_work_3:


    File : :file:`davidson/u0_h_u0.irp.f_template_708`

    .. code:: fortran

        subroutine H_S2_u_0_nstates_openmp_work_3(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t \rangle$ and $s_t = S^2 | u_t\rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_3`
       * :c:func:`get_all_spin_singles_and_doubles_3`
       * :c:func:`get_s2`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`i_h_j_single_spin`

 
.. c:function:: h_s2_u_0_nstates_openmp_work_4:


    File : :file:`davidson/u0_h_u0.irp.f_template_708`

    .. code:: fortran

        subroutine H_S2_u_0_nstates_openmp_work_4(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t \rangle$ and $s_t = S^2 | u_t\rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_4`
       * :c:func:`get_all_spin_singles_and_doubles_4`
       * :c:func:`get_s2`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`i_h_j_single_spin`

 
.. c:function:: h_s2_u_0_nstates_openmp_work_n_int:


    File : :file:`davidson/u0_h_u0.irp.f_template_708`

    .. code:: fortran

        subroutine H_S2_u_0_nstates_openmp_work_N_int(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t \rangle$ and $s_t = S^2 | u_t\rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_and_doubles_n_int`
       * :c:func:`get_all_spin_singles_n_int`
       * :c:func:`get_s2`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`i_h_j_single_spin`

 
.. c:function:: h_s2_u_0_nstates_zmq:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        subroutine H_S2_u_0_nstates_zmq(v_0,s_0,u_0,N_st,sze)


    Computes $v_0 = H | u_0\rangle$ and $s_0 = S^2  | u_0\rangle$
    
    n : number of determinants
    
    H_jj : array of $\langle j | H | j \rangle$
    
    S2_jj : array of $\langle j | S^2 | j \rangle$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_initialized`
       * :c:data:`n_det`
       * :c:data:`n_states_diag`
       * :c:data:`nproc`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`ref_bitmask_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_diag_hjj_sjj`
       * :c:func:`u_0_h_u_0`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_collector`
       * :c:func:`davidson_slave_inproc`
       * :c:func:`dset_order`
       * :c:func:`dtranspose`
       * :c:func:`end_parallel_job`
       * :c:func:`new_parallel_job`
       * :c:func:`omp_set_nested`

 
.. c:function:: h_s2_u_0_two_e_nstates_openmp:


    File : :file:`davidson/u0_wee_u0.irp.f`

    .. code:: fortran

        subroutine H_S2_u_0_two_e_nstates_openmp(v_0,s_0,u_0,N_st,sze)


    Computes $v_0 = H | u_0\rangle$ and $s_0 = S^2 | u_0\rangle$
    
    Assumes that the determinants are in psi_det
    
    istart, iend, ishift, istep are used in ZMQ parallelization.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`psi_bilinear_matrix_order_reverse`
       * :c:data:`psi_bilinear_matrix_values`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`u_0_h_u_0_two_e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dset_order`
       * :c:func:`dtranspose`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work`

 
.. c:function:: h_s2_u_0_two_e_nstates_openmp_work:


    File : :file:`davidson/u0_wee_u0.irp.f`

    .. code:: fortran

        subroutine H_S2_u_0_two_e_nstates_openmp_work(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t\rangle$ and $s_t = S^2 | u_t\rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_two_e_nstates_openmp`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_1`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_2`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_3`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_4`
       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work_n_int`

 
.. c:function:: h_s2_u_0_two_e_nstates_openmp_work_1:


    File : :file:`davidson/u0_wee_u0.irp.f_template_457`

    .. code:: fortran

        subroutine H_S2_u_0_two_e_nstates_openmp_work_1(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t \rangle$ and $s_t = S^2  | u_t \rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_1`
       * :c:func:`get_all_spin_singles_and_doubles_1`
       * :c:func:`get_s2`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`i_wee_j_single`

 
.. c:function:: h_s2_u_0_two_e_nstates_openmp_work_2:


    File : :file:`davidson/u0_wee_u0.irp.f_template_457`

    .. code:: fortran

        subroutine H_S2_u_0_two_e_nstates_openmp_work_2(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t \rangle$ and $s_t = S^2  | u_t \rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_2`
       * :c:func:`get_all_spin_singles_and_doubles_2`
       * :c:func:`get_s2`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`i_wee_j_single`

 
.. c:function:: h_s2_u_0_two_e_nstates_openmp_work_3:


    File : :file:`davidson/u0_wee_u0.irp.f_template_457`

    .. code:: fortran

        subroutine H_S2_u_0_two_e_nstates_openmp_work_3(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t \rangle$ and $s_t = S^2  | u_t \rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_3`
       * :c:func:`get_all_spin_singles_and_doubles_3`
       * :c:func:`get_s2`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`i_wee_j_single`

 
.. c:function:: h_s2_u_0_two_e_nstates_openmp_work_4:


    File : :file:`davidson/u0_wee_u0.irp.f_template_457`

    .. code:: fortran

        subroutine H_S2_u_0_two_e_nstates_openmp_work_4(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t \rangle$ and $s_t = S^2  | u_t \rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_4`
       * :c:func:`get_all_spin_singles_and_doubles_4`
       * :c:func:`get_s2`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`i_wee_j_single`

 
.. c:function:: h_s2_u_0_two_e_nstates_openmp_work_n_int:


    File : :file:`davidson/u0_wee_u0.irp.f_template_457`

    .. code:: fortran

        subroutine H_S2_u_0_two_e_nstates_openmp_work_N_int(v_t,s_t,u_t,N_st,sze,istart,iend,ishift,istep)


    Computes $v_t = H | u_t \rangle$ and $s_t = S^2  | u_t \rangle$
    
    Default should be 1,N_det,0,1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`nthreads_davidson`
       * :c:data:`psi_bilinear_matrix_columns_loc`
       * :c:data:`psi_bilinear_matrix_order_transp_reverse`
       * :c:data:`psi_bilinear_matrix_transp_rows_loc`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_two_e_nstates_openmp_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_all_spin_singles_and_doubles_n_int`
       * :c:func:`get_all_spin_singles_n_int`
       * :c:func:`get_s2`
       * :c:func:`i_h_j_double_alpha_beta`
       * :c:func:`i_h_j_double_spin`
       * :c:func:`i_wee_j_single`

 
.. c:function:: print_energy_components:


    File : :file:`davidson/print_e_components.irp.f`

    .. code:: fortran

        subroutine print_energy_components()


    Prints the different components of the energy.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`mo_pseudo_integrals`
       * :c:data:`n_states`
       * :c:data:`nuclear_repulsion`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_energy`

 
.. c:function:: u_0_h_u_0:


    File : :file:`davidson/u0_h_u0.irp.f`

    .. code:: fortran

        subroutine u_0_H_u_0(e_0,s_0,u_0,n,keys_tmp,Nint,N_st,sze)


    Computes $E_0 = \frac{\langle u_0 | H | u_0 \rangle}{\langle u_0 | u_0 \rangle}$
    
    and      $S_0 = \frac{\langle u_0 | S^2 | u_0 \rangle}{\langle u_0 | u_0 \rangle}$
    
    n : number of determinants
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`distributed_davidson`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`n_det_max_full`
       * :c:data:`n_states_diag`
       * :c:data:`s2_matrix_all_dets`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_energy`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_nstates_openmp`
       * :c:func:`h_s2_u_0_nstates_zmq`

 
.. c:function:: u_0_h_u_0_two_e:


    File : :file:`davidson/u0_wee_u0.irp.f`

    .. code:: fortran

        subroutine u_0_H_u_0_two_e(e_0,u_0,n,keys_tmp,Nint,N_st,sze)


    Computes $E_0 = \frac{ \langle u_0 | H | u_0\rangle}{\langle u_0 | u_0 \rangle}$.
    
    n : number of determinants
    

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_energy_two_e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`h_s2_u_0_two_e_nstates_openmp`

 
.. c:function:: zmq_get_n_states_diag:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        integer function zmq_get_N_states_diag(zmq_to_qp_run_socket, worker_id)


    Get N_states_diag from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_states_diag`
       * :c:data:`zmq_state`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`n_states_diag`

 
.. c:function:: zmq_put_n_states_diag:


    File : :file:`davidson/davidson_parallel.irp.f`

    .. code:: fortran

        integer function zmq_put_N_states_diag(zmq_to_qp_run_socket,worker_id)


    Put N_states_diag on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states_diag`
       * :c:data:`zmq_state`

