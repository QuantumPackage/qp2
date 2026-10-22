.. _module_ccsd: 
 
.. program:: ccsd 
 
.. default-role:: option 
 
ccsd
====

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
EZFIO parameters 
---------------- 
 
.. option:: energy
 
    CCSD energy
 
 
.. option:: energy_t
 
    CCSD(T) energy
 
 
 
Programs 
-------- 
 
 * :ref:`ccsd` 
 * :ref:`ccsd_space_orb` 
 * :ref:`ccsd_spin_orb` 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: binary_search:


    File : :file:`ccsd_t_space_orb_stoch.irp.f`

    .. code:: fortran

        integer*8 function binary_search(arr, key, sze)


    Searches the key in array arr(1:sze) between l_in and r_in, and returns its index

 
.. c:macro:: ccsd_energy_space_x:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine ccsd_energy_space_x(nO,nV,d_cc_space_v_oovv,d_cc_space_f_vo,tau_x,t1,energy)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`gpu_default_stream`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_ddot`
       * :c:func:`gpu_ddot_64`
       * :c:func:`gpu_set_stream`
       * :c:func:`gpu_stream_create`
       * :c:func:`gpu_stream_destroy`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: ccsd_energy_spin:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine ccsd_energy_spin(nO,nV,t1,t2,Fov,v_oovv,energy)


    CCSD energy in spin orbitals

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

 
.. c:macro:: ccsd_par_t_space:


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_space(nO,nV,t1,t2,energy)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_v`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`form_v`
       * :c:func:`form_w`

 
.. c:macro:: ccsd_par_t_space_stoch:


    File : :file:`ccsd_t_space_orb_stoch.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_space_stoch(nO,nV,t1,t2,f_o,f_v,v_vvvo,v_vvoo,v_vooo,energy,t_error)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_par_t_stop`
       * :c:data:`nthreads_pt2`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsort_big`
       * :c:func:`i8set_order_big`
       * :c:func:`random_number`
       * :c:func:`wall_time`

 
.. c:macro:: ccsd_par_t_space_v2:


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_space_v2(nO,nV,t1,t2,f_o,f_v,v_vvvo,v_vvoo,v_vooo,energy)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`form_v_ijk`
       * :c:func:`form_w_ijk`
       * :c:func:`wall_time`

 
.. c:macro:: ccsd_par_t_space_v3:


    File : :file:`ccsd_t_space_orb_abc.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_space_v3(nO,nV,t1,t2,f_o,f_v,v_vvvo,v_vvoo,v_vooo,energy)



 
.. c:macro:: ccsd_par_t_spin:


    File : :file:`ccsd_t_spin_orb.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_spin(nO,nV,t1,t2,f_o,f_v,f_ov,v_ooov,v_vvoo,v_vvvo,energy)



 
.. c:macro:: ccsd_par_t_spin_v2:


    File : :file:`ccsd_t_spin_orb.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_spin_v2(nO,nV,t1,t2,f_o,f_v,f_ov,v_ooov,v_vvoo,energy)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gen_v_spin`
       * :c:func:`wall_time`

 
.. c:macro:: ccsd_t_task_aba:


    File : :file:`ccsd_t_space_orb_abc.irp.f`

    .. code:: fortran

        double precision function ccsd_t_task_aba(a,b,nO,nV,t1,T_oovv,T_voov,&
      X_ooov,X_oovv,X_vovv,f_o,f_v) result(e)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`form_v_abc`
       * :c:func:`form_w_abc`

 
.. c:macro:: ccsd_t_task_abc:


    File : :file:`ccsd_t_space_orb_abc.irp.f`

    .. code:: fortran

        double precision function ccsd_t_task_abc(a,b,c,nO,nV,t1,T_oovv,T_voov,&
      X_ooov,X_oovv,X_vovv,f_o,f_v) result(e)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`form_v_abc`
       * :c:func:`form_w_abc`

 
.. c:macro:: compute_a1:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_A1(nO,nV,t1,t2,tau,A1)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_vvoo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_a1_chol:


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_A1_chol(nO,nV,t1,t2,tau,d_cc_space_v_vooo, &
  d_cc_space_v_oooo, d_cc_space_v_vvoo, A1)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`gpu_default_stream`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space_chol`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_allocate`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_set_stream`
       * :c:func:`gpu_stream_create`
       * :c:func:`gpu_stream_destroy`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: compute_b1:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_B1(nO,nV,t1,t2,B1)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_b1_gam:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_B1_gam(nO,nV,t1,t2,B1,gam)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`gen_v_space`

 
.. c:macro:: compute_cf_oo:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cF_oo(nO,nV,t1,tau_t,Foo,Fov,v_ooov,v_oovv,cF_oo)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`
       * :c:func:`use_cf_oo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_cf_ov:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cF_ov(nO,nV,t1,Fov,v_oovv,cF_ov)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`
       * :c:func:`use_cf_ov`

 
.. c:macro:: compute_cf_vv:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cF_vv(nO,nV,t1,tau_t,Fov,Fvv,v_oovv,cF_vv)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dgemv`
       * :c:func:`gen_v_spin_3idx_ij_l`

 
.. c:macro:: compute_cw_oooo:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cW_oooo(nO,nV,t1,t2,tau,v_oooo,v_ooov,v_oovv,cW_oooo)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_cw_ovvo:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cW_ovvo(nO,nV,t1,t2,tau,v_ovvo,v_oovo,v_oovv,cW_ovvo)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`gen_v_spin_3idx_ij_l`

 
.. c:macro:: compute_cw_vvvf:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cW_vvvf(nO,nV,t1,t2,tau,f,v_vvvf,v_vovf,v_oovv,cW_vvvf)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`use_cw_vvvf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_cw_vvvv:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cW_vvvv(nO,nV,t1,t2,tau,v_vvvv,v_vovv,v_oovv,cW_vvvv)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_g_occ:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_g_occ(nO,nV,t1,t2,H_oo,g_occ)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_v_ovoo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_g_occ_chol:


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_g_occ_chol(nO,nV,t1,t2,H_oo, &
  d_cc_space_f_vo, d_cc_space_v_ov_chol, d_cc_space_v_oo_chol, d_cc_space_v_ovoo, g_occ)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space_chol`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_allocate`
       * :c:func:`gpu_copy`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_dgemv`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: compute_g_vir:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_g_vir(nO,nV,t1,t2,H_vv,g_vir)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_v_vvvo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_g_vir_chol:


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_g_vir_chol(nO,nV,t1,t2,H_vv,d_cc_space_f_vo, &
  d_cc_space_v_ov_chol, d_cc_space_v_vv_chol, g_vir)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`
       * :c:data:`gpu_default_stream`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space_chol`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_allocate`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_set_stream`
       * :c:func:`gpu_stream_create`
       * :c:func:`gpu_stream_destroy`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: compute_h_oo:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_H_oo(nO,nV,t1,t2,tau,H_oo)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_w_oovv`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_h_oo_chol:


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_H_oo_chol(nO,nV,tau_x,d_cc_space_f_oo, &
    d_cc_space_v_ov_chol,d_cc_space_v_vo_chol,H_oo)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_allocate`
       * :c:func:`gpu_blas_create`
       * :c:func:`gpu_blas_destroy`
       * :c:func:`gpu_copy`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: compute_h_vo:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_H_vo(nO,nV,t1,t2,H_vo)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_w_vvoo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemv`

 
.. c:macro:: compute_h_vo_chol:


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_H_vo_chol(nO,nV,t1,d_cc_space_f_vo, &
         d_cc_space_v_ov_chol,d_cc_space_v_vo_chol, H_vo)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`
       * :c:data:`gpu_default_stream`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_allocate`
       * :c:func:`gpu_copy`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_set_stream`
       * :c:func:`gpu_stream_create`
       * :c:func:`gpu_stream_destroy`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: compute_h_vv:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_H_vv(nO,nV,t1,t2,tau,H_vv)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_w_vvoo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_h_vv_chol:


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_H_vv_chol(nO,nV,tau_x,d_cc_space_f_vv, &
         d_cc_space_v_ov_chol,H_vv)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_allocate`
       * :c:func:`gpu_blas_create`
       * :c:func:`gpu_blas_destroy`
       * :c:func:`gpu_copy`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: compute_j1:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_J1(nO,nV,t1,t2,v_ovvo,v_ovoo,v_vvvo,v_vvoo,J1)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_j1_chol:


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_J1_chol(nO,nV,t1,t2,v_ovvo,v_ovoo,v_vvoo,d_cc_space_v_vo_chol,d_cc_space_v_vv_chol,J1)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`
       * :c:data:`gpu_default_stream`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space_chol`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_allocate`
       * :c:func:`gpu_copy`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_set_stream`
       * :c:func:`gpu_stream_create`
       * :c:func:`gpu_stream_destroy`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: compute_k1:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_K1(nO,nV,t1,t2,v_ovoo,v_vvoo,v_ovov,v_vvov,K1)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: compute_k1_chol:


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_K1_chol(nO,nV,t1,t2,v_ovoo,v_vvoo,v_ovov, &
   d_cc_space_v_ov_chol,d_cc_space_v_vv_chol,K1)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`
       * :c:data:`gpu_default_stream`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space_chol`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_allocate`
       * :c:func:`gpu_copy`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_set_stream`
       * :c:func:`gpu_stream_create`
       * :c:func:`gpu_stream_destroy`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: compute_r1_space:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_r1_space(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r1,max_r1)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vvov`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dgemv`

 
.. c:macro:: compute_r1_space_chol:


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_r1_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r1,max_r1,d_cc_space_f_ov,d_cc_space_f_vo, &
    d_cc_space_v_voov, d_cc_space_v_ovov, d_cc_space_v_oovo, d_cc_space_v_vo_chol, d_cc_space_v_vv_chol)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`cholesky_mo_num`
       * :c:data:`gpu_default_stream`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_allocate`
       * :c:func:`gpu_copy`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_dgemv`
       * :c:func:`gpu_set_stream`
       * :c:func:`gpu_stream_create`
       * :c:func:`gpu_stream_destroy`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: compute_r1_spin:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_r1_spin(nO,nV,t1,t2,f_o,f_v,Fov,cF_oo,cF_ov,cF_vv,v_oovo,v_ovov,r1)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`gen_v_spin_3idx`

 
.. c:macro:: compute_r2_space:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_r2_space(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r2,max_r2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_a1`
       * :c:func:`compute_b1_gam`
       * :c:func:`compute_g_occ`
       * :c:func:`compute_g_vir`
       * :c:func:`compute_j1`
       * :c:func:`compute_k1`
       * :c:func:`dgemm`

 
.. c:macro:: compute_r2_space_chol:


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_r2_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv, &
    d_cc_space_v_oovv, d_cc_space_v_vooo, d_cc_space_v_oooo, d_cc_space_v_oovo, d_cc_space_v_ovvo, d_cc_space_v_ovoo, &
    d_cc_space_v_ovov, d_cc_space_v_vvoo, d_cc_space_v_oo_chol, d_cc_space_v_ov_chol, d_cc_space_v_vo_chol, d_cc_space_v_vv_chol, &
    d_cc_space_f_vo, &
    r2,max_r2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cholesky_mo_num`
       * :c:data:`gpu_default_stream`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_a1_chol`
       * :c:func:`compute_g_occ_chol`
       * :c:func:`compute_g_vir_chol`
       * :c:func:`compute_j1_chol`
       * :c:func:`compute_k1_chol`
       * :c:func:`gpu_allocate`
       * :c:func:`gpu_blas_create`
       * :c:func:`gpu_blas_destroy`
       * :c:func:`gpu_copy`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_dgemm`
       * :c:func:`gpu_set_stream`
       * :c:func:`gpu_stream_create`
       * :c:func:`gpu_stream_destroy`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: compute_r2_spin:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_r2_spin(nO,nV,t1,t2,tau,f_o,f_v,cF_oo,cF_ov,cF_vv,cW_oooo,cW_ovvo,v_ovoo,v_oovv,v_ovvo,r2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_dev`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`gen_v_spin_3idx_i_kl`
       * :c:func:`use_cw_vvvf`
       * :c:func:`wall_time`

 
.. c:macro:: compute_tau_spin:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_tau_spin(nO,nV,t1,t2,tau)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

 
.. c:macro:: compute_tau_t_spin:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_tau_t_spin(nO,nV,t1,t2,tau_t)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

 
.. c:macro:: form_v:


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine form_v(nO,nV,t1,w,v)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_vvoo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_par_t_space`

 
.. c:macro:: form_v_abc:


    File : :file:`ccsd_t_space_orb_abc.irp.f`

    .. code:: fortran

        subroutine form_v_abc(nO,nV,a,b,c,T_ov,X_oovv,W_abc,V_abc,W_cba,V_cba,W_bca,V_bca,W_cab,V_cab,W_bac,V_bac,W_acb,V_acb)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_t_task_aba`
       * :c:func:`ccsd_t_task_abc`

 
.. c:macro:: form_v_ijk:


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine form_v_ijk(nO,nV,i,j,k,T_vo,X_vvoo,w,v)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_par_t_space_v2`

 
.. c:macro:: form_w:


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine form_w(nO,nV,t2,W)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_vvvo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_par_t_space`

 
.. c:macro:: form_w_abc:


    File : :file:`ccsd_t_space_orb_abc.irp.f`

    .. code:: fortran

        subroutine form_w_abc(nO,nV,a,b,c,T_voov,T_oovv,X_vovv,X_ooov,W_abc,W_cba,W_bca,W_cab,W_bac,W_acb)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_t_task_aba`
       * :c:func:`ccsd_t_task_abc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: form_w_ijk:


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine form_w_ijk(nO,nV,i,j,k,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_par_t_space_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:macro:: run_ccsd_space_orb:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine run_ccsd_space_orb



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_diis_depth`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_max_iter`
       * :c:data:`cc_no_m`
       * :c:data:`cc_par_t`
       * :c:data:`cc_ref`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_thresh_conv`
       * :c:data:`cc_update_method`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`qp_max_mem`
       * :c:data:`read_amplitudes`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd`
       * :c:func:`ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`abort`
       * :c:func:`ccsd_energy_space_x`
       * :c:func:`ccsd_par_t_space_stoch`
       * :c:func:`compute_h_oo`
       * :c:func:`compute_h_oo_chol`
       * :c:func:`compute_h_vo`
       * :c:func:`compute_h_vo_chol`
       * :c:func:`compute_h_vv`
       * :c:func:`compute_h_vv_chol`
       * :c:func:`compute_r1_space`
       * :c:func:`compute_r1_space_chol`
       * :c:func:`compute_r2_space`
       * :c:func:`compute_r2_space_chol`
       * :c:func:`det_energy`
       * :c:func:`ezfio_set_utils_cc_io_amplitudes`
       * :c:func:`gpu_allocate`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_upload`
       * :c:func:`guess_t1`
       * :c:func:`guess_t2`
       * :c:func:`print_det`
       * :c:func:`resident_memory`
       * :c:func:`save_energy`
       * :c:func:`update_t1`
       * :c:func:`update_t2`
       * :c:func:`update_t_ccsd_diis_v3`
       * :c:func:`update_tau_space`
       * :c:func:`update_tau_x_space`
       * :c:func:`wall_time`
       * :c:func:`write_t1`
       * :c:func:`write_t2`

 
.. c:macro:: run_ccsd_spin_orb:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine run_ccsd_spin_orb


    CCSD in spin orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_dev`
       * :c:data:`cc_diis_depth`
       * :c:data:`cc_max_iter`
       * :c:data:`cc_no_m`
       * :c:data:`cc_par_t`
       * :c:data:`cc_ref`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`
       * :c:data:`cc_thresh_conv`
       * :c:data:`cc_update_method`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`
       * :c:data:`read_amplitudes`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd`
       * :c:func:`ccsd_spin_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_energy_spin`
       * :c:func:`ccsd_par_t_spin_v2`
       * :c:func:`compute_cf_oo`
       * :c:func:`compute_cf_ov`
       * :c:func:`compute_cf_vv`
       * :c:func:`compute_cw_oooo`
       * :c:func:`compute_cw_ovvo`
       * :c:func:`compute_r1_spin`
       * :c:func:`compute_r2_spin`
       * :c:func:`compute_tau_spin`
       * :c:func:`compute_tau_t_spin`
       * :c:func:`det_energy`
       * :c:func:`ezfio_set_utils_cc_io_amplitudes`
       * :c:func:`gpu_allocate`
       * :c:func:`gpu_deallocate`
       * :c:func:`gpu_upload`
       * :c:func:`guess_t1`
       * :c:func:`guess_t2`
       * :c:func:`print_det`
       * :c:func:`save_energy`
       * :c:func:`update_t1`
       * :c:func:`update_t2`
       * :c:func:`update_t_ccsd_diis_v3`
       * :c:func:`wall_time`
       * :c:func:`write_t1`
       * :c:func:`write_t2`

 
.. c:macro:: save_energy:


    File : :file:`save_energy.irp.f`

    .. code:: fortran

        subroutine save_energy(E,ET)


    Saves the energy in |EZFIO|.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`
       * :c:func:`run_ccsd_spin_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_ccsd_energy`
       * :c:func:`ezfio_set_ccsd_energy_t`

 
.. c:macro:: update_tau_space:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine update_tau_space(nO,nV,h_t1,t1,t2,tau)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`gpu_default_stream`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_set_stream`
       * :c:func:`gpu_stream_create`
       * :c:func:`gpu_stream_destroy`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: update_tau_x_space:


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine update_tau_x_space(nO,nV,tau,tau_x)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`blas_handle`
       * :c:data:`gpu_default_stream`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_space_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`gpu_dgeam`
       * :c:func:`gpu_set_stream`
       * :c:func:`gpu_stream_create`
       * :c:func:`gpu_stream_destroy`
       * :c:func:`gpu_synchronize`

 
.. c:macro:: use_cf_oo:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine use_cF_oo(nO,nV,t1,t2,tau_t,F_oo,F_ov,v_ooov,v_oovv,r1,r2)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_cf_oo`
       * :c:func:`dgemm`

 
.. c:macro:: use_cf_ov:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine use_cF_ov(nO,nV,t1,t2,F_ov,v_oovv,r1,r2)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_cf_ov`
       * :c:func:`dgemm`

 
.. c:macro:: use_cf_vv:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine use_cF_vv(nO,nV,t1,t2,r1,r2)



 
.. c:macro:: use_cw_vvvf:


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine use_cW_vvvf(nO,nV,t1,t2,tau,v_oovv,r2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_dev`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_spin`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_cw_vvvf`
       * :c:func:`dgemm`
       * :c:func:`gen_v_spin_3idx`
       * :c:func:`wall_time`

