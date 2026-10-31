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
 
 * :ref:`ccsd_space_orb` 
 * :ref:`ccsd_spin_orb` 
 * :ref:`ccsd` 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void binary_search(void)


    File : :file:`ccsd_t_space_orb_stoch.irp.f`

    .. code:: fortran

        integer*8 function binary_search(arr, key, sze)


    Searches the key in array arr(1:sze) between l_in and r_in, and returns its index

 
.. c:function:: void ccsd_energy_space_x(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine ccsd_energy_space_x(nO,nV,d_cc_space_v_oovv,d_cc_space_f_vo,tau_x,t1,energy)


    CCSD energy: energy = tau_x : v_oovv + 2 * f_vo : t1, with
    f_vo and v_oovv the vo Fock block and oovv two-electron
    integrals, nO, nV the occupied and virtual spatial orbital
    counts. It is added to the reference determinant energy.
    tau_x may be contracted with v_oovv in place of tau because
    v_oovv(i,j,a,b) = <ij|ab> is symmetric in (a, b) and
    update_tau_x_space computes tau_x(i,j,a,b) = 2 * tau(i,j,a,b)
    - tau(i,j,b,a), so the contraction v_oovv : tau_x is exactly
    equal to v_oovv : tau.

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

 
.. c:function:: void ccsd_energy_spin(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine ccsd_energy_spin(nO,nV,t1,t2,Fov,v_oovv,energy)


    CCSD energy in spin orbitals

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

 
.. c:function:: void ccsd_par_t_space(void)


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_space(nO,nV,t1,t2,energy)


    (T) correction to the CCSD energy (CCSD(T)), in spatial orbitals
    (closed shell systems). form_w and form_v build the W and V
    triples from T1 and T2, integrals and f_o, f_v from the
    cc_space_* module providers. energy is the sum over i,j,k,a,b,c
    of (4W(a,b,c)+W(b,c,a)+W(c,a,b))*(V(a,b,c)-V(c,b,a)) /
    (f_o(i)+f_o(j)+f_o(k)-f_v(a)-f_v(b)-f_v(c)), divided by 3.

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

 
.. c:function:: void ccsd_par_t_space_stoch(void)


    File : :file:`ccsd_t_space_orb_stoch.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_space_stoch(nO,nV,t1,t2,f_o,f_v,v_vvvo,v_vvoo,v_vooo,energy,t_error)


    Stochastic (T) correction to the CCSD energy (CCSD(T)), in
    spatial orbitals (closed shell systems). Virtual triplets are
    sampled with probability proportional to 1/max(0.2,
    f_v(a)+f_v(b)+f_v(c)). On entry energy is the CCSD total energy,
    on exit the (T) correction; t_error is the statistical error bar
    and the loop stops when it is below cc_par_t_stop.

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

 
.. c:function:: void ccsd_par_t_space_v2(void)


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_space_v2(nO,nV,t1,t2,f_o,f_v,v_vvvo,v_vvoo,v_vooo,energy)


    (T) correction to the CCSD energy (CCSD(T)), in spatial orbitals
    (closed shell systems), OpenMP version of ccsd_par_t_space. The
    integrals and f_o, f_v are passed as arguments and permuted into
    work arrays; the W and V slices for one (i,j,k) are built with
    form_w_ijk and form_v_ijk. energy is the same (T) sum, divided
    by 3.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`form_v_ijk`
       * :c:func:`form_w_ijk`
       * :c:func:`wall_time`

 
.. c:function:: void ccsd_par_t_space_v3(void)


    File : :file:`ccsd_t_space_orb_abc.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_space_v3(nO,nV,t1,t2,f_o,f_v,v_vvvo,v_vvoo,v_vooo,energy)


    (T) correction to the CCSD energy (CCSD(T)), in spatial orbitals
    (closed shell systems). Loops over the virtual triples: for each
    a<b<c it calls ccsd_t_task_abc, and for each a<b it calls
    ccsd_t_task_aba(a,b) and ccsd_t_task_aba(b,a) (two equal virtual
    orbitals), OpenMP parallel over the triples. energy is the sum,
    divided by 3.

 
.. c:function:: void ccsd_par_t_spin(void)


    File : :file:`ccsd_t_spin_orb.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_spin(nO,nV,t1,t2,f_o,f_v,f_ov,v_ooov,v_vvoo,v_vvvo,energy)


    (T) correction to the CCSD energy (CCSD(T)), in spin orbitals (v1).
    Builds the triple amplitudes t3 from T2 against the vvvo and
    ooov integrals, divided by f_o(i)+f_o(j)+f_o(k)-f_v(a)-f_v(b)-
    f_v(c), then adds to energy the terms e_t (t3*delta*t3, /36),
    e_st (T1 against v_vvoo acting on t3, *1/4) and e_dt (T2 against
    f_ov times t3, *1/4).

 
.. c:function:: void ccsd_par_t_spin_v2(void)


    File : :file:`ccsd_t_spin_orb.irp.f`

    .. code:: fortran

        subroutine ccsd_par_t_spin_v2(nO,nV,t1,t2,f_o,f_v,f_ov,v_ooov,v_vvoo,energy)


    (T) correction to the CCSD energy (CCSD(T)), in spin orbitals
    (v2), OpenMP version of ccsd_par_t_spin: the triple amplitudes
    t3 are built for one (b,c) virtual pair at a time, and the vvvo
    integrals are generated in the module with gen_v_spin instead of
    being passed as an argument.

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

 
.. c:function:: void ccsd_t_task_aba(void)


    File : :file:`ccsd_t_space_orb_abc.irp.f`

    .. code:: fortran

        double precision function ccsd_t_task_aba(a,b,nO,nV,t1,T_oovv,T_voov,&
              X_ooov,X_oovv,X_vovv,f_o,f_v) result(e)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`form_v_abc`
       * :c:func:`form_w_abc`

 
.. c:function:: void ccsd_t_task_abc(void)


    File : :file:`ccsd_t_space_orb_abc.irp.f`

    .. code:: fortran

        double precision function ccsd_t_task_abc(a,b,c,nO,nV,t1,T_oovv,T_voov,&
              X_ooov,X_oovv,X_vovv,f_o,f_v) result(e)



    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`form_v_abc`
       * :c:func:`form_w_abc`

 
.. c:function:: void compute_a1(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_A1(nO,nV,t1,t2,tau,A1)


    Computes A1(u,v,i,j), an intermediate of the T2 residual that
    is contracted with tau over the occupied indices (i, j):
    A1(u,v,i,j) = v_oooo(u,v,i,j) + sum over a of
    v_ovoo(u,a,i,j) * t1(v,a) + sum over a of v_vooo(a,v,i,j)
    * t1(u,a) + sum over (a, b) of v_vvoo(a,b,i,j) * tau(u,v,a,b).

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

 
.. c:function:: void compute_a1_chol(void)


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_A1_chol(nO,nV,t1,t2,tau,d_cc_space_v_vooo, &
          d_cc_space_v_oooo, d_cc_space_v_vvoo, A1)


    Computes A1(u,v,i,j), an intermediate of the T2 residual,
    the GPU/Cholesky counterpart of compute_A1, called from
    compute_r2_space_chol. A1 is built from d_cc_space_v_oooo,
    d_cc_space_v_vooo and d_cc_space_v_vvoo contracted with t1
    and tau, as sketched in the in-code comments; the v_ovoo
    term of compute_A1 is formed from d_cc_space_v_vooo with a
    transposed index.
    The BLAS formulation reproduces the four-term formula of
    compute_A1 (author-confirmed 2026-10-31).

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

 
.. c:function:: void compute_b1(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_B1(nO,nV,t1,t2,B1)


    Computes the full B1 intermediate of the T2 residual:
    B1(a,b,beta,gam) = v_vvvv(a,b,beta,gam) - sum over i of
    v_vvvo(a,b,beta,i) * t1(i,gam)
    - sum over i of v_vvov(a,b,i,gam) * t1(i,beta).
    Legacy alternative to compute_B1_gam: its only call in this
    file is commented out, and the subroutine is kept as
    documentation (author, 2026-10-31).

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

 
.. c:function:: void compute_b1_gam(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_B1_gam(nO,nV,t1,t2,B1,gam)


    Computes the slice at a fixed virtual index gam of the B1
    intermediate of the T2 residual:
    B1(a,b,beta) = v_vvvv(a,b,beta,gam) - sum over i of
    v_vvvo(a,b,beta,i) * t1(i,gam)
    - sum over i of v_vvov(a,b,i,gam) * t1(i,beta).

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

 
.. c:function:: void compute_cf_oo(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cF_oo(nO,nV,t1,tau_t,Foo,Fov,v_ooov,v_oovv,cF_oo)


    CCSD cF_oo intermediate in spin orbitals (i,m,n occupied,
    e,f virtual spin orbitals):
    cF_oo(m,i) = Foo(m,i) with zero diagonal, plus
    sum_e t1(n,e)*v_ooov(m,n,i,e) + 0.5*sum_e Fov(m,e)*t1(i,e)
    + 0.5*sum v_oovv(m,n,e,f)*tau_t(i,n,e,f).

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`
       * :c:func:`use_cf_oo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void compute_cf_ov(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cF_ov(nO,nV,t1,Fov,v_oovv,cF_ov)


    CCSD cF_ov intermediate in spin orbitals (m occupied,
    e virtual spin orbital):
    cF_ov(m,e) = Fov(m,e) + sum t1(n,f)*v_oovv(m,n,e,f).

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`
       * :c:func:`use_cf_ov`

 
.. c:function:: void compute_cf_vv(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cF_vv(nO,nV,t1,tau_t,Fov,Fvv,v_oovv,cF_vv)


    CCSD cF_vv intermediate in spin orbitals (a,e virtual spin
    orbitals): cF_vv(a,e) = Fvv(a,e) with zero diagonal, minus
    0.5*sum_m t1(m,a)*Fov(m,e), plus sum_m t1(m,f)*v_ovvv(m,a,f,e),
    minus 0.5*sum tau_t(m,n,a,f)*v_oovv(m,n,e,f).

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

 
.. c:function:: void compute_cw_oooo(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cW_oooo(nO,nV,t1,t2,tau,v_oooo,v_ooov,v_oovv,cW_oooo)


    CCSD cW_oooo intermediate in spin orbitals (i,j,m,n
    occupied):
    cW_oooo(m,n,i,j) = v_oooo(m,n,i,j)
      + sum_e t1(j,e)*v_ooov(m,n,i,e) - t1(i,e)*v_ooov(m,n,j,e)
      + 0.25*sum tau(i,j,e,f)*v_oovv(m,n,e,f).

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void compute_cw_ovvo(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cW_ovvo(nO,nV,t1,t2,tau,v_ovvo,v_oovo,v_oovv,cW_ovvo)


    CCSD cW_ovvo intermediate in spin orbitals (i,j occupied,
    a,b virtual):
    cW_ovvo(i,a,b,j) = v_ovvo(i,a,b,j)
      + sum_f t1(j,f)*v_ovvv(i,a,b,f)
      - sum_n t1(n,b)*v_oovo(i,n,a,j)
      - sum (0.5*t2(j,n,f,b) + t1(j,f)*t1(n,b))*v_oovv(i,n,a,f).

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

 
.. c:function:: void compute_cw_vvvf(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cW_vvvf(nO,nV,t1,t2,tau,f,v_vvvf,v_vovf,v_oovv,cW_vvvf)


    CCSD cW_vvvv intermediate in spin orbitals, restricted to a
    fixed last index f (a,b,e virtual spin orbitals):
    cW_vvvf(a,b,e) = v_vvvf(a,b,e)
      - sum_m t1(m,b)*v_vovf(a,m,e) + t1(m,a)*v_vovf(b,m,e)
      + 0.25*sum tau(m,n,a,b)*v_oovv(m,n,e,f).

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`use_cw_vvvf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void compute_cw_vvvv(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_cW_vvvv(nO,nV,t1,t2,tau,v_vvvv,v_vovv,v_oovv,cW_vvvv)


    CCSD cW_vvvv intermediate in spin orbitals (a,b,c,d virtual
    spin orbitals). Not called by run_ccsd_spin_orb, whose
    cW_vvvv contribution to r2 is applied through use_cW_vvvf.
    cW_vvvv(a,b,c,d) = v_vvvv(a,b,c,d)
      - sum_m t1(m,b)*v_vovv(a,m,c,d) + t1(m,a)*v_vovv(b,m,c,d)
      + 0.25*sum tau(m,n,a,b)*v_oovv(m,n,c,d).

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void compute_g_occ(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_g_occ(nO,nV,t1,t2,H_oo,g_occ)


    Computes g_occ(u,i), an occupied-space intermediate of the T2
    residual:
    g_occ(u,i) = H_oo(u,i) + sum over a of f_vo(a,i) * t1(u,a)
    + sum over (j, a) of (2 * v_ovoo(u,a,i,j) - v_ovoo(u,a,j,i))
    * t1(j,a).

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

 
.. c:function:: void compute_g_occ_chol(void)


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_g_occ_chol(nO,nV,t1,t2,H_oo, &
          d_cc_space_f_vo, d_cc_space_v_ov_chol, d_cc_space_v_oo_chol, d_cc_space_v_ovoo, g_occ)


    Computes g_occ(u,i), an occupied-space intermediate of the
    T2 residual, the GPU/Cholesky counterpart of compute_g_occ,
    called from compute_r2_space_chol. g_occ is initialized to
    H_oo and updated with the t1 and d_cc_space_f_vo product, a
    t1 contraction over the Cholesky vectors d_cc_space_v_ov_chol
    and d_cc_space_v_oo_chol, and a d_cc_space_v_ovoo contraction
    with t1.
    The BLAS sequence reproduces the compute_g_occ formula
    H_oo(u,i) + sum f_vo(a,i) * t1(u,a) + sum (2 *
    v_ovoo(u,a,i,j) - v_ovoo(u,a,j,i)) * t1(j,a)
    (author-confirmed 2026-10-31).

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

 
.. c:function:: void compute_g_vir(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_g_vir(nO,nV,t1,t2,H_vv,g_vir)


    Computes g_vir(a,beta), a virtual-space intermediate of the
    T2 residual:
    g_vir(a,beta) = H_vv(a,beta) - sum over i of f_vo(a,i)
    * t1(i,beta) + sum over (i, b) of
    (2 * v_vvvo(a,b,beta,i) - v_vvvo(b,a,beta,i)) * t1(i,b).

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

 
.. c:function:: void compute_g_vir_chol(void)


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_g_vir_chol(nO,nV,t1,t2,H_vv,d_cc_space_f_vo, &
          d_cc_space_v_ov_chol, d_cc_space_v_vv_chol, g_vir)


    Computes g_vir(a,beta), a virtual-space intermediate of the
    T2 residual, the GPU/Cholesky counterpart of compute_g_vir,
    called from compute_r2_space_chol. g_vir is built from
    -d_cc_space_f_vo : t1, t1 contractions over the Cholesky
    vectors d_cc_space_v_ov_chol and d_cc_space_v_vv_chol, and
    H_vv.
    The BLAS sequence reproduces the compute_g_vir formula
    (H_vv - f_vo : t1 + the v_vvvo terms); the pre-GPU signed
    transpose of tmp_vo (commit fa6d1419~1) is one gpu_dgeam
    (author-confirmed 2026-10-31).

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

 
.. c:function:: void compute_h_oo(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_H_oo(nO,nV,t1,t2,tau,H_oo)


    Computes the block H_oo used in the T1 residual:
    H_oo(u,i) = f_oo(u,i) + sum over (j, a, b) of
    w_oovv(i,j,a,b) * tau(u,j,a,b).

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

 
.. c:function:: void compute_h_oo_chol(void)


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_H_oo_chol(nO,nV,tau_x,d_cc_space_f_oo, &
            d_cc_space_v_ov_chol,d_cc_space_v_vo_chol,H_oo)


    Computes the block H_oo used in the T1 residual, the
    GPU/Cholesky counterpart of compute_H_oo, called from
    run_ccsd_space_orb when do_mo_cholesky is set. H_oo is
    initialized to d_cc_space_f_oo and updated with tau_x
    contractions over d_cc_space_v_ov_chol and d_cc_space_v_vo_chol.
    The BLAS sequence reproduces the compute_H_oo formula
    f_oo(u,i) + sum w_oovv(i,j,a,b) * tau(u,j,a,b), the
    two-electron contraction formed through the Cholesky vectors
    (author-confirmed 2026-10-31).

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

 
.. c:function:: void compute_h_vo(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_H_vo(nO,nV,t1,t2,H_vo)


    Computes the block H_vo used in the T1 residual:
    H_vo(a,i) = f_vo(a,i) + sum over (j, b) of
    w_vvoo(a,b,i,j) * t1(j,b).

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

 
.. c:function:: void compute_h_vo_chol(void)


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_H_vo_chol(nO,nV,t1,d_cc_space_f_vo, &
                 d_cc_space_v_ov_chol,d_cc_space_v_vo_chol, H_vo)


    Computes the block H_vo used in the T1 residual, the
    GPU/Cholesky counterpart of compute_H_vo, called from
    run_ccsd_space_orb when do_mo_cholesky is set. H_vo is
    initialized to d_cc_space_f_vo and updated with t1
    contractions over d_cc_space_v_ov_chol and d_cc_space_v_vo_chol.
    The BLAS sequence reproduces the compute_H_vo formula
    f_vo(a,i) + sum w_vvoo(a,b,i,j) * t1(j,b); the pre-GPU
    transpose loop (commit fa6d1419~1) is a gpu_dgeam that
    copies each column in transposed order (author-confirmed
    2026-10-31).

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

 
.. c:function:: void compute_h_vv(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_H_vv(nO,nV,t1,t2,tau,H_vv)


    Computes the block H_vv used in the T1 residual:
    H_vv(a,beta) = f_vv(a,beta) - sum over (i, j, b) of
    w_vvoo(a,b,i,j) * tau(i,j,beta,b).

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

 
.. c:function:: void compute_h_vv_chol(void)


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_H_vv_chol(nO,nV,tau_x,d_cc_space_f_vv, &
                 d_cc_space_v_ov_chol,H_vv)


    Computes the block H_vv used in the T1 residual, the
    GPU/Cholesky counterpart of compute_H_vv, called from
    run_ccsd_space_orb when do_mo_cholesky is set. H_vv is
    initialized to d_cc_space_f_vv and updated with tau_x
    contractions over d_cc_space_v_ov_chol.
    The BLAS sequence reproduces the compute_H_vv formula
    f_vv(a,beta) - sum w_vvoo(a,b,i,j) * tau(i,j,beta,b), the
    contraction formed through d_cc_space_v_ov_chol
    (author-confirmed 2026-10-31).

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

 
.. c:function:: void compute_j1(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_J1(nO,nV,t1,t2,v_ovvo,v_ovoo,v_vvvo,v_vvoo,J1)


    Computes J1(u,a,beta,i), an intermediate of the T2 residual,
    as a sum of five terms over v_ovvo, v_ovoo, v_vvvo and v_vvoo
    contracted with t1 and t2; the full index structure is in the
    commented reference loops at the top of the body.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void compute_j1_chol(void)


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_J1_chol(nO,nV,t1,t2,v_ovvo,v_ovoo,v_vvoo,d_cc_space_v_vo_chol,d_cc_space_v_vv_chol,J1)


    Computes J1(u,a,beta,i), an intermediate of the T2 residual,
    the GPU/Cholesky counterpart of compute_J1, called from
    compute_r2_space_chol. J1 is initialized to d_cc_space_v_ovvo
    and updated from d_cc_space_v_ovoo, d_cc_space_v_vvoo, t1,
    t2 and the Cholesky vectors d_cc_space_v_vo_chol and
    d_cc_space_v_vv_chol.
    The BLAS sequence reproduces the compute_J1 formula; the
    pre-GPU explicit expressions for Y_ovov and X_voov (commit
    fa6d1419~1) are gpu_dgeam calls on t2 and v_vvoo slices
    (author-confirmed 2026-10-31).

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

 
.. c:function:: void compute_k1(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_K1(nO,nV,t1,t2,v_ovoo,v_vvoo,v_ovov,v_vvov,K1)


    Computes K1(u,a,i,beta), an intermediate of the T2 residual,
    as a sum of four terms over v_ovov, v_ovoo, v_vvov and v_vvoo
    contracted with t1 and t2; the full index structure is in the
    commented reference loops at the top of the body.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`compute_r2_space`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void compute_k1_chol(void)


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_K1_chol(nO,nV,t1,t2,v_ovoo,v_vvoo,v_ovov, &
           d_cc_space_v_ov_chol,d_cc_space_v_vv_chol,K1)


    Computes K1(u,a,i,beta), an intermediate of the T2 residual,
    the GPU/Cholesky counterpart of compute_K1, called from
    compute_r2_space_chol. K1 is initialized to d_cc_space_v_ovov
    and updated from d_cc_space_v_ovoo, d_cc_space_v_vvoo, t1,
    t2 and the Cholesky vectors d_cc_space_v_ov_chol and
    d_cc_space_v_vv_chol.
    The BLAS sequence reproduces the compute_K1 formula; the
    pre-GPU explicit loops for Y and for the K1tmp accumulation
    (commit fa6d1419~1) are gpu_dgeam calls that a standard
    BLAS reading of the arguments does not reproduce
    (author-confirmed 2026-10-31).

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

 
.. c:function:: void compute_r1_space(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_r1_space(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r1,max_r1)


    Computes the residual r1 of the T1 equation, from the Fock
    blocks cc_space_f_ov, cc_space_f_vo, the two-electron
    integrals (global cc_space_v_*), t1, t2, tau, and the blocks
    H_oo, H_vv, H_vo.
    r1 is sign-flipped at the end for consistency with the spin-
    orbital code; max_r1 is the largest |r1| element.

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

 
.. c:function:: void compute_r1_space_chol(void)


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_r1_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r1,max_r1,d_cc_space_f_ov,d_cc_space_f_vo, &
            d_cc_space_v_voov, d_cc_space_v_ovov, d_cc_space_v_oovo, d_cc_space_v_vo_chol, d_cc_space_v_vv_chol)


    Computes the residual r1 of the T1 equation, the
    GPU/Cholesky counterpart of compute_r1_space, called from
    run_ccsd_space_orb when do_mo_cholesky is set. r1 is
    initialized to d_cc_space_f_ov and updated from t1, t2, tau,
    H_oo, H_vv, H_vo, d_cc_space_f_vo, the integrals v_voov,
    v_ovov, v_oovo, and the Cholesky vectors v_vo_chol and
    v_vv_chol. r1 is sign-flipped at the end for consistency
    with the spin-orbital code; max_r1 is the largest |r1|.
    The BLAS sequence reproduces the non-Cholesky r1 formula of
    compute_r1_space; the pre-GPU X_voov expression (commit
    fa6d1419~1) is two gpu_dgeam calls on t2 and t1 slices
    (author-confirmed 2026-10-31).

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

 
.. c:function:: void compute_r1_spin(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_r1_spin(nO,nV,t1,t2,f_o,f_v,Fov,cF_oo,cF_ov,cF_vv,v_oovo,v_ovov,r1)


    CCSD spin-orbital residual r1(i,a) for the T1 amplitudes,
    from t1, t2, the Fock matrix diagonal f_o, f_v, Fov and the
    cF_oo, cF_ov, cF_vv intermediates. The result is
    (f_o(i)-f_v(a))*t1(i,a) minus the accumulated terms,
    which start from Fov(i,a).

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

 
.. c:function:: void compute_r2_space(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_r2_space(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r2,max_r2)


    Computes the residual r2 of the T2 equation, from t1, t2, tau,
    H_oo, H_vv, H_vo, and the intermediates g_occ, g_vir, A1, B1
    (via compute_B1_gam), J1, K1.
    r2 is sign-flipped at the end for consistency with the spin-
    orbital code; max_r2 is the largest absolute value of r2 over
    its (i,j,a,b) indices, i.e. max_r2 = max |r2|, consistent with
    max_r1 (dabs) and with the spin-orbital reference
    ccsd_spin_orb_sub.irp.f, which uses maxval(abs) for both
    residuals. Standard CC practice monitors a norm of the
    absolute residual (CFOUR keyword CC_CONV, Q-Chem keyword
    CC_T_CONV), which is what max_r2 does. The missing-dabs defect
    that let max_r2 come out smaller than max |r2|, or even
    negative, was fixed in issue #396 (2026-10-31): the dabs was
    added to the accumulation.

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

 
.. c:function:: void compute_r2_space_chol(void)


    File : :file:`ccsd_space_orb_sub_chol.irp.f`

    .. code:: fortran

        subroutine compute_r2_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv, &
            d_cc_space_v_oovv, d_cc_space_v_vooo, d_cc_space_v_oooo, d_cc_space_v_oovo, d_cc_space_v_ovvo, d_cc_space_v_ovoo, &
            d_cc_space_v_ovov, d_cc_space_v_vvoo, d_cc_space_v_oo_chol, d_cc_space_v_ov_chol, d_cc_space_v_vo_chol, d_cc_space_v_vv_chol, &
            d_cc_space_f_vo, &
            r2,max_r2)


    Computes the residual r2 of the T2 equation, the
    GPU/Cholesky counterpart of compute_r2_space, called from
    run_ccsd_space_orb when do_mo_cholesky is set. r2 is
    initialized to the integral d_cc_space_v_oovv and updated
    from t1, t2, tau, the global integrals, and the Cholesky
    vectors, using compute_A1_chol, compute_g_occ_chol,
    compute_g_vir_chol, compute_J1_chol and compute_K1_chol.
    r2 is sign-flipped at the end for consistency with the
    spin-orbital code; max_r2 is the largest absolute value of r2
    over its indices, i.e. max_r2 = max |r2|, consistent with
    max_r1 (dabs) and with the spin-orbital reference. The
    dabs is taken in the per-thread leaf max_r2_local (see
    compute_r2_space), so it cannot be moved to the reduction:
    max(|a|,|b|) != |max(a,b)|. The missing-dabs defect was fixed
    in issue #396 (2026-10-31): the dabs was added to the leaf
    accumulation.

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

 
.. c:function:: void compute_r2_spin(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_r2_spin(nO,nV,t1,t2,tau,f_o,f_v,cF_oo,cF_ov,cF_vv,cW_oooo,cW_ovvo,v_ovoo,v_oovv,v_ovvo,r2)


    CCSD spin-orbital residual r2(i,j,a,b) for the T2
    amplitudes, from t1, t2, tau, the Fock matrix diagonal f_o,
    f_v and the cF_oo, cF_ov, cF_vv, cW_oooo and cW_ovvo
    intermediates. The result is
    (f_o(i)+f_o(j)-f_v(a)-f_v(b))*t2(i,j,a,b) minus the
    accumulated terms, which start from v_oovv(i,j,a,b).

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

 
.. c:function:: void compute_tau_spin(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_tau_spin(nO,nV,t1,t2,tau)


    CCSD tau intermediate in spin orbitals:
    tau(i,j,a,b) = t2(i,j,a,b) + t1(i,a)*t1(j,b)
                  - t1(i,b)*t1(j,a),
    i,j occupied and a,b virtual spin orbitals.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

 
.. c:function:: void compute_tau_t_spin(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine compute_tau_t_spin(nO,nV,t1,t2,tau_t)


    CCSD tau_t intermediate in spin orbitals:
    tau_t(i,j,a,b) = t2(i,j,a,b) + 0.5*(t1(i,a)*t1(j,b)
                     - t1(i,b)*t1(j,a)),
    i,j occupied and a,b virtual spin orbitals.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_ccsd_spin_orb`

 
.. c:function:: void form_v(void)


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine form_v(nO,nV,t1,w,v)


    Builds V(i,j,k,a,b,c) for the (T) correction: W(i,j,k,a,b,c)
    plus the three terms coupling T1 to the vvoo integrals,
    v_vvoo(b,c,j,k)*t1(i,a) + v_vvoo(a,c,i,k)*t1(j,b) +
    v_vvoo(a,b,i,j)*t1(k,c). Integrals from the cc_space_* providers.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_vvoo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_par_t_space`

 
.. c:function:: void form_v_abc(void)


    File : :file:`ccsd_t_space_orb_abc.irp.f`

    .. code:: fortran

        subroutine form_v_abc(nO,nV,a,b,c,T_ov,X_oovv,W_abc,V_abc,W_cba,V_cba,W_bca,V_bca,W_cab,V_cab,W_bac,V_bac,W_acb,V_acb)


    Computes the six V(i,j,k) arrays for the six virtual permutations
    of one virtual triplet (a,b,c) of the (T) correction: each V is
    the corresponding W plus three terms coupling T1 (T_ov) to the
    vvoo integrals (X_oovv), e.g. V_abc(i,j,k) = W_abc(i,j,k) +
    X_oovv(j,k,b,c)*T_ov(i,a) + X_oovv(i,k,a,c)*T_ov(j,b) +
    X_oovv(i,j,a,b)*T_ov(k,c).

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_t_task_aba`
       * :c:func:`ccsd_t_task_abc`

 
.. c:function:: void form_v_ijk(void)


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine form_v_ijk(nO,nV,i,j,k,T_vo,X_vvoo,w,v)


    Computes the V(a,b,c) slice for one (i,j,k) ooo triple of the
    (T) correction: V(a,b,c) = W(a,b,c) + X_vvoo(b,c,k,j)*T_vo(a,i)
    + X_vvoo(a,c,k,i)*T_vo(b,j) + X_vvoo(a,b,j,i)*T_vo(c,k), OpenMP
    parallel over a,b,c.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_par_t_space_v2`

 
.. c:function:: void form_w(void)


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine form_w(nO,nV,t2,W)


    Builds W(i,j,k,a,b,c) for the (T) correction: sums T2 amplitudes
    against the vvvo integrals (d over virtuals, positive sign) and
    against the vooo integrals (l over occupied orbitals, negative
    sign), six index patterns each. Integrals from the cc_space_*
    module providers.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_vvvo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_par_t_space`

 
.. c:function:: void form_w_abc(void)


    File : :file:`ccsd_t_space_orb_abc.irp.f`

    .. code:: fortran

        subroutine form_w_abc(nO,nV,a,b,c,T_voov,T_oovv,X_vovv,X_ooov,W_abc,W_cba,W_bca,W_cab,W_bac,W_acb)


    Computes the six W(i,j,k) arrays, for the six virtual permutations
    abc, cba, bca, cab, bac, acb, of one virtual triplet (a,b,c) of
    the (T) correction, from the permuted T2 work arrays T_voov,
    T_oovv and the integral work arrays X_vovv, X_ooov, using BLAS
    dgemm.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_t_task_aba`
       * :c:func:`ccsd_t_task_abc`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void form_w_ijk(void)


    File : :file:`ccsd_t_space_orb.irp.f`

    .. code:: fortran

        subroutine form_w_ijk(nO,nV,i,j,k,T_vvoo,T_ovvo,X_vvvo,X_ovoo,W)


    Computes the W(a,b,c) slice for one (i,j,k) ooo triple of the
    (T) correction, i.e. W(i,j,k,::,::,::) of form_w, from the
    permuted T2 work arrays T_vvoo, T_ovvo and the integral work
    arrays X_vvvo, X_ovoo, using BLAS dgemm.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`ccsd_par_t_space_v2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void run_ccsd_space_orb(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine run_ccsd_space_orb


    Closed-shell CCSD in spatial orbitals (aborts otherwise), with
    the amplitudes and the two-electron integrals on the GPU.
    Iterates T1 and T2 with a DIIS or a simple update until the
    largest residual element is below cc_thresh_conv, saves the
    energies through save_energy, and evaluates the (T) correction
    with ccsd_par_t_space_stoch when cc_par_t is set.

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

 
.. c:function:: void run_ccsd_spin_orb(void)


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

 
.. c:function:: void save_energy(void)


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

 
.. c:function:: void update_tau_space(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine update_tau_space(nO,nV,h_t1,t1,t2,tau)


    Computes tau(i,j,a,b) = t2(i,j,a,b) + h_t1(j,b) * t1(i,a),
    with i, j occupied and a, b virtual spatial orbitals.

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

 
.. c:function:: void update_tau_x_space(void)


    File : :file:`ccsd_space_orb_sub.irp.f`

    .. code:: fortran

        subroutine update_tau_x_space(nO,nV,tau,tau_x)


    Computes tau_x(i,j,a,b) = 2 * tau(i,j,a,b) - tau(i,j,b,a).

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

 
.. c:function:: void use_cf_oo(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine use_cF_oo(nO,nV,t1,t2,tau_t,F_oo,F_ov,v_ooov,v_oovv,r1,r2)


    Apply the CCSD cF_oo intermediate to the spin-orbital
    residuals r1 and r2, modified in place:
    r1(i,a) =! r1(i,a) -( sum_m cF_oo(m,i)*t1(m,a),)
    r2(i,j,a,b) =! r2(i,j,a,b) +( sum_m cF_oo(m,i)*t2(j,m,a,b))
                    - cF_oo(m,j)*t2(i,m,a,b).

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_cf_oo`
       * :c:func:`dgemm`

 
.. c:function:: void use_cf_ov(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine use_cF_ov(nO,nV,t1,t2,F_ov,v_oovv,r1,r2)


    Apply the CCSD cF_ov intermediate to the spin-orbital
    residuals r1 and r2, modified in place:
    r1(i,a) =! r1(i,a) +( sum t2(i,m,a,e)*cF_ov(m,e),)
    r2(i,j,a,b) =! r2(i,j,a,b) -( 0.5*sum [t2(i,j,a,e)*t1(m,b))
      - t2(i,j,b,e)*t1(m,a) + t2(i,m,a,b)*t1(j,e)
      - t2(j,m,a,b)*t1(i,e)]*cF_ov(m,e).

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_cf_ov`
       * :c:func:`dgemm`

 
.. c:function:: void use_cf_vv(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine use_cF_vv(nO,nV,t1,t2,r1,r2)


    Intended to apply the cF_vv intermediate to the
    spin-orbital residuals r1 and r2. The compute_cF_vv
    call is commented out, so r1 and r2 are left unchanged.

 
.. c:function:: void use_cw_vvvf(void)


    File : :file:`ccsd_spin_orb_sub.irp.f`

    .. code:: fortran

        subroutine use_cW_vvvf(nO,nV,t1,t2,tau,v_oovv,r2)


    Apply the CCSD cW_vvvv intermediate to the spin-orbital
    residual r2, modified in place: for each f, the f-slice
    cW_vvvf(a,b,e) = cW_vvvv(a,b,e,f) is built by
    compute_cW_vvvf from the integral slices v_vvvf and v_vovf,
    and r2(i,j,a,b) =! and r2(i,j,a,b) +( 0.5*sum_e tau(i,j,e,f)*cW_vvvf(a,b,e).)

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

