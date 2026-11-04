.. _module_utils_cc: 
 
.. program:: utils_cc 
 
.. default-role:: option 
 
utils_cc
========

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
EZFIO parameters 
---------------- 
 
.. option:: cc_thresh_conv
 
    Threshold for the convergence of the residual equations.
 
    Default: 1e-6
 
.. option:: cc_max_iter
 
    Maximum number of iterations.
 
    Default: 100
 
.. option:: cc_diis_depth
 
    Maximum depth of the DIIS, i.e., maximum number of iterations that the DIIS keeps in memory. Warning, we allocate matrices with the diis depth at the beginning without update. If you don't have enough memory it should crash in memory.
 
    Default: 8
 
.. option:: cc_level_shift
 
    Level shift for the CC
 
    Default: 0.0
 
.. option:: cc_level_shift_guess
 
    Level shift for the guess of the CC amplitudes
 
    Default: 0.0
 
.. option:: cc_update_method
 
    Method used to update the CC amplitudes. none -> normal, diis -> with diis.
 
    Default: diis
 
.. option:: cc_guess_t1
 
    Guess used to initialize the T1 amplitudes. none -> 0, MP -> perturbation theory, read -> read from disk.
 
    Default: MP
 
.. option:: cc_guess_t2
 
    Guess used to initialize the T2 amplitudes. none -> 0, MP -> perturbation theory, read -> read from disk.
 
    Default: MP
 
.. option:: io_amplitudes
 
    Read/Write |CCSD| amplitudes from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: cc_par_t
 
    If true, the CCSD(T) will be computed.
 
    Default: False
 
.. option:: cc_par_t_stop
 
    Stops the calculation when the statistical error bar is below the given value.
 
    Default: 1.e-5
 
.. option:: cc_dev
 
    Only for dev purposes.
 
    Default: False
 
.. option:: cc_ref
 
    Index of the reference determinant in psi_det for CC calculation.
 
    Default: 1
 
 
Programs 
-------- 
 
 * :ref:`print_wf_qp_edit` 
 
Providers 
--------- 
 
.. c:var:: int cc_list_gen


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer, allocatable	:: cc_list_gen	(cc_n_mo)


    List of general orbitals without core and deleted ones

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`dim_list_core_orb`
       * :c:data:`dim_list_del_orb`
       * :c:data:`list_core`
       * :c:data:`list_del`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_ppqq`

 
.. c:var:: int cc_list_occ


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer, allocatable	:: cc_list_occ	(cc_nOa)
        integer, allocatable	:: cc_list_vir	(cc_nVa)


    List of occupied and virtual spatial orbitals without core and deleted ones

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`

 
.. c:var:: int cc_list_occ_spin


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer, allocatable	:: cc_list_occ_spin	(cc_nO_m,2)
        integer, allocatable	:: cc_list_vir_spin	(cc_nV_m,2)
        logical	:: cc_ref_is_open_shell	


    List of occupied and virtual spin orbitals without core and deleted ones

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_list_vir


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer, allocatable	:: cc_list_occ	(cc_nOa)
        integer, allocatable	:: cc_list_vir	(cc_nVa)


    List of occupied and virtual spatial orbitals without core and deleted ones

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`

 
.. c:var:: int cc_list_vir_spin


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer, allocatable	:: cc_list_occ_spin	(cc_nO_m,2)
        integer, allocatable	:: cc_list_vir_spin	(cc_nV_m,2)
        logical	:: cc_ref_is_open_shell	


    List of occupied and virtual spin orbitals without core and deleted ones

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_n_mo


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_no_m


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_no_s


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_noa


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_noab


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_nob


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_nv_m


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_nv_s


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_nva


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_nvab


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: int cc_nvb


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer	:: cc_no_m	
        integer	:: cc_noa	
        integer	:: cc_nob	
        integer	:: cc_noab	
        integer	:: cc_nv_m	
        integer	:: cc_nva	
        integer	:: cc_nvb	
        integer	:: cc_nvab	
        integer	:: cc_n_mo	
        integer, allocatable	:: cc_no_s	(2)
        integer, allocatable	:: cc_nv_s	(2)


    Number of orbitals without core and deleted ones of the cc_ref det in psi_det
    a: alpha, b: beta
    nO_m: max(a,b) occupied
    nOa: nb a occupied
    nOb: nb b occupied
    nOab: nb a+b occupied
    nV_m: max(a,b) virtual
    nVa: nb a virtual
    nVb: nb b virtual
    nVab: nb a+b virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_gen`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_space_f_o`
       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_v`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`
       * :c:data:`cc_space_v_aabb`
       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_aiia`
       * :c:data:`cc_space_v_iaai`
       * :c:data:`cc_space_v_iaia`
       * :c:data:`cc_space_v_iiaa`
       * :c:data:`cc_space_v_iijj`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`
       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: bool cc_ref_is_open_shell


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        integer, allocatable	:: cc_list_occ_spin	(cc_nO_m,2)
        integer, allocatable	:: cc_list_vir_spin	(cc_nV_m,2)
        logical	:: cc_ref_is_open_shell	


    List of occupied and virtual spin orbitals without core and deleted ones

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

 
.. c:var:: double cc_space_f_o


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_f_o	(cc_nOa)


    Diagonal of the OO Fock block of the reference
    determinant cc_ref: cc_space_f_o(i) = cc_space_f_oo(i,i)
    for i = 1,cc_nOa, i.e. F(i,i) for the occupied MOs
    (positions in cc_list_occ).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_space_f_oo`


 
.. c:var:: double cc_space_f_oo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_f_oo	(cc_nOa,cc_nOa)


    OO block of the alpha-spin Fock matrix of the reference
    determinant cc_ref: cc_space_f_oo(i1,i2) is the element
    F(idx1,idx2) with idx1 = cc_list_occ(i1) and idx2 =
    cc_list_occ(i2), extracted by gen_f_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_o`

 
.. c:var:: double cc_space_f_ov


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_f_ov	(cc_nOa,cc_nVa)


    OV block of the alpha-spin Fock matrix of the reference
    determinant cc_ref: cc_space_f_ov(i1,i2) is the element
    F(idx1,idx2) with idx1 = cc_list_occ(i1) and idx2 =
    cc_list_vir(i2), extracted by gen_f_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`


 
.. c:var:: double cc_space_f_v


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_f_v	(cc_nVa)


    Diagonal of the VV Fock block of the reference
    determinant cc_ref: cc_space_f_v(i) = cc_space_f_vv(i,i)
    for i = 1,cc_nVa, i.e. F(i,i) for the virtual MOs
    (positions in cc_list_vir).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_space_f_vv`


 
.. c:var:: double cc_space_f_vo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_f_vo	(cc_nVa,cc_nOa)


    VO block of the alpha-spin Fock matrix of the reference
    determinant cc_ref: cc_space_f_vo(i1,i2) is the element
    F(idx1,idx2) with idx1 = cc_list_vir(i1) and idx2 =
    cc_list_occ(i2), extracted by gen_f_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`


 
.. c:var:: double cc_space_f_vv


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_f_vv	(cc_nVa,cc_nVa)


    VV block of the alpha-spin Fock matrix of the reference
    determinant cc_ref: cc_space_f_vv(i1,i2) is the element
    F(idx1,idx2) with idx1 = cc_list_vir(i1) and idx2 =
    cc_list_vir(i2), extracted by gen_f_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_v`

 
.. c:var:: double cc_space_v


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v	(mo_num,mo_num,mo_num,mo_num)


    Two-electron MO integrals over the full MO space:
    cc_space_v(i,j,k,l) = <ij|kl> in MOs, i,j,k,l = 1,mo_num.
    If do_mo_cholesky is set, built from the Cholesky vectors
    (cholesky_mo_transp) with dgemm, otherwise from
    get_two_e_integral.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`banned_excitation`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_cache_shift`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`


 
.. c:var:: double cc_space_v_aabb


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_aabb	(cc_nVa,cc_nVa)


    <aa|bb> integrals
    a,b: virtual MO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_vvvv`


 
.. c:var:: double cc_space_v_aaii


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_aaii	(cc_nVa,cc_nOa)


    <aa|ii> integrals
    a: virtual MO
    i: occupied MO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_vvoo`


 
.. c:var:: double cc_space_v_aiia


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_aiia	(cc_nVa,cc_nOa)


    <ai|ia> integrals
    a: virtual MO
    i: inactive MO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_voov`


 
.. c:var:: double cc_space_v_iaai


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_iaai	(cc_nOa,cc_nVa)


    <ia|ai> integrals
    a: virtual MO
    i: inactive MO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_ovvo`


 
.. c:var:: double cc_space_v_iaia


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_iaia	(cc_nOa,cc_nVa)


    <ia|ia> integrals
    a: virtual MO
    i: occupied MO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_ovov`


 
.. c:var:: double cc_space_v_iiaa


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_iiaa	(cc_nOa,cc_nVa)


    <ii|aa> integrals
    a: virtual MO
    i: occupied MO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_oovv`


 
.. c:var:: double cc_space_v_iijj


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_iijj	(cc_nOa,cc_nOa)


    <ii|jj> integrals
    i,j: occupied MO

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_oooo`


 
.. c:var:: double cc_space_v_oo_chol


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_oo_chol	(cholesky_mo_num,cc_nOa,cc_nOa)


    Cholesky vectors of the OO block:
    cc_space_v_oo_chol(k,i1,i2) = cholesky_mo_transp(k,idx1,idx2)
    with idx1 = cc_list_occ(i1), idx2 = cc_list_occ(i2),
    k = 1,cholesky_mo_num, filled by gen_v_space_chol.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_vovo`

 
.. c:var:: double cc_space_v_oooo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_oooo	(cc_nOa,cc_nOa,cc_nOa,cc_nOa)


    Two-electron MO integrals of the OOOO block:
    cc_space_v_oooo(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1..i4 positions in cc_list_occ (occupied MOs, core and
    deleted ones excluded). If do_mo_cholesky is set, built from
    cc_space_v_oo_chol with dgemm, otherwise by gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_iijj`

 
.. c:var:: double cc_space_v_ooov


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_ooov	(cc_nOa,cc_nOa,cc_nOa,cc_nVa)


    OOOV block of two-electron MO integrals: cc_space_v_ooov
    (i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs, i4 in cc_list_vir,
    i1..i3 in cc_list_occ. Computed by gen_v_space, or copied
    from cc_space_v_ovoo if do_mo_cholesky is set; the loop
    bounds come from cc_space_v_ooov (once cc_space_v_oovo,
    wrong when cc_nVa != cc_nOa: bug, author-confirmed
    2026-10-31, now fixed).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`


 
.. c:var:: double cc_space_v_oovo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_oovo	(cc_nOa,cc_nOa,cc_nVa,cc_nOa)


    Two-electron MO integrals of the OOVO block:
    cc_space_v_oovo(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i3 in cc_list_vir, i1,i2,i4 in cc_list_occ. If
    do_mo_cholesky is set, copied from cc_space_v_vooo,
    otherwise computed by gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`


 
.. c:var:: double cc_space_v_oovv


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_oovv	(cc_nOa,cc_nOa,cc_nVa,cc_nVa)


    Two-electron MO integrals of the OOVV block:
    cc_space_v_oovv(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1,i2 in cc_list_occ, i3,i4 in cc_list_vir. If
    do_mo_cholesky is set, copied from
    cc_space_v_vvoo(i3,i4,i1,i2), otherwise computed by
    gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_iiaa`

 
.. c:var:: double cc_space_v_ov_chol


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_ov_chol	(cholesky_mo_num,cc_nOa,cc_nVa)


    Cholesky vectors of the OV block:
    cc_space_v_ov_chol(k,i1,i2) = cholesky_mo_transp(k,idx1,idx2)
    with idx1 = cc_list_occ(i1), idx2 = cc_list_vir(i2),
    k = 1,cholesky_mo_num, filled by gen_v_space_chol.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`


 
.. c:var:: double cc_space_v_ovoo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_ovoo	(cc_nOa,cc_nVa,cc_nOa,cc_nOa)


    Two-electron MO integrals of the OVOO block:
    cc_space_v_ovoo(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i2 in cc_list_vir, i1,i3,i4 in cc_list_occ. If
    do_mo_cholesky is set, copied from cc_space_v_vooo,
    otherwise computed by gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_ooov`

 
.. c:var:: double cc_space_v_ovov


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_ovov	(cc_nOa,cc_nVa,cc_nOa,cc_nVa)


    Two-electron MO integrals of the OVOV block:
    cc_space_v_ovov(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1,i3 in cc_list_occ, i2,i4 in cc_list_vir. If
    do_mo_cholesky is set, copied from
    cc_space_v_vovo(i2,i1,i4,i3), otherwise computed by
    gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_iaia`

 
.. c:var:: double cc_space_v_ovvo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_ovvo	(cc_nOa,cc_nVa,cc_nVa,cc_nOa)


    Two-electron MO integrals of the OVVO block:
    cc_space_v_ovvo(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i2,i3 in cc_list_vir, i1,i4 in cc_list_occ. If
    do_mo_cholesky is set, copied from
    cc_space_v_vvoo(i3,i2,i1,i4), otherwise computed by
    gen_v_space.
    Verified numerically against the direct integrals (H2O,
    cc-pVTZ, nOa=5, nVa=53, 1110 Cholesky vectors, 2026-10-31):
    the copy agrees to ~3e-14 relative. In the code's physicist
    convention <i1 i2 | i3 i4> = <i3 i2 | i1 i4> (8-fold
    symmetry, indices (i1,i3) and (i2,i4) on electrons 1,2).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_iaai`

 
.. c:var:: double cc_space_v_ovvv


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_ovvv	(cc_nOa,cc_nVa,cc_nVa,cc_nVa)


    Two-electron MO integrals of the OVVV block:
    cc_space_v_ovvv(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i2,i3,i4 in cc_list_vir, i1 in cc_list_occ. Computed
    by gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`


 
.. c:var:: double cc_space_v_ppqq


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_ppqq	(cc_n_mo,cc_n_mo)


    <pp|qq> integrals for general MOs (excepted core and deleted ones)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_gen`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`


 
.. c:var:: double cc_space_v_vo_chol


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_vo_chol	(cholesky_mo_num,cc_nVa,cc_nOa)


    Cholesky vectors of the VO block:
    cc_space_v_vo_chol(k,i1,i2) = cholesky_mo_transp(k,idx1,idx2)
    with idx1 = cc_list_vir(i1), idx2 = cc_list_occ(i2),
    k = 1,cholesky_mo_num, filled by gen_v_space_chol.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_vvoo`

 
.. c:var:: double cc_space_v_vooo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_vooo	(cc_nVa,cc_nOa,cc_nOa,cc_nOa)


    Two-electron MO integrals of the VOOO block:
    cc_space_v_vooo(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1 in cc_list_vir, i2..i4 in cc_list_occ. If do_mo_cholesky
    is set, built from cc_space_v_vo_chol and cc_space_v_oo_chol
    with dgemm, otherwise by gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_ovoo`

 
.. c:var:: double cc_space_v_voov


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_voov	(cc_nVa,cc_nOa,cc_nOa,cc_nVa)


    Two-electron MO integrals of the VOOV block:
    cc_space_v_voov(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1,i4 in cc_list_vir, i2,i3 in cc_list_occ. If
    do_mo_cholesky is set, copied from
    cc_space_v_vvoo(i1,i4,i3,i2), otherwise computed by
    gen_v_space.
    Verified numerically against the direct integrals (H2O,
    cc-pVTZ, nOa=5, nVa=53, 1110 Cholesky vectors, 2026-10-31):
    the copy agrees to ~3e-14 relative. In the code's physicist
    convention <i1 i2 | i3 i4> = <i1 i4 | i3 i2> (8-fold
    symmetry, indices (i1,i3) and (i2,i4) on electrons 1,2).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_aiia`

 
.. c:var:: double cc_space_v_vovo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_vovo	(cc_nVa,cc_nOa,cc_nVa,cc_nOa)


    Two-electron MO integrals of the VOVO block:
    cc_space_v_vovo(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1,i3 in cc_list_vir, i2,i4 in cc_list_occ. If
    do_mo_cholesky is set, built from cc_space_v_vv_chol
    and cc_space_v_oo_chol with dgemm, otherwise by
    gen_v_space.
    The dgemm contracts the Cholesky vectors over the
    (list1,list3) and (list2,list4) orbital pairs, T(k,i1,i3)
    and T(k,i2,i4), which reproduces <i1 i2 | i3 i4> in the
    code's physicist convention (= (i1 i2 | i3 i4) in chemist
    notation); the Cholesky branch agrees with the direct one,
    up to the Cholesky threshold (author, 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_ovov`

 
.. c:var:: double cc_space_v_vovv


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_vovv	(cc_nVa,cc_nOa,cc_nVa,cc_nVa)


    Two-electron MO integrals of the VOVV block:
    cc_space_v_vovv(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1,i3,i4 in cc_list_vir, i2 in cc_list_occ. Computed
    by gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`


 
.. c:var:: double cc_space_v_vv_chol


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_vv_chol	(cholesky_mo_num,cc_nVa,cc_nVa)


    Cholesky vectors of the VV block:
    cc_space_v_vv_chol(k,i1,i2) = cholesky_mo_transp(k,idx1,idx2)
    with idx1 = cc_list_vir(i1), idx2 = cc_list_vir(i2),
    k = 1,cholesky_mo_num, filled by gen_v_space_chol.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_vovo`

 
.. c:var:: double cc_space_v_vvoo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_vvoo	(cc_nVa,cc_nVa,cc_nOa,cc_nOa)


    Two-electron MO integrals of the VVOO block:
    cc_space_v_vvoo(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1,i2 in cc_list_vir, i3,i4 in cc_list_occ. If
    do_mo_cholesky is set, built from cc_space_v_vo_chol
    with dgemm, otherwise by gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_aaii`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_voov`

 
.. c:var:: double cc_space_v_vvov


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_vvov	(cc_nVa,cc_nVa,cc_nOa,cc_nVa)


    Two-electron MO integrals of the VVOV block:
    cc_space_v_vvov(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1,i2,i4 in cc_list_vir, i3 in cc_list_occ. Computed
    by gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`


 
.. c:var:: double cc_space_v_vvvo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_vvvo	(cc_nVa,cc_nVa,cc_nVa,cc_nOa)


    Two-electron MO integrals of the VVVO block:
    cc_space_v_vvvo(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1,i2,i3 in cc_list_vir, i4 in cc_list_occ. Computed
    by gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`


 
.. c:var:: double cc_space_v_vvvv


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_v_vvvv	(cc_nVa,cc_nVa,cc_nVa,cc_nVa)


    Two-electron MO integrals of the VVVV block:
    cc_space_v_vvvv(i1,i2,i3,i4) = <i1 i2 | i3 i4> in MOs,
    i1..i4 in cc_list_vir. Computed by gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_aabb`

 
.. c:var:: double cc_space_w_oovv


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_w_oovv	(cc_nOa,cc_nOa,cc_nVa,cc_nVa)


    Antisymmetrized OOVV block: cc_space_w_oovv(i,j,a,b) =
    2d0*<ij|ab> - <ji|ab>, antisymmetric in i,j. i,j are
    positions in cc_list_occ, a,b in cc_list_vir; <pq|rs>
    are the two-electron MO integrals from gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`


 
.. c:var:: double cc_space_w_vvoo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_space_w_vvoo	(cc_nVa,cc_nVa,cc_nOa,cc_nOa)


    Antisymmetrized VVOO block: cc_space_w_vvoo(a,b,i,j) =
    2d0*<ab|ij> - <ba|ij>, antisymmetric in a,b. a,b are
    positions in cc_list_vir, i,j in cc_list_occ; <pq|rs>
    are the two-electron MO integrals from gen_v_space.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cc_list_occ`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`


 
.. c:var:: double cc_spin_f_oo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_f_oo	(cc_nOab,cc_nOab)


    Spin-orbital OO block of the Fock matrix of the reference
    determinant cc_ref: p,q over the occupied spin orbitals,
    index order alpha then beta, cross-spin elements zero,
    as built by gen_f_spin.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`


 
.. c:var:: double cc_spin_f_ov


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_f_ov	(cc_nOab,cc_nVab)


    Spin-orbital OV block of the Fock matrix of the reference
    determinant cc_ref: p over the occupied and q over the
    virtual spin orbitals, index order alpha then beta,
    cross-spin elements zero, as built by gen_f_spin.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`


 
.. c:var:: double cc_spin_f_vv


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_f_vv	(cc_nVab,cc_nVab)


    Spin-orbital VV block of the Fock matrix of the reference
    determinant cc_ref: p,q over the virtual spin orbitals,
    index order alpha then beta, cross-spin elements zero,
    as built by gen_f_spin.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cc_ref`
       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`psi_det`


 
.. c:var:: double cc_spin_v_oo_chol


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_v_oo_chol	(cholesky_mo_num,cc_nOab,cc_nOab)


    Cholesky vectors of the OOOO block in spin orbitals:
    cc_spin_v_oo_chol(k,i1,i2) = cholesky_mo_transp(k,idx1,
    idx2) with idx1,idx2 positions in list_occ (occupied
    spin orbitals, alpha then beta), k = 1,cholesky_mo_num,
    filled by gen_v_space_chol.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`


 
.. c:var:: double cc_spin_v_oooo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_v_oooo	(cc_nOab,cc_nOab,cc_nOab,cc_nOab)


    Spin-orbital two-electron integrals of the OOOO block:
    cc_spin_v_oooo(p,q,r,s) = <pq||rs> in spin orbitals,
    p..s over the occupied spin orbitals, index order alpha
    then beta (1..cc_nOab), as built by gen_v_spin.
    The Cholesky path exists for open-shell CCSD as well: open-shell
    references are routed to run_ccsd_spin_orb (src/ccsd/ccsd.irp.f),
    whose driver loads cholesky_mo_transp when do_mo_cholesky is set
    (src/ccsd/ccsd_spin_orb_sub.irp.f), and gen_v_spin has a
    do_mo_cholesky branch for all four spin combinations
    (author-confirmed 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`


 
.. c:var:: double cc_spin_v_ooov


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_v_ooov	(cc_nOab,cc_nOab,cc_nOab,cc_nVab)


    Spin-orbital two-electron integrals of the OOOV block,
    <pq||rs> in spin orbitals, p,r,s occupied, q virtual,
    index order alpha then beta, as built by gen_v_spin.
    The Cholesky path exists for open-shell CCSD as well: open-shell
    references are routed to run_ccsd_spin_orb (src/ccsd/ccsd.irp.f),
    whose driver loads cholesky_mo_transp when do_mo_cholesky is set
    (src/ccsd/ccsd_spin_orb_sub.irp.f), and gen_v_spin has a
    do_mo_cholesky branch for all four spin combinations
    (author-confirmed 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`


 
.. c:var:: double cc_spin_v_oovo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_v_oovo	(cc_nOab,cc_nOab,cc_nVab,cc_nOab)


    Spin-orbital two-electron integrals of the OOVO block,
    <pq||rs> in spin orbitals, p,r,s occupied, q virtual,
    index order alpha then beta, as built by gen_v_spin.
    The Cholesky path exists for open-shell CCSD as well: open-shell
    references are routed to run_ccsd_spin_orb (src/ccsd/ccsd.irp.f),
    whose driver loads cholesky_mo_transp when do_mo_cholesky is set
    (src/ccsd/ccsd_spin_orb_sub.irp.f), and gen_v_spin has a
    do_mo_cholesky branch for all four spin combinations
    (author-confirmed 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`


 
.. c:var:: double cc_spin_v_oovv


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_v_oovv	(cc_nOab,cc_nOab,cc_nVab,cc_nVab)


    Spin-orbital two-electron integrals of the OOVV block,
    <pq||rs> in spin orbitals, p,r occupied, q,s virtual,
    index order alpha then beta, as built by gen_v_spin.
    The Cholesky path exists for open-shell CCSD as well: open-shell
    references are routed to run_ccsd_spin_orb (src/ccsd/ccsd.irp.f),
    whose driver loads cholesky_mo_transp when do_mo_cholesky is set
    (src/ccsd/ccsd_spin_orb_sub.irp.f), and gen_v_spin has a
    do_mo_cholesky branch for all four spin combinations
    (author-confirmed 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`


 
.. c:var:: double cc_spin_v_ovoo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_v_ovoo	(cc_nOab,cc_nVab,cc_nOab,cc_nOab)


    Spin-orbital two-electron integrals of the OVOO block,
    <pq||rs> in spin orbitals, p,r,s occupied, q virtual,
    index order alpha then beta, as built by gen_v_spin.
    The Cholesky path exists for open-shell CCSD as well: open-shell
    references are routed to run_ccsd_spin_orb (src/ccsd/ccsd.irp.f),
    whose driver loads cholesky_mo_transp when do_mo_cholesky is set
    (src/ccsd/ccsd_spin_orb_sub.irp.f), and gen_v_spin has a
    do_mo_cholesky branch for all four spin combinations
    (author-confirmed 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`


 
.. c:var:: double cc_spin_v_ovov


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_v_ovov	(cc_nOab,cc_nVab,cc_nOab,cc_nVab)


    Spin-orbital two-electron integrals of the OVOV block,
    <pq||rs> in spin orbitals, p,r occupied, q,s virtual,
    index order alpha then beta, as built by gen_v_spin.
    The Cholesky path exists for open-shell CCSD as well: open-shell
    references are routed to run_ccsd_spin_orb (src/ccsd/ccsd.irp.f),
    whose driver loads cholesky_mo_transp when do_mo_cholesky is set
    (src/ccsd/ccsd_spin_orb_sub.irp.f), and gen_v_spin has a
    do_mo_cholesky branch for all four spin combinations
    (author-confirmed 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`


 
.. c:var:: double cc_spin_v_ovvo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_v_ovvo	(cc_nOab,cc_nVab,cc_nVab,cc_nOab)


    Spin-orbital two-electron integrals of the OVVO block,
    <pq||rs> in spin orbitals, p,s occupied, q,r virtual,
    index order alpha then beta, as built by gen_v_spin.
    The Cholesky path exists for open-shell CCSD as well: open-shell
    references are routed to run_ccsd_spin_orb (src/ccsd/ccsd.irp.f),
    whose driver loads cholesky_mo_transp when do_mo_cholesky is set
    (src/ccsd/ccsd_spin_orb_sub.irp.f), and gen_v_spin has a
    do_mo_cholesky branch for all four spin combinations
    (author-confirmed 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`


 
.. c:var:: double cc_spin_v_vvoo


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cc_spin_v_vvoo	(cc_nVab,cc_nVab,cc_nOab,cc_nOab)


    Spin-orbital two-electron integrals of the VVOO block,
    <pq||rs> in spin orbitals, p,q virtual, r,s occupied,
    index order alpha then beta, as built by gen_v_spin.
    The Cholesky path exists for open-shell CCSD as well: open-shell
    references are routed to run_ccsd_spin_orb (src/ccsd/ccsd.irp.f),
    whose driver loads cholesky_mo_transp when do_mo_cholesky is set
    (src/ccsd/ccsd_spin_orb_sub.irp.f), and gen_v_spin has a
    do_mo_cholesky branch for all four spin combinations
    (author-confirmed 2026-10-31).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`
       * :c:data:`cc_no_m`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void compute_err1(void)


    File : :file:`diis.irp.f`

    .. code:: fortran

        subroutine compute_err1(nO,nV,f_o,f_v,r1,err1)


    Compute the error vector for the t1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_level_shift`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`update_t_ccsd_diis`
       * :c:func:`update_t_ccsd_diis_v3`

 
.. c:function:: void compute_err2(void)


    File : :file:`diis.irp.f`

    .. code:: fortran

        subroutine compute_err2(nO,nV,f_o,f_v,r2,err2)


    Compute the error vector for the t2

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_level_shift`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`update_t_ccsd`
       * :c:func:`update_t_ccsd_diis`
       * :c:func:`update_t_ccsd_diis_v3`

 
.. c:function:: void det_energy(void)


    File : :file:`energy.irp.f`

    .. code:: fortran

        subroutine det_energy(det,energy)


    Energy of the input determinant: the diagonal Hamiltonian matrix
    element $\langle i|H|i \rangle$ plus the nuclear repulsion.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`nuclear_repulsion`

 
.. c:function:: void diis_cc(void)


    File : :file:`diis.irp.f`

    .. code:: fortran

        subroutine diis_cc(all_err,all_t,sze,m,iter,t)


    DIIS. Take the error vectors and the amplitudes of the previous
    iterations to compute the new amplitudes

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`update_t_ccsd`
       * :c:func:`update_t_ccsd_diis`
       * :c:func:`update_t_ccsd_diis_v3`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`abort`
       * :c:func:`dgesv`

 
.. c:function:: void extract_list_orb_space(void)


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        subroutine extract_list_orb_space(det,nO,nV,list_occ,list_vir)


    Returns the the list of occupied and virtual alpha spin orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`abort`
       * :c:func:`apply_hole`

 
.. c:function:: void extract_list_orb_spin(void)


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        subroutine extract_list_orb_spin(det,nO_m,nV_m,list_occ,list_vir)


    Returns the the list of occupied alpha/beta, virtual alpha/beta spin orbitals
    size(nO_m,1) must be max(nOa,nOb) and size(nV_m,1) must be max(nVa,nVb)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_list_occ_spin`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_hole`

 
.. c:function:: void extract_n_spin(void)


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        subroutine extract_n_spin(det,n)


    Returns the number of occupied alpha, occupied beta, virtual alpha, virtual beta spin orbitals
    in det without counting the core and deleted orbitals in the format n(nOa,nOb,nVa,nVb)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_no_m`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_hole`

 
.. c:function:: void gen_f_space(void)


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        subroutine gen_f_space(det,n1,n2,list1,list2,f)


    Extract a block of the Fock matrix of determinant det.
    The alpha-spin Fock matrix is computed by
    get_fock_matrix_spin, and f(i1,i2) is the element
    F(list1(i1),list2(i2)) of it.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_f_oo`
       * :c:data:`cc_space_f_ov`
       * :c:data:`cc_space_f_vo`
       * :c:data:`cc_space_f_vv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_fock_matrix_spin`

 
.. c:function:: void gen_f_spin(void)


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        subroutine gen_f_spin(det, n1,n2, n1_S,n2_S, list1,list2, dim1,dim2, f)


    Compute the Fock matrix corresponding to two lists of spin orbitals.
    Ex: occ/occ, occ/vir,...

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_spin_f_oo`
       * :c:data:`cc_spin_f_ov`
       * :c:data:`cc_spin_f_vv`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_fock_matrix_spin`
       * :c:func:`shift_idx_spin`

 
.. c:function:: void gen_v_space(void)


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        subroutine gen_v_space(n1,n2,n3,n4,list1,list2,list3,list4,v)


    Extract two-electron MO integrals over four lists of MOs:
    v(i1,i2,i3,i4) = <idx1 idx2 | idx3 idx4> in MOs, with
    idx1 = list1(i1), idx2 = list2(i2), idx3 = list3(i3) and
    idx4 = list4(i4). If do_mo_cholesky is set, the block is
    built from the Cholesky vectors (gen_v_space_chol),
    otherwise from get_two_e_integral.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`
       * :c:data:`mo_integrals_map`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_oooo`
       * :c:data:`cc_space_v_ooov`
       * :c:data:`cc_space_v_oovo`
       * :c:data:`cc_space_v_oovv`
       * :c:data:`cc_space_v_ovoo`
       * :c:data:`cc_space_v_ovov`
       * :c:data:`cc_space_v_ovvo`
       * :c:data:`cc_space_v_ovvv`
       * :c:data:`cc_space_v_ppqq`
       * :c:data:`cc_space_v_vooo`
       * :c:data:`cc_space_v_voov`
       * :c:data:`cc_space_v_vovo`
       * :c:data:`cc_space_v_vovv`
       * :c:data:`cc_space_v_vvoo`
       * :c:data:`cc_space_v_vvov`
       * :c:data:`cc_space_v_vvvo`
       * :c:data:`cc_space_v_vvvv`
       * :c:data:`cc_space_w_oovv`
       * :c:data:`cc_space_w_vvoo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`gen_v_space_chol`

 
.. c:function:: void gen_v_space_chol(void)


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        subroutine gen_v_space_chol(n1,n3,list1,list3,v,ldv)


    Fill a buffer with Cholesky vectors: v(k,i1,i3) =
    cholesky_mo_transp(k,idx1,idx3), with idx1 = list1(i1),
    idx3 = list3(i3) and k = 1,ldv. The Cholesky index k is
    the first index (leading dimension ldv) so that the buffer
    can be used directly in a dgemm call.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v_oo_chol`
       * :c:data:`cc_space_v_ov_chol`
       * :c:data:`cc_space_v_vo_chol`
       * :c:data:`cc_space_v_vv_chol`
       * :c:data:`cc_spin_v_oo_chol`
       * :c:func:`gen_v_space`
       * :c:func:`gen_v_spin`
       * :c:func:`gen_v_spin_3idx`
       * :c:func:`gen_v_spin_3idx_i_kl`
       * :c:func:`gen_v_spin_3idx_ij_l`

 
.. c:function:: void gen_v_spin(void)


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        subroutine gen_v_spin(n1,n2,n3,n4, n1_S,n2_S,n3_S,n4_S, list1,list2,list3,list4, dim1,dim2,dim3,dim4, v)


    Compute the 2e-integrals corresponding to four lists of spin orbitals.
    Ex: occ/occ/occ/occ, occ/vir/occ/vir, ...

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`cc_spin_v_oooo`
       * :c:data:`cc_spin_v_ooov`
       * :c:data:`cc_spin_v_oovo`
       * :c:data:`cc_spin_v_oovv`
       * :c:data:`cc_spin_v_ovoo`
       * :c:data:`cc_spin_v_ovov`
       * :c:data:`cc_spin_v_ovvo`
       * :c:data:`cc_spin_v_vvoo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`gen_v_space_chol`
       * :c:func:`shift_idx_spin`

 
.. c:function:: void gen_v_spin_3idx(void)


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        subroutine gen_v_spin_3idx(n1,n2,n3,n4, idx_l, n1_S,n2_S,n3_S,n4_S, list1,list2,list3,list4, dim1,dim2,dim3, v_l)


    Compute the bi electronic integrals corresponding to four lists of spin orbitals.
    Ex: occ/occ/occ/occ, occ/vir/occ/vir, ...

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`gen_v_space_chol`
       * :c:func:`shift_idx_spin`

 
.. c:function:: void gen_v_spin_3idx_i_kl(void)


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        subroutine gen_v_spin_3idx_i_kl(n1,n2,n3,n4, idx_j, n1_S,n2_S,n3_S,n4_S, list1,list2,list3,list4, dim1,dim2,dim3, v_j)


    Compute the bi electronic integrals corresponding to four lists of spin orbitals.
    Ex: occ/occ/occ/occ, occ/vir/occ/vir, ...

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`gen_v_space_chol`
       * :c:func:`shift_idx_spin`

 
.. c:function:: void gen_v_spin_3idx_ij_l(void)


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        subroutine gen_v_spin_3idx_ij_l(n1,n2,n3,n4, idx_k, n1_S,n2_S,n3_S,n4_S, list1,list2,list3,list4, dim1,dim2,dim3, v_k)


    Compute the bi electronic integrals corresponding to four lists of spin orbitals.
    Ex: occ/occ/occ/occ, occ/vir/occ/vir, ...

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_space_v`
       * :c:data:`cholesky_mo_num`
       * :c:data:`do_mo_cholesky`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`gen_v_space_chol`
       * :c:func:`shift_idx_spin`

 
.. c:function:: void get_excitation_general(void)


    File : :file:`phase.irp.f`

    .. code:: fortran

        subroutine get_excitation_general(det1,det2,degree,n,list_anni,list_crea,phase,Nint)


    Excitation that maps det1 onto det2. list_anni (list_crea) are
    the 1-based orbital indices, increasing per spin, of the orbitals
    occupied in det1 (det2) but not in det2 (det1). n and degree are
    the excitation degrees per spin and total. phase is $(-1)^d$, d
    the sum, per spin, of the distance between the i-th annihilation
    and creation positions in the occupied lists of det1 (det2).

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_phase_general`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`abort`
       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: void get_fock_matrix_spin(void)


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        subroutine get_fock_matrix_spin(det,s,f)


    Fock matrix alpha or beta of an arbitrary det

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gen_f_space`
       * :c:func:`gen_f_spin`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`apply_hole`

 
.. c:function:: void get_phase_general(void)


    File : :file:`phase.irp.f`

    .. code:: fortran

        subroutine get_phase_general(det1,det2,phase,degree,Nint)


    Phase and degree of the excitation that maps det1 onto det2.
    Thin wrapper around get_excitation_general: the orbital lists
    are discarded, only phase and degree are returned. det1 and
    det2 must have the same number of electrons per spin,
    otherwise the subroutine aborts.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_general`

 
.. c:function:: void guess_t1(void)


    File : :file:`guess_t.irp.f`

    .. code:: fortran

        subroutine guess_t1(nO,nV,f_o,f_v,f_ov,t1)


    Update the T1 amplitudes for CC

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_guess_t1`
       * :c:data:`cc_level_shift_guess`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`abort`
       * :c:func:`read_t1`

 
.. c:function:: void guess_t2(void)


    File : :file:`guess_t.irp.f`

    .. code:: fortran

        subroutine guess_t2(nO,nV,f_o,f_v,v_oovv,t2)


    Update the T2 amplitudes for CC

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_guess_t2`
       * :c:data:`cc_level_shift_guess`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`abort`
       * :c:func:`read_t2`

 
.. c:function:: void is_core(void)


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        function is_core(i)


    True if the orbital i is a core orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`list_core`

 
.. c:function:: void is_del(void)


    File : :file:`occupancy.irp.f`

    .. code:: fortran

        function is_del(i)


    True if the orbital i is a deleted orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_del_orb`
       * :c:data:`list_del`

 
.. c:function:: void print_wf_qp_edit_subroutine(void)


    File : :file:`print_wf_qp_edit.irp.f`

    .. code:: fortran

        subroutine print_wf_qp_edit_subroutine()


    Prints the first n_det_qp_edit determinants of the wave
    function stored in psi_coef and psi_det: for each determinant
    its index, its coefficients psi_coef(i,:) and the determinant
    itself, via print_det.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_qp_edit`
       * :c:data:`n_int`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_wf_qp_edit`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_det`

 
.. c:function:: void read_t1(void)


    File : :file:`guess_t.irp.f`

    .. code:: fortran

        subroutine read_t1(nO,nV,t1)


    Read the T1 amplitudes for CC

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`read_amplitudes`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`guess_t1`

 
.. c:function:: void read_t2(void)


    File : :file:`guess_t.irp.f`

    .. code:: fortran

        subroutine read_t2(nO,nV,t2)


    Read the T2 amplitudes for CC

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`read_amplitudes`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`guess_t2`

 
.. c:function:: void shift_idx_spin(void)


    File : :file:`mo_integrals_cc.irp.f`

    .. code:: fortran

        subroutine shift_idx_spin(s,n_S,shift)


    Shift for the partitionning alpha/beta of the spin orbitals
    n_S(1): number of spin alpha in the correspondong list
    n_S(2): number of spin beta in the correspondong list

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`gen_f_spin`
       * :c:func:`gen_v_spin`
       * :c:func:`gen_v_spin_3idx`
       * :c:func:`gen_v_spin_3idx_i_kl`
       * :c:func:`gen_v_spin_3idx_ij_l`

 
.. c:function:: void update_all_err(void)


    File : :file:`diis.irp.f`

    .. code:: fortran

        subroutine update_all_err(err,all_err,sze,m,iter)


    Shift all the err vectors of the previous iterations to add the new one
    The last err vector is placed in the last position and all the others are
    moved toward the first one.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`update_t_ccsd`
       * :c:func:`update_t_ccsd_diis`
       * :c:func:`update_t_ccsd_diis_v3`

 
.. c:function:: void update_all_t(void)


    File : :file:`diis.irp.f`

    .. code:: fortran

        subroutine update_all_t(t,all_t,sze,m,iter)


    Shift all the t vectors of the previous iterations to add the new one
    The last t vector is placed in the last position and all the others are
    moved toward the first one.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`update_t_ccsd`
       * :c:func:`update_t_ccsd_diis`
       * :c:func:`update_t_ccsd_diis_v3`

 
.. c:function:: void update_t1(void)


    File : :file:`update_t.irp.f`

    .. code:: fortran

        subroutine update_t1(nO,nV,f_o,f_v,r1,t1)


    Update the T1 amplitudes for CC

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_level_shift`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`update_t_ccsd`

 
.. c:function:: void update_t2(void)


    File : :file:`update_t.irp.f`

    .. code:: fortran

        subroutine update_t2(nO,nV,f_o,f_v,r2,t2)


    Update the T2 amplitudes for CC

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_level_shift`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`update_t_ccsd`

 
.. c:function:: void update_t_ccsd(void)


    File : :file:`diis.irp.f`

    .. code:: fortran

        subroutine update_t_ccsd(nO,nV,nb_iter,f_o,f_v,r1,r2,t1,t2,all_err1,all_err2,all_t1,all_t2)


    Update the CCSD amplitudes t1 and t2 in place from the residuals
    r1 and r2. If cc_update_method is 'diis', t1 gets the standard
    update and t2 is extrapolated with DIIS on the histories all_err2
    and all_t2. If it is 'none', both get the standard update. Any
    other value prints a message and leaves t1 and t2 unchanged.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_diis_depth`
       * :c:data:`cc_update_method`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_err2`
       * :c:func:`diis_cc`
       * :c:func:`update_all_err`
       * :c:func:`update_all_t`
       * :c:func:`update_t1`
       * :c:func:`update_t2`

 
.. c:function:: void update_t_ccsd_diis(void)


    File : :file:`diis.irp.f`

    .. code:: fortran

        subroutine update_t_ccsd_diis(nO,nV,nb_iter,f_o,f_v,r1,r2,t1,t2,all_err1,all_err2,all_t1,all_t2)


    Update the CCSD amplitudes t1 and t2 in place with a single DIIS
    extrapolation on the concatenated t1/t2 amplitude vector. The
    error vectors computed from r1 and r2, and the amplitudes, are
    first appended to the per-amplitude histories all_err1, all_err2,
    all_t1 and all_t2, then the combined histories are extrapolated
    and the result is split back into t1 and t2.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_diis_depth`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_err1`
       * :c:func:`compute_err2`
       * :c:func:`diis_cc`
       * :c:func:`update_all_err`
       * :c:func:`update_all_t`

 
.. c:function:: void update_t_ccsd_diis_v3(void)


    File : :file:`diis.irp.f`

    .. code:: fortran

        subroutine update_t_ccsd_diis_v3(nO,nV,nb_iter,f_o,f_v,r1,r2,t1,t2,all_err,all_t)


    Update the CCSD amplitudes t1 and t2 in place with DIIS on the
    concatenated t1/t2 amplitude vector, like update_t_ccsd_diis,
    but with a single history buffer: all_err and all_t already hold
    the combined t1/t2 errors and amplitudes, and t1 and t2 are flat
    vectors.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cc_diis_depth`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`compute_err1`
       * :c:func:`compute_err2`
       * :c:func:`diis_cc`
       * :c:func:`update_all_err`
       * :c:func:`update_all_t`

 
.. c:function:: void write_t1(void)


    File : :file:`guess_t.irp.f`

    .. code:: fortran

        subroutine write_t1(nO,nV,t1)


    Write the T1 amplitudes for CC

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`read_amplitudes`

 
.. c:function:: void write_t2(void)


    File : :file:`guess_t.irp.f`

    .. code:: fortran

        subroutine write_t2(nO,nV,t2)


    Write the T2 amplitudes for CC

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`read_amplitudes`

