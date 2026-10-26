.. _module_casscf_cipsi: 
 
.. program:: casscf_cipsi 
 
.. default-role:: option 
 
======
casscf
======

|CASSCF| program with the CIPSI algorithm.


Example of inputs for GROUND STATE calculations
-----------------------------------------------
NOTICE :: FOR EXCITED STATES CALCULATIONS SEE THE FILE "example_casscf_multistate.sh"

a) Small active space : standard CASSCF 
---------------------------------------
Let's do O2 (triplet) in aug-cc-pvdz with the following geometry (xyz format, Bohr units)
2

 O           0.0000000000        0.0000000000       -1.1408000000
 O           0.0000000000        0.0000000000        1.1408000000

# Create the ezfio folder 
qp create_ezfio -b aug-cc-pvdz O2.xyz -m 3 -a -o O2_avdz

# Start with an ROHF guess 
qp run scf | tee ${EZFIO_FILE}.rohf.out

# Get the ROHF energy for check 
qp get hartree_fock energy # should be -149.4684509

# Define the full valence active space: the two 1s are doubly occupied, the other 8 valence orbitals are active 
# CASSCF(12e,10orb) 
qp set_mo_class -c "[1-2]" -a "[3-10]" -v "[11-46]"

# Specify that you want an near exact CASSCF, i.e. the CIPSI selection will stop at pt2_max = 10^-10
qp set casscf_cipsi small_active_space True 
# RUN THE CASSCF 
qp run casscf | tee ${EZFIO_FILE}.casscf.out
# you should find around -149.7243542


b) Large active space : Exploit the selected CI in the active space 
-------------------------------------------------------------------
#Let us start from the small active space calculation orbitals and add another 10 virtuals: CASSCF(12e,20orb)
qp set_mo_class -c "[1-2]" -a "[3-20]" -v "[21-46]"
# As this active space is larger, you unset the small_active_space feature 
qp set casscf_cipsi small_active_space False
# As it is a large active space, the energy convergence thereshold is set to be 0.0001
qp run casscf | tee ${EZFIO_FILE}.casscf_large.out
# you should find around -149.9046


 
 
 
EZFIO parameters 
---------------- 
 
.. option:: energy
 
    Calculated Selected |FCI| energy
 
 
.. option:: energy_pt2
 
    Calculated |FCI| energy + |PT2|
 
 
.. option:: state_following_casscf
 
    If |true|, the CASSCF will try to follow the guess CI vector and orbitals
 
    Default: False
 
.. option:: diag_hess_cas
 
    If |true|, only the DIAGONAL part of the hessian is retained for the CASSCF
 
    Default: False
 
.. option:: hess_cv_cv
 
    If |true|, the core-virtual - core-virtual part of the hessian is computed
 
    Default: True
 
.. option:: level_shift_casscf
 
    Energy shift on the virtual MOs to improve SCF convergence
 
    Default: 0.005
 
.. option:: fast_2rdm
 
    If true, the two-rdm are computed with a fast algo
 
    Default: True
 
.. option:: criterion_casscf
 
    choice of the criterion for the convergence of the casscf: can be  energy  or  gradients  or e_pt2
 
    Default: e_pt2
 
.. option:: thresh_casscf
 
    Threshold on the convergence of the CASCF energy.
 
    Default: 1.e-06
 
.. option:: pt2_min_casscf
 
    Minimum value of the pt2_max parameter for the CIPSI in the CASSCF iterations.
 
    Default: 1.e-04
 
.. option:: n_big_act_orb
 
    Number of active orbitals from which the active space is considered as large, and therefore pt2_min_casscf is activated.
 
    Default: 16
 
.. option:: adaptive_pt2_max
 
    If |true|, the pt2_max value in the CIPSI iterations will automatically adapt, otherwise it is fixed at the value given in the EZFIO folder
 
    Default: True
 
.. option:: small_active_space
 
    If |true|, the pt2_max value in the CIPSI is set to 10-10 and will not change
 
    Default: False
 
.. option:: act_mos_opt
 
    If |true|, the active orbitals are also optimized variationally
 
    Default: False
 
 
Programs 
-------- 
 
 * :ref:`casscf` 
 * :ref:`get_energy` 
 * :ref:`test_chol` 
 
Providers 
--------- 
 
.. c:var:: bool bavard


    File : :file:`bavard.irp.f`

    .. code:: fortran

        logical	:: bavard	



    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`best_vector_ovrlp_casscf`
       * :c:data:`excit`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`gradvec2`
       * :c:data:`gradvec_old`
       * :c:data:`hessmat_old`
       * :c:data:`hessmat_peter`
       * :c:data:`natorbsci`
       * :c:data:`occnum`
       * :c:data:`p0tuvx`
       * :c:data:`p0tuvx_peter`
       * :c:data:`sxeigenvec`
       * :c:data:`sxmatrix`
       * :c:data:`umat`

 
.. c:var:: double best_overlap_casscf


    File : :file:`neworbs.irp.f`

    .. code:: fortran

        integer	:: best_vector_ovrlp_casscf	
        double precision	:: best_overlap_casscf	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`nmonoex`
       * :c:data:`sxeigenvec`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_improvement`
       * :c:data:`neworbs`
       * :c:data:`sxvector`

 
.. c:var:: int best_vector_ovrlp_casscf


    File : :file:`neworbs.irp.f`

    .. code:: fortran

        integer	:: best_vector_ovrlp_casscf	
        double precision	:: best_overlap_casscf	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`nmonoex`
       * :c:data:`sxeigenvec`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`energy_improvement`
       * :c:data:`neworbs`
       * :c:data:`sxvector`

 
.. c:var:: double bielec_pqxx_array


    File : :file:`bielec.irp.f`

    .. code:: fortran

        real*8, allocatable	:: bielec_pqxx_array	(mo_num,mo_num,n_core_inact_act_orb,n_core_inact_act_orb)


    WARNING !!! Old version !!! NOT USED ANYMORE IN THE PROGRAM !!! TOO BIG TO BE STORED ON LARGE SYSTEMS !!!
    
    Replaced by the Cholesky-based function bielec_PQxx
    
    bielec_PQxx_array : integral (pq|xx) with p,q arbitrary, x core or active
    indices are unshifted orbital numbers

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`gpu_num`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`nproc_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_no_array`

 
.. c:var:: double bielec_pqxx_no_array


    File : :file:`bielec_natorb.irp.f`

    .. code:: fortran

        real*8, allocatable	:: bielec_pqxx_no_array	(mo_num,mo_num,n_core_inact_act_orb,n_core_inact_act_orb)


    WARNING !!! Old version !!! NOT USED ANYMORE IN THE PROGRAM !!! TOO BIG TO BE STORED ON LARGE SYSTEMS !!!
    
    Replaced by the Cholesky-based function bielec_PQxx_no
    
    integral (pq|xx) in the basis of natural MOs
    indices are unshifted orbital numbers
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`natorbsci`


 
.. c:var:: double bielec_pxxq_array


    File : :file:`bielec.irp.f`

    .. code:: fortran

        real*8, allocatable	:: bielec_pxxq_array	(mo_num,n_core_inact_act_orb,n_core_inact_act_orb,mo_num)


    WARNING !!! Old version !!! NOT USED ANYMORE IN THE PROGRAM !!! TOO BIG TO BE STORED ON LARGE SYSTEMS !!!
    
    Replaced by the Cholesky-based function bielec_PxxQ
    
    bielec_PxxQ_array : integral (px|xq) with p,q arbitrary, x core or active
    indices are unshifted orbital numbers

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_mo_transp_sp`
       * :c:data:`cholesky_mo_transp_sp_d`
       * :c:data:`do_mo_cholesky`
       * :c:data:`gpu_busy`
       * :c:data:`gpu_integral_buffer1`
       * :c:data:`gpu_num`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`mo_cholesky_double`
       * :c:data:`mo_integrals_cache_min`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`nproc_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pxxq_no_array`

 
.. c:var:: double bielec_pxxq_no_array


    File : :file:`bielec_natorb.irp.f`

    .. code:: fortran

        real*8, allocatable	:: bielec_pxxq_no_array	(mo_num,n_core_inact_act_orb,n_core_inact_act_orb,mo_num)


    WARNING !!! Old version !!! NOT USED ANYMORE IN THE PROGRAM !!! TOO BIG TO BE STORED ON LARGE SYSTEMS !!!
    
    Replaced by the Cholesky-based function bielec_PxxQ_no
    
    integral (px|xq) in the basis of natural MOs
    indices are unshifted orbital numbers

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pxxq_array`
       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`natorbsci`


 
.. c:var:: double bielecci


    File : :file:`bielec.irp.f`

    .. code:: fortran

        real*8, allocatable	:: bielecci	(n_act_orb,n_act_orb,n_act_orb,mo_num)


    bielecCI : integrals (tu|vp) with p arbitrary, tuv active
    index p runs over the whole basis, t,u,v only over the active orbitals
    
    This array can be stored anyway. Ex: 50 active orbitals, 1500 MOs ==> 8x50^3x1500 = 1.5 Gb

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`list_act`
       * :c:data:`mo_integrals_map`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci_no`
       * :c:data:`etwo`

 
.. c:var:: double bielecci_no


    File : :file:`bielec_natorb.irp.f`

    .. code:: fortran

        real*8, allocatable	:: bielecci_no	(n_act_orb,n_act_orb,n_act_orb,mo_num)


    integrals (tu|vp) in the basis of natural MOs
    index p runs over the whole basis, t,u,v only over the active orbitals
    
    This array can be stored anyway. Ex: 50 active orbitals, 1500 MOs ==> 8x50^3x1500 = 1.5 Gb

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci`
       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`natorbsci`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`

 
.. c:var:: double cholesky_no_1_idx_transp


    File : :file:`chol_bielec.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cholesky_no_1_idx_transp	(cholesky_mo_num,n_act_orb,mo_num)


    Cholesky vectors with ONE orbital on the active natural orbital basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`natorbsci`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`cholesky_no_total_transp`

 
.. c:var:: double cholesky_no_2_idx_transp


    File : :file:`chol_bielec.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cholesky_no_2_idx_transp	(cholesky_mo_num,n_act_orb,n_act_orb)


    Cholesky vectors with TWO orbital on the active natural orbital basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`natorbsci`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_no_total_transp`

 
.. c:var:: double cholesky_no_total_transp


    File : :file:`chol_bielec.irp.f`

    .. code:: fortran

        double precision, allocatable	:: cholesky_no_total_transp	(cholesky_mo_num,mo_num,mo_num)


    Cholesky vectors defined on all basis including the NO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_mo_transp`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`

 
.. c:var:: double d0tu


    File : :file:`densities.irp.f`

    .. code:: fortran

        real*8, allocatable	:: d0tu	(n_act_orb,n_act_orb)


    the first-order density matrix in the basis of the starting MOs.
    matrix is state averaged.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`one_e_dm_mo_alpha_average`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`etwo`
       * :c:data:`natorbsci`

 
.. c:var:: double d0tu_alpha_ao


    File : :file:`densities.irp.f`

    .. code:: fortran

        double precision, allocatable	:: d0tu_alpha_ao	(ao_num,ao_num)
        double precision, allocatable	:: d0tu_beta_ao	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`one_e_dm_mo_alpha_average`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mcscf_fock_alpha_ao`

 
.. c:var:: double d0tu_beta_ao


    File : :file:`densities.irp.f`

    .. code:: fortran

        double precision, allocatable	:: d0tu_alpha_ao	(ao_num,ao_num)
        double precision, allocatable	:: d0tu_beta_ao	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`one_e_dm_mo_alpha_average`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mcscf_fock_alpha_ao`

 
.. c:var:: int dim_n_max_overlap


    File : :file:`swap_orb.irp.f`

    .. code:: fortran

        integer, allocatable	:: max_overlap	(nMonoEx)
        integer	:: n_max_overlap	
        integer	:: dim_n_max_overlap	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nmonoex`
       * :c:data:`sxvector_lowest`
       * :c:data:`thresh_overlap_switch`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`orb_swap`

 
.. c:var:: bool do_ddci


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the CAS case, all those are always false except do_only_cas


 
.. c:var:: bool do_only_1h1p


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the CAS case, all those are always false except do_only_cas


 
.. c:var:: bool do_only_cas


    File : :file:`class.irp.f`

    .. code:: fortran

        logical	:: do_only_1h1p	
        logical	:: do_only_cas	
        logical	:: do_ddci	


    In the CAS case, all those are always false except do_only_cas


 
.. c:var:: double ecore


    File : :file:`tot_en.irp.f`

    .. code:: fortran

        real*8	:: etwo	
        real*8	:: eone	
        real*8	:: eone_bis	
        real*8	:: etwo_bis	
        real*8	:: etwo_ter	
        real*8	:: ecore	
        real*8	:: ecore_bis	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci`
       * :c:data:`d0tu`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`nuclear_repulsion`
       * :c:data:`p0tuvx`


 
.. c:var:: double ecore_bis


    File : :file:`tot_en.irp.f`

    .. code:: fortran

        real*8	:: etwo	
        real*8	:: eone	
        real*8	:: eone_bis	
        real*8	:: etwo_bis	
        real*8	:: etwo_ter	
        real*8	:: ecore	
        real*8	:: ecore_bis	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci`
       * :c:data:`d0tu`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`nuclear_repulsion`
       * :c:data:`p0tuvx`


 
.. c:var:: double energy_improvement


    File : :file:`neworbs.irp.f`

    .. code:: fortran

        real*8	:: energy_improvement	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`best_vector_ovrlp_casscf`
       * :c:data:`state_following_casscf`
       * :c:data:`sxeigenvec`


 
.. c:var:: double eone


    File : :file:`tot_en.irp.f`

    .. code:: fortran

        real*8	:: etwo	
        real*8	:: eone	
        real*8	:: eone_bis	
        real*8	:: etwo_bis	
        real*8	:: etwo_ter	
        real*8	:: ecore	
        real*8	:: ecore_bis	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci`
       * :c:data:`d0tu`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`nuclear_repulsion`
       * :c:data:`p0tuvx`


 
.. c:var:: double eone_bis


    File : :file:`tot_en.irp.f`

    .. code:: fortran

        real*8	:: etwo	
        real*8	:: eone	
        real*8	:: eone_bis	
        real*8	:: etwo_bis	
        real*8	:: etwo_ter	
        real*8	:: ecore	
        real*8	:: ecore_bis	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci`
       * :c:data:`d0tu`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`nuclear_repulsion`
       * :c:data:`p0tuvx`


 
.. c:var:: double etwo


    File : :file:`tot_en.irp.f`

    .. code:: fortran

        real*8	:: etwo	
        real*8	:: eone	
        real*8	:: eone_bis	
        real*8	:: etwo_bis	
        real*8	:: etwo_ter	
        real*8	:: ecore	
        real*8	:: ecore_bis	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci`
       * :c:data:`d0tu`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`nuclear_repulsion`
       * :c:data:`p0tuvx`


 
.. c:var:: double etwo_bis


    File : :file:`tot_en.irp.f`

    .. code:: fortran

        real*8	:: etwo	
        real*8	:: eone	
        real*8	:: eone_bis	
        real*8	:: etwo_bis	
        real*8	:: etwo_ter	
        real*8	:: ecore	
        real*8	:: ecore_bis	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci`
       * :c:data:`d0tu`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`nuclear_repulsion`
       * :c:data:`p0tuvx`


 
.. c:var:: double etwo_ter


    File : :file:`tot_en.irp.f`

    .. code:: fortran

        real*8	:: etwo	
        real*8	:: eone	
        real*8	:: eone_bis	
        real*8	:: etwo_bis	
        real*8	:: etwo_ter	
        real*8	:: ecore	
        real*8	:: ecore_bis	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci`
       * :c:data:`d0tu`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`nuclear_repulsion`
       * :c:data:`p0tuvx`


 
.. c:var:: int excit


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer, allocatable	:: excit	(2,nMonoEx)
        character*3, allocatable	:: excit_class	(nMonoEx)
        integer, allocatable	:: list_idx_c_a	(3,n_c_a_prov)
        integer, allocatable	:: list_idx_c_v	(3,n_c_v_prov)
        integer, allocatable	:: list_idx_a_v	(3,n_a_v_prov)
        integer, allocatable	:: mat_idx_c_a	(n_core_inact_orb,n_act_orb)
        integer, allocatable	:: mat_idx_c_v	(n_core_inact_orb,n_virt_orb)
        integer, allocatable	:: mat_idx_a_v	(n_act_orb,n_virt_orb)


    a list of the orbitals involved in the excitation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_c_a_prov`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec_old`
       * :c:data:`hessmat`
       * :c:data:`hessmat_old`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`orb_swap`

 
.. c:var:: char excit_class


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer, allocatable	:: excit	(2,nMonoEx)
        character*3, allocatable	:: excit_class	(nMonoEx)
        integer, allocatable	:: list_idx_c_a	(3,n_c_a_prov)
        integer, allocatable	:: list_idx_c_v	(3,n_c_v_prov)
        integer, allocatable	:: list_idx_a_v	(3,n_a_v_prov)
        integer, allocatable	:: mat_idx_c_a	(n_core_inact_orb,n_act_orb)
        integer, allocatable	:: mat_idx_c_v	(n_core_inact_orb,n_virt_orb)
        integer, allocatable	:: mat_idx_a_v	(n_act_orb,n_virt_orb)


    a list of the orbitals involved in the excitation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_c_a_prov`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec_old`
       * :c:data:`hessmat`
       * :c:data:`hessmat_old`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`orb_swap`

 
.. c:var:: double fapq


    File : :file:`mcscf_fock.irp.f`

    .. code:: fortran

        real*8, allocatable	:: fapq	(mo_num,mo_num)


    the active active Fock matrix, in molecular orbitals
    we create them in MOs, quite expensive
    
    for an implementation in AOs we need first the natural orbitals
    for forming an active density matrix in AOs
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`occnum`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`

 
.. c:var:: double fipq


    File : :file:`mcscf_fock.irp.f`

    .. code:: fortran

        real*8, allocatable	:: fipq	(mo_num,mo_num)


    the inactive Fock matrix, in molecular orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`list_core_inact`
       * :c:data:`mo_num`
       * :c:data:`n_core_inact_orb`
       * :c:data:`one_ints_no`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fapq`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`

 
.. c:var:: double gradvec2


    File : :file:`gradient.irp.f`

    .. code:: fortran

        real*8, allocatable	:: gradvec2	(nMonoEx)
        real*8	:: norm_grad_vec2	
        real*8, allocatable	:: norm_grad_vec2_tab	(3)


    calculate the orbital gradient <Psi\| H E_pq \|Psi> from density
    matrices and integrals; Siegbahn et al, Phys Scr 1980
    eqs 14 a,b,c

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`sxmatrix`

 
.. c:var:: double gradvec_old


    File : :file:`grad_old.irp.f`

    .. code:: fortran

        real*8, allocatable	:: gradvec_old	(nMonoEx)


    calculate the orbital gradient <Psi\| H E_pq \|Psi> by hand, i.e. for
    each determinant I we determine the string E_pq \|I> (alpha and beta
    separately) and generate <Psi\|H E_pq \|I>
    sum_I c_I <Psi\|H E_pq \|I> is then the pq component of the orbital
    gradient
    E_pq = a^+_pa_q + a^+_Pa_Q

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`excit`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nmonoex`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`


 
.. c:var:: double hessdiag


    File : :file:`hessian.irp.f`

    .. code:: fortran

        real*8, allocatable	:: hessdiag	(nMonoEx)


    the diagonal of the Hessian, needed for the Davidson procedure

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci_no`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`sxmatrix`

 
.. c:var:: double hessmat


    File : :file:`hessian.irp.f`

    .. code:: fortran

        double precision, allocatable	:: hessmat	(nMonoEx,nMonoEx)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`excit`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`hess_cv_cv`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_c_a_prov`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`orb_swap`
       * :c:data:`sxmatrix`

 
.. c:var:: double hessmat_old


    File : :file:`hessian_old.irp.f`

    .. code:: fortran

        real*8, allocatable	:: hessmat_old	(nMonoEx,nMonoEx)


    calculate the orbital hessian 2 <Psi\| E_pq H E_rs \|Psi>
    + <Psi\| E_pq E_rs H \|Psi> + <Psi\| E_rs E_pq H \|Psi> by hand,
    determinant per determinant, as for the gradient
    
    we assume that we have natural active orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`excit`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nmonoex`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`


 
.. c:var:: double hessmat_peter


    File : :file:`hessian_old.irp.f`

    .. code:: fortran

        real*8, allocatable	:: hessmat_peter	(nMonoEx,nMonoEx)


    explicit hessian matrix from density matrices and integrals
    of course, this will be used for a direct Davidson procedure later
    we will not store the matrix in real life
    formulas are broken down as functions for the 6 classes of matrix elements
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`bavard`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`


 
.. c:var:: int index_orb_swap


    File : :file:`swap_orb.irp.f`

    .. code:: fortran

        integer, allocatable	:: orb_swap	(2,dim_n_max_overlap)
        integer, allocatable	:: index_orb_swap	(dim_n_max_overlap)
        integer	:: n_orb_swap	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`excit`
       * :c:data:`hessmat`
       * :c:data:`max_overlap`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`neworbs`
       * :c:data:`switch_mo_coef`

 
.. c:var:: int list_idx_a_v


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer, allocatable	:: excit	(2,nMonoEx)
        character*3, allocatable	:: excit_class	(nMonoEx)
        integer, allocatable	:: list_idx_c_a	(3,n_c_a_prov)
        integer, allocatable	:: list_idx_c_v	(3,n_c_v_prov)
        integer, allocatable	:: list_idx_a_v	(3,n_a_v_prov)
        integer, allocatable	:: mat_idx_c_a	(n_core_inact_orb,n_act_orb)
        integer, allocatable	:: mat_idx_c_v	(n_core_inact_orb,n_virt_orb)
        integer, allocatable	:: mat_idx_a_v	(n_act_orb,n_virt_orb)


    a list of the orbitals involved in the excitation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_c_a_prov`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec_old`
       * :c:data:`hessmat`
       * :c:data:`hessmat_old`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`orb_swap`

 
.. c:var:: int list_idx_c_a


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer, allocatable	:: excit	(2,nMonoEx)
        character*3, allocatable	:: excit_class	(nMonoEx)
        integer, allocatable	:: list_idx_c_a	(3,n_c_a_prov)
        integer, allocatable	:: list_idx_c_v	(3,n_c_v_prov)
        integer, allocatable	:: list_idx_a_v	(3,n_a_v_prov)
        integer, allocatable	:: mat_idx_c_a	(n_core_inact_orb,n_act_orb)
        integer, allocatable	:: mat_idx_c_v	(n_core_inact_orb,n_virt_orb)
        integer, allocatable	:: mat_idx_a_v	(n_act_orb,n_virt_orb)


    a list of the orbitals involved in the excitation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_c_a_prov`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec_old`
       * :c:data:`hessmat`
       * :c:data:`hessmat_old`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`orb_swap`

 
.. c:var:: int list_idx_c_v


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer, allocatable	:: excit	(2,nMonoEx)
        character*3, allocatable	:: excit_class	(nMonoEx)
        integer, allocatable	:: list_idx_c_a	(3,n_c_a_prov)
        integer, allocatable	:: list_idx_c_v	(3,n_c_v_prov)
        integer, allocatable	:: list_idx_a_v	(3,n_a_v_prov)
        integer, allocatable	:: mat_idx_c_a	(n_core_inact_orb,n_act_orb)
        integer, allocatable	:: mat_idx_c_v	(n_core_inact_orb,n_virt_orb)
        integer, allocatable	:: mat_idx_a_v	(n_act_orb,n_virt_orb)


    a list of the orbitals involved in the excitation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_c_a_prov`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec_old`
       * :c:data:`hessmat`
       * :c:data:`hessmat_old`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`orb_swap`

 
.. c:var:: double lowest_super_ci_coef_mo


    File : :file:`superci_dm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: lowest_super_ci_coef_mo	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`excit`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`natorbsci`
       * :c:data:`nmonoex`
       * :c:data:`sxeigenvec`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`super_ci_dm`

 
.. c:var:: int mat_idx_a_v


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer, allocatable	:: excit	(2,nMonoEx)
        character*3, allocatable	:: excit_class	(nMonoEx)
        integer, allocatable	:: list_idx_c_a	(3,n_c_a_prov)
        integer, allocatable	:: list_idx_c_v	(3,n_c_v_prov)
        integer, allocatable	:: list_idx_a_v	(3,n_a_v_prov)
        integer, allocatable	:: mat_idx_c_a	(n_core_inact_orb,n_act_orb)
        integer, allocatable	:: mat_idx_c_v	(n_core_inact_orb,n_virt_orb)
        integer, allocatable	:: mat_idx_a_v	(n_act_orb,n_virt_orb)


    a list of the orbitals involved in the excitation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_c_a_prov`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec_old`
       * :c:data:`hessmat`
       * :c:data:`hessmat_old`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`orb_swap`

 
.. c:var:: int mat_idx_c_a


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer, allocatable	:: excit	(2,nMonoEx)
        character*3, allocatable	:: excit_class	(nMonoEx)
        integer, allocatable	:: list_idx_c_a	(3,n_c_a_prov)
        integer, allocatable	:: list_idx_c_v	(3,n_c_v_prov)
        integer, allocatable	:: list_idx_a_v	(3,n_a_v_prov)
        integer, allocatable	:: mat_idx_c_a	(n_core_inact_orb,n_act_orb)
        integer, allocatable	:: mat_idx_c_v	(n_core_inact_orb,n_virt_orb)
        integer, allocatable	:: mat_idx_a_v	(n_act_orb,n_virt_orb)


    a list of the orbitals involved in the excitation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_c_a_prov`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec_old`
       * :c:data:`hessmat`
       * :c:data:`hessmat_old`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`orb_swap`

 
.. c:var:: int mat_idx_c_v


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer, allocatable	:: excit	(2,nMonoEx)
        character*3, allocatable	:: excit_class	(nMonoEx)
        integer, allocatable	:: list_idx_c_a	(3,n_c_a_prov)
        integer, allocatable	:: list_idx_c_v	(3,n_c_v_prov)
        integer, allocatable	:: list_idx_a_v	(3,n_a_v_prov)
        integer, allocatable	:: mat_idx_c_a	(n_core_inact_orb,n_act_orb)
        integer, allocatable	:: mat_idx_c_v	(n_core_inact_orb,n_virt_orb)
        integer, allocatable	:: mat_idx_a_v	(n_act_orb,n_virt_orb)


    a list of the orbitals involved in the excitation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_c_a_prov`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec_old`
       * :c:data:`hessmat`
       * :c:data:`hessmat_old`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`orb_swap`

 
.. c:var:: double mat_tmp_dm_super_ci


    File : :file:`superci_dm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mat_tmp_dm_super_ci	(n_act_orb,n_act_orb)


    computation of the term in [ ] in the equation B3.d of Roos et. al. Chemical Physics 48 (1980) 157-173
    
    !!!!! WARNING !!!!!! there is a TYPO: a MINUS SIGN SHOULD APPEAR in that term

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`super_ci_dm`

 
.. c:var:: int max_overlap


    File : :file:`swap_orb.irp.f`

    .. code:: fortran

        integer, allocatable	:: max_overlap	(nMonoEx)
        integer	:: n_max_overlap	
        integer	:: dim_n_max_overlap	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nmonoex`
       * :c:data:`sxvector_lowest`
       * :c:data:`thresh_overlap_switch`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`orb_swap`

 
.. c:var:: double mcscf_fock_alpha_ao


    File : :file:`mcscf_fock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mcscf_fock_alpha_ao	(ao_num,ao_num)
        double precision, allocatable	:: mcscf_fock_beta_ao	(ao_num,ao_num)


    mcscf_fock_alpha_ao are set to usual Fock like operator but computed with the MCSCF densities on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`d0tu_alpha_ao`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mcscf_fock_alpha_mo`

 
.. c:var:: double mcscf_fock_alpha_mo


    File : :file:`mcscf_fock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mcscf_fock_alpha_mo	(mo_num,mo_num)
        double precision, allocatable	:: mcscf_fock_beta_mo	(mo_num,mo_num)


    Mo_mcscf_fock_alpha are set to usual Fock like operator but computed with the MCSCF densities on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mcscf_fock_alpha_ao`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mcscf_fock_mo`

 
.. c:var:: double mcscf_fock_beta_ao


    File : :file:`mcscf_fock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mcscf_fock_alpha_ao	(ao_num,ao_num)
        double precision, allocatable	:: mcscf_fock_beta_ao	(ao_num,ao_num)


    mcscf_fock_alpha_ao are set to usual Fock like operator but computed with the MCSCF densities on the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`d0tu_alpha_ao`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mcscf_fock_alpha_mo`

 
.. c:var:: double mcscf_fock_beta_mo


    File : :file:`mcscf_fock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mcscf_fock_alpha_mo	(mo_num,mo_num)
        double precision, allocatable	:: mcscf_fock_beta_mo	(mo_num,mo_num)


    Mo_mcscf_fock_alpha are set to usual Fock like operator but computed with the MCSCF densities on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mcscf_fock_alpha_ao`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mcscf_fock_mo`

 
.. c:var:: double mcscf_fock_diag_mo


    File : :file:`mcscf_fock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mcscf_fock_mo	(mo_num,mo_num)
        double precision, allocatable	:: mcscf_fock_diag_mo	(mo_num)


    MCSF Fock matrix on the MO basis.
    For open shells, the ROHF Fock Matrix is ::
    
          |  Rcc  |  F^b  |  Fcv  |
          |-----------------------|
          |  F^b  |  Roo  |  F^a  |
          |-----------------------|
          |  Fcv  |  F^a  |  Rvv  |
    
    C: Core, O: Open, V: Virtual
    
    Rcc = Acc Fcc^a + Bcc Fcc^b
    Roo = Aoo Foo^a + Boo Foo^b
    Rvv = Avv Fvv^a + Bvv Fvv^b
    Fcv = (F^a + F^b)/2
    
    F^a: Fock matrix alpha (MO), F^b: Fock matrix beta (MO)
    A,B: Coupling parameters
    
    J. Chem. Phys. 133, 141102 (2010), https://doi.org/10.1063/1.3503173
    Coupling parameters from J. Chem. Phys. 125, 204110 (2006); https://doi.org/10.1063/1.2393223.
    
    ::
    
     cc   oo   vv
     A  -0.5  0.5  1.5
     B   1.5  0.5 -0.5
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mcscf_fock_alpha_mo`
       * :c:data:`mo_num`


 
.. c:var:: double mcscf_fock_mo


    File : :file:`mcscf_fock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mcscf_fock_mo	(mo_num,mo_num)
        double precision, allocatable	:: mcscf_fock_diag_mo	(mo_num)


    MCSF Fock matrix on the MO basis.
    For open shells, the ROHF Fock Matrix is ::
    
          |  Rcc  |  F^b  |  Fcv  |
          |-----------------------|
          |  F^b  |  Roo  |  F^a  |
          |-----------------------|
          |  Fcv  |  F^a  |  Rvv  |
    
    C: Core, O: Open, V: Virtual
    
    Rcc = Acc Fcc^a + Bcc Fcc^b
    Roo = Aoo Foo^a + Boo Foo^b
    Rvv = Avv Fvv^a + Bvv Fvv^b
    Fcv = (F^a + F^b)/2
    
    F^a: Fock matrix alpha (MO), F^b: Fock matrix beta (MO)
    A,B: Coupling parameters
    
    J. Chem. Phys. 133, 141102 (2010), https://doi.org/10.1063/1.3503173
    Coupling parameters from J. Chem. Phys. 125, 204110 (2006); https://doi.org/10.1063/1.2393223.
    
    ::
    
     cc   oo   vv
     A  -0.5  0.5  1.5
     B   1.5  0.5 -0.5
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mcscf_fock_alpha_mo`
       * :c:data:`mo_num`


 
.. c:var:: int n_a_v_prov


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer	:: n_c_a_prov	
        integer	:: n_c_v_prov	
        integer	:: n_a_v_prov	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`excit`
       * :c:data:`hessmat`

 
.. c:var:: int n_c_a_prov


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer	:: n_c_a_prov	
        integer	:: n_c_v_prov	
        integer	:: n_a_v_prov	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`excit`
       * :c:data:`hessmat`

 
.. c:var:: int n_c_v_prov


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer	:: n_c_a_prov	
        integer	:: n_c_v_prov	
        integer	:: n_a_v_prov	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`excit`
       * :c:data:`hessmat`

 
.. c:var:: int n_guess_sx_mat


    File : :file:`neworbs.irp.f`

    .. code:: fortran

        real*8, allocatable	:: sxmatrix	(nMonoEx+1,nMonoEx+1)
        integer	:: n_guess_sx_mat	


    Single-excitation matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`diag_hess_cas`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`level_shift_casscf`
       * :c:data:`nmonoex`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`sxeigenvec`

 
.. c:var:: int n_max_overlap


    File : :file:`swap_orb.irp.f`

    .. code:: fortran

        integer, allocatable	:: max_overlap	(nMonoEx)
        integer	:: n_max_overlap	
        integer	:: dim_n_max_overlap	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nmonoex`
       * :c:data:`sxvector_lowest`
       * :c:data:`thresh_overlap_switch`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`orb_swap`

 
.. c:var:: int n_orb_swap


    File : :file:`swap_orb.irp.f`

    .. code:: fortran

        integer, allocatable	:: orb_swap	(2,dim_n_max_overlap)
        integer, allocatable	:: index_orb_swap	(dim_n_max_overlap)
        integer	:: n_orb_swap	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`excit`
       * :c:data:`hessmat`
       * :c:data:`max_overlap`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`neworbs`
       * :c:data:`switch_mo_coef`

 
.. c:var:: double natorbsci


    File : :file:`natorb.irp.f`

    .. code:: fortran

        real*8, allocatable	:: natorbsci	(n_act_orb,n_act_orb)
        real*8, allocatable	:: occ_act	(n_act_orb)


    Natural orbitals of CI

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`d0tu`
       * :c:data:`n_act_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_no_array`
       * :c:data:`bielec_pxxq_no_array`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`natorbsci_mos`
       * :c:data:`occnum`
       * :c:data:`one_ints_no`
       * :c:data:`p0tuvx_no`
       * :c:data:`super_ci_dm`

 
.. c:var:: double natorbsci_mos


    File : :file:`natorb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: natorbsci_mos	(mo_num,mo_num)


    Rotation matrix from current MOs to the CI natural MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`natorbsci`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`natorbsfci`

 
.. c:var:: double natorbsfci


    File : :file:`natorb.irp.f`

    .. code:: fortran

        real*8, allocatable	:: natorbsfci	(ao_num,mo_num)


    FCI natural orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`natorbsci_mos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`neworbs`
       * :c:data:`superci_natorb`
       * :c:data:`switch_mo_coef`

 
.. c:var:: double neworbs


    File : :file:`neworbs.irp.f`

    .. code:: fortran

        double precision, allocatable	:: neworbs	(ao_num,mo_num)


    Updated orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`best_vector_ovrlp_casscf`
       * :c:data:`level_shift_casscf`
       * :c:data:`mo_num`
       * :c:data:`natorbsfci`
       * :c:data:`orb_swap`
       * :c:data:`state_following_casscf`
       * :c:data:`switch_mo_coef`
       * :c:data:`umat`


 
.. c:var:: int nmonoex


    File : :file:`gradient.irp.f`

    .. code:: fortran

        integer	:: nmonoex	


    Number of single excitations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`best_vector_ovrlp_casscf`
       * :c:data:`excit`
       * :c:data:`gradvec2`
       * :c:data:`gradvec_old`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_old`
       * :c:data:`hessmat_peter`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`max_overlap`
       * :c:data:`sxeigenvec`
       * :c:data:`sxmatrix`
       * :c:data:`sxvector`
       * :c:data:`sxvector_lowest`

 
.. c:var:: double norm_grad_vec2


    File : :file:`gradient.irp.f`

    .. code:: fortran

        real*8, allocatable	:: gradvec2	(nMonoEx)
        real*8	:: norm_grad_vec2	
        real*8, allocatable	:: norm_grad_vec2_tab	(3)


    calculate the orbital gradient <Psi\| H E_pq \|Psi> from density
    matrices and integrals; Siegbahn et al, Phys Scr 1980
    eqs 14 a,b,c

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`sxmatrix`

 
.. c:var:: double norm_grad_vec2_tab


    File : :file:`gradient.irp.f`

    .. code:: fortran

        real*8, allocatable	:: gradvec2	(nMonoEx)
        real*8	:: norm_grad_vec2	
        real*8, allocatable	:: norm_grad_vec2_tab	(3)


    calculate the orbital gradient <Psi\| H E_pq \|Psi> from density
    matrices and integrals; Siegbahn et al, Phys Scr 1980
    eqs 14 a,b,c

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`nmonoex`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`sxmatrix`

 
.. c:var:: double occ_act


    File : :file:`natorb.irp.f`

    .. code:: fortran

        real*8, allocatable	:: natorbsci	(n_act_orb,n_act_orb)
        real*8, allocatable	:: occ_act	(n_act_orb)


    Natural orbitals of CI

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`d0tu`
       * :c:data:`n_act_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_no_array`
       * :c:data:`bielec_pxxq_no_array`
       * :c:data:`bielecci_no`
       * :c:data:`cholesky_no_1_idx_transp`
       * :c:data:`cholesky_no_2_idx_transp`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`natorbsci_mos`
       * :c:data:`occnum`
       * :c:data:`one_ints_no`
       * :c:data:`p0tuvx_no`
       * :c:data:`super_ci_dm`

 
.. c:var:: double occnum


    File : :file:`natorb.irp.f`

    .. code:: fortran

        real*8, allocatable	:: occnum	(mo_num)


    MO occupation numbers

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`natorbsci`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fapq`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`

 
.. c:var:: double one_ints_no


    File : :file:`natorb.irp.f`

    .. code:: fortran

        real*8, allocatable	:: one_ints_no	(mo_num,mo_num)


    Transformed one-e integrals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`n_act_orb`
       * :c:data:`natorbsci`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fipq`

 
.. c:var:: int orb_swap


    File : :file:`swap_orb.irp.f`

    .. code:: fortran

        integer, allocatable	:: orb_swap	(2,dim_n_max_overlap)
        integer, allocatable	:: index_orb_swap	(dim_n_max_overlap)
        integer	:: n_orb_swap	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`excit`
       * :c:data:`hessmat`
       * :c:data:`max_overlap`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`neworbs`
       * :c:data:`switch_mo_coef`

 
.. c:var:: double p0tuvx


    File : :file:`densities.irp.f`

    .. code:: fortran

        real*8, allocatable	:: p0tuvx	(n_act_orb,n_act_orb,n_act_orb,n_act_orb)


    The second-order density matrix in the basis of the starting MOs ONLY IN THE RANGE OF ACTIVE MOS
    The values are state averaged
    
    We use the spin-free generators of mono-excitations
    E_pq destroys q and creates p
    D_pq   =     <0|E_pq|0> = D_qp
    P_pqrs = 1/2 <0|E_pq E_rs - delta_qr E_ps|0>
    
    P0tuvx(p,q,r,s) = chemist notation : 1/2 <0|E_pq E_rs - delta_qr E_ps|0>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`fast_2rdm`
       * :c:data:`n_act_orb`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`p0tuvx_peter`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`etwo`
       * :c:data:`p0tuvx_no`

 
.. c:var:: double p0tuvx_no


    File : :file:`natorb.irp.f`

    .. code:: fortran

        real*8, allocatable	:: p0tuvx_no	(n_act_orb,n_act_orb,n_act_orb,n_act_orb)


    4-index transformation of 2part matrices

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_act_orb`
       * :c:data:`natorbsci`
       * :c:data:`p0tuvx`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`hessmat_peter`
       * :c:data:`super_ci_dm`

 
.. c:var:: double p0tuvx_peter


    File : :file:`densities_peter.irp.f`

    .. code:: fortran

        real*8, allocatable	:: p0tuvx_peter	(n_act_orb,n_act_orb,n_act_orb,n_act_orb)


    the second-order density matrix in the basis of the starting MOs
    matrices are state averaged
    
    we use the spin-free generators of mono-excitations
    E_pq destroys q and creates p
    D_pq   =     <0|E_pq|0> = D_qp
    P_pqrs = 1/2 <0|E_pq E_rs - delta_qr E_ps|0>
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`n_act_orb`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`p0tuvx`

 
.. c:var:: double super_ci_dm


    File : :file:`superci_dm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: super_ci_dm	(mo_num,mo_num)


    density matrix of the super CI matrix, in the basis of NATURAL ORBITALS OF THE CASCI WF
    
    This is obtained from annex B of Roos et. al. Chemical Physics 48 (1980) 157-173
    
    WARNING ::: in the equation B3.d there is a TYPO with a forgotten MINUS SIGN (see variable mat_tmp_dm_super_ci )

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`mat_tmp_dm_super_ci`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`natorbsci`
       * :c:data:`p0tuvx_no`
       * :c:data:`sxeigenvec`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`superci_natorb`

 
.. c:var:: double superci_nat_occ


    File : :file:`superci_dm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: superci_natorb	(ao_num,mo_num)
        double precision, allocatable	:: superci_nat_occ	(mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_num`
       * :c:data:`natorbsfci`
       * :c:data:`super_ci_dm`


 
.. c:var:: double superci_natorb


    File : :file:`superci_dm.irp.f`

    .. code:: fortran

        double precision, allocatable	:: superci_natorb	(ao_num,mo_num)
        double precision, allocatable	:: superci_nat_occ	(mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_num`
       * :c:data:`natorbsfci`
       * :c:data:`super_ci_dm`


 
.. c:var:: double switch_mo_coef


    File : :file:`swap_orb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: switch_mo_coef	(ao_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_num`
       * :c:data:`natorbsfci`
       * :c:data:`orb_swap`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`neworbs`

 
.. c:var:: double sxeigenval


    File : :file:`neworbs.irp.f`

    .. code:: fortran

        real*8, allocatable	:: sxeigenvec	(nMonoEx+1,nMonoEx+1)
        real*8, allocatable	:: sxeigenval	(nMonoEx+1)


    Eigenvectors/eigenvalues of the single-excitation matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`n_det_max_full`
       * :c:data:`n_states_diag`
       * :c:data:`nmonoex`
       * :c:data:`sxmatrix`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`best_vector_ovrlp_casscf`
       * :c:data:`energy_improvement`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`super_ci_dm`
       * :c:data:`sxvector`
       * :c:data:`sxvector_lowest`

 
.. c:var:: double sxeigenvec


    File : :file:`neworbs.irp.f`

    .. code:: fortran

        real*8, allocatable	:: sxeigenvec	(nMonoEx+1,nMonoEx+1)
        real*8, allocatable	:: sxeigenval	(nMonoEx+1)


    Eigenvectors/eigenvalues of the single-excitation matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`n_det_max_full`
       * :c:data:`n_states_diag`
       * :c:data:`nmonoex`
       * :c:data:`sxmatrix`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`best_vector_ovrlp_casscf`
       * :c:data:`energy_improvement`
       * :c:data:`lowest_super_ci_coef_mo`
       * :c:data:`super_ci_dm`
       * :c:data:`sxvector`
       * :c:data:`sxvector_lowest`

 
.. c:var:: double sxmatrix


    File : :file:`neworbs.irp.f`

    .. code:: fortran

        real*8, allocatable	:: sxmatrix	(nMonoEx+1,nMonoEx+1)
        integer	:: n_guess_sx_mat	


    Single-excitation matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`diag_hess_cas`
       * :c:data:`gradvec2`
       * :c:data:`hessdiag`
       * :c:data:`hessmat`
       * :c:data:`level_shift_casscf`
       * :c:data:`nmonoex`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`sxeigenvec`

 
.. c:var:: double sxvector


    File : :file:`neworbs.irp.f`

    .. code:: fortran

        double precision, allocatable	:: sxvector	(nMonoEx+1)


    Best eigenvector of the single-excitation matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`best_vector_ovrlp_casscf`
       * :c:data:`nmonoex`
       * :c:data:`sxeigenvec`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`umat`

 
.. c:var:: double sxvector_lowest


    File : :file:`swap_orb.irp.f`

    .. code:: fortran

        double precision, allocatable	:: sxvector_lowest	(nMonoEx)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nmonoex`
       * :c:data:`sxeigenvec`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`max_overlap`

 
.. c:var:: double thresh_overlap_switch


    File : :file:`swap_orb.irp.f`

    .. code:: fortran

        double precision	:: thresh_overlap_switch	



    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`max_overlap`

 
.. c:var:: double umat


    File : :file:`neworbs.irp.f`

    .. code:: fortran

        real*8, allocatable	:: umat	(mo_num,mo_num)


    Orbital rotation matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bavard`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`sxvector`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`neworbs`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void bielec_no_basis(void)


    File : :file:`chol_bielec.irp.f`

    .. code:: fortran

        double precision function bielec_no_basis(i_1,j_1,i_2,j_2)


    integral (i_1 j_1|i_2 j_2) in the mixed basis of both MOs and natural MOs
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`

 
.. c:function:: void bielec_pqxx(void)


    File : :file:`chol_bielec.irp.f`

    .. code:: fortran

        double precision function bielec_PQxx(i_mo, j_mo, i_ca, j_ca)


    function that computes (i_mo j_mo \|i_ca j_ca) with Cholesky decomposition
    
    indices are unshifted orbital numbers
    
    where i_ca, j_ca are in [1:n_core_inact_act_orb]

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core_inact_act`

 
.. c:function:: void bielec_pqxx_no(void)


    File : :file:`chol_bielec.irp.f`

    .. code:: fortran

        double precision function bielec_PQxx_no(i_mo, j_mo, i_ca, j_ca)


    function that computes (i_mo j_mo| i_ca j_ca) with Cholesky decomposition  on the NO basis for active orbitals
    
    where i_ca, j_ca are in [1:n_core_inact_act_orb]

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`list_core_inact_act`

 
.. c:function:: void bielec_pxxq(void)


    File : :file:`chol_bielec.irp.f`

    .. code:: fortran

        double precision function bielec_PxxQ(i_mo, i_ca, j_ca, j_mo)


    function that computes (i_mo j_mo \|i_ca j_ca) with Cholesky decomposition
    
    where i_ca, j_ca are in [1:n_core_inact_act_orb]

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_core_inact_act`

 
.. c:function:: void bielec_pxxq_no(void)


    File : :file:`chol_bielec.irp.f`

    .. code:: fortran

        double precision function bielec_PxxQ_no(i_mo, j_ca, i_ca, j_mo)


    function that computes (i_mo j_ca \|i_ca j_mo) with Cholesky decomposition  on the NO basis for active orbitals
    
    where i_ca, j_ca are in [1:n_core_inact_act_orb]

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`list_core_inact_act`

 
.. c:function:: void calc_grad_elem(void)


    File : :file:`grad_old.irp.f`

    .. code:: fortran

        subroutine calc_grad_elem(ihole,ipart,res)


    eq 18 of Siegbahn et al, Physica Scripta 1980
    we calculate 2 <Psi\| H E_pq \| Psi>, q=hole, p=particle

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`gradvec_old`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`det_copy`
       * :c:func:`det_extract`
       * :c:func:`do_signed_mono_excitation`
       * :c:func:`i_h_psi`

 
.. c:function:: void calc_hess_elem(void)


    File : :file:`hessian_old.irp.f`

    .. code:: fortran

        subroutine calc_hess_elem(ihole,ipart,jhole,jpart,res)


    eq 19 of Siegbahn et al, Physica Scripta 1980
    we calculate  2 <Psi\| E_pq H E_rs \|Psi>
    + <Psi\| E_pq E_rs H \|Psi> + <Psi\| E_rs E_pq H \|Psi>
    average over all states is performed.
    no transition between states.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`hessmat_old`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`det_copy`
       * :c:func:`det_extract`
       * :c:func:`do_signed_mono_excitation`
       * :c:func:`i_h_j`
       * :c:func:`i_h_psi`

 
.. c:function:: void davidson_diag_sx_mat(void)


    File : :file:`dav_sx_mat.irp.f`

    .. code:: fortran

        subroutine davidson_diag_sx_mat(N_st, u_in, energies)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states_diag`
       * :c:data:`nmonoex`
       * :c:data:`sxmatrix`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`sxeigenvec`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`davidson_general`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`nthreads_davidson`

 
.. c:function:: void det_copy(void)


    File : :file:`det_manip.irp.f`

    .. code:: fortran

        subroutine det_copy(key1,key2,Nint)


    copy a determinant from key1 to key2

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`calc_grad_elem`
       * :c:func:`calc_hess_elem`
       * :c:func:`do_signed_mono_excitation`
       * :c:data:`p0tuvx_peter`

 
.. c:function:: void det_extract(void)


    File : :file:`det_manip.irp.f`

    .. code:: fortran

        subroutine det_extract(key,nu,Nint)


    extract a determinant from the list of determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`calc_grad_elem`
       * :c:func:`calc_hess_elem`
       * :c:func:`do_signed_mono_excitation`
       * :c:data:`p0tuvx_peter`

 
.. c:function:: void do_signed_mono_excitation(void)


    File : :file:`det_manip.irp.f`

    .. code:: fortran

        subroutine do_signed_mono_excitation(key1,key2,nu,ihole,ipart,       &
              ispin,phase,ierr)


    we create the mono-excitation, and determine, if possible,
    the phase and the number in the list of determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`calc_grad_elem`
       * :c:func:`calc_hess_elem`
       * :c:func:`do_spinfree_mono_excitation`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`det_copy`
       * :c:func:`det_extract`
       * :c:func:`do_single_excitation`
       * :c:func:`get_single_excitation`

 
.. c:function:: void do_spinfree_mono_excitation(void)


    File : :file:`det_manip.irp.f`

    .. code:: fortran

        subroutine do_spinfree_mono_excitation(key_in,key_out1,key_out2      &
              ,nu1,nu2,ihole,ipart,phase1,phase2,ierr,jerr)


    we create the spin-free mono-excitation E_pq=(a^+_p a_q + a^+_P a_Q)
    we may create two determinants as result
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`p0tuvx_peter`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`do_signed_mono_excitation`

 
.. c:function:: void driver_optorb(void)


    File : :file:`driver_optorb.irp.f`

    .. code:: fortran

        subroutine driver_optorb



 
.. c:function:: void gradvec_ia(void)


    File : :file:`gradient.irp.f`

    .. code:: fortran

        real*8 function gradvec_ia(i,a)


    the orbital gradient core/inactive -> virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`

 
.. c:function:: void gradvec_it(void)


    File : :file:`gradient.irp.f`

    .. code:: fortran

        real*8 function gradvec_it(i,t)


    the orbital gradient core/inactive -> active
    we assume natural orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`n_act_orb`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

 
.. c:function:: void gradvec_ta(void)


    File : :file:`gradient.irp.f`

    .. code:: fortran

        real*8 function gradvec_ta(t,a)


    the orbital gradient active -> virtual
    we assume natural orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci_no`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

 
.. c:function:: void hessmat_iajb(void)


    File : :file:`hessian.irp.f`

    .. code:: fortran

        real*8 function hessmat_iajb(i,a,j,b)


    the orbital hessian for core/inactive -> virtual, core/inactive -> virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`

 
.. c:function:: void hessmat_iatb(void)


    File : :file:`hessian.irp.f`

    .. code:: fortran

        real*8 function hessmat_iatb(i,a,t,b)


    the orbital hessian for core/inactive -> virtual, active -> virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci_no`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

 
.. c:function:: void hessmat_itja(void)


    File : :file:`hessian.irp.f`

    .. code:: fortran

        real*8 function hessmat_itja(i,t,j,a)


    the orbital hessian for core/inactive -> active, core/inactive -> virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci_no`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

 
.. c:function:: void hessmat_itju(void)


    File : :file:`hessian.irp.f`

    .. code:: fortran

        real*8 function hessmat_itju(i,t,j,u)


    the orbital hessian for core/inactive -> active, core/inactive -> active
    i, t, j, u are list indices, the corresponding orbitals are ii,tt,jj,uu
    
    we assume natural orbitals

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci_no`
       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`n_act_orb`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

 
.. c:function:: void hessmat_itua(void)


    File : :file:`hessian.irp.f`

    .. code:: fortran

        real*8 function hessmat_itua(i,t,u,a)


    the orbital hessian for core/inactive -> active, active -> virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fapq`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_core_inact`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

 
.. c:function:: void hessmat_taub(void)


    File : :file:`hessian.irp.f`

    .. code:: fortran

        real*8 function hessmat_taub(t,a,u,b)


    the orbital hessian for act->virt,act->virt

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielecci_no`
       * :c:data:`cholesky_mo_num`
       * :c:data:`cholesky_no_total_transp`
       * :c:data:`fipq`
       * :c:data:`list_act`
       * :c:data:`list_virt`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`occnum`
       * :c:data:`p0tuvx_no`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dgemv`

 
.. c:function:: void print_grad(void)


    File : :file:`get_energy.irp.f`

    .. code:: fortran

        subroutine print_grad



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`excit`
       * :c:data:`gradvec2`
       * :c:data:`nmonoex`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`get_energy`

 
.. c:function:: void reorder_orbitals_for_casscf(void)


    File : :file:`reorder_orb.irp.f`

    .. code:: fortran

        subroutine reorder_orbitals_for_casscf


    routine that reorders the orbitals of the CASSCF in terms block of core, active and virtual

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`list_act`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_core`
       * :c:data:`list_del`
       * :c:data:`list_inact`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`list_virt`
       * :c:data:`mo_class`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_del_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`casscf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`isort`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`list_virt`
       * :c:data:`mo_coef`

 
.. c:function:: void routine(void)


    File : :file:`get_energy.irp.f`

    .. code:: fortran

        subroutine routine



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_act`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_act_orb`
       * :c:data:`state_av_act_2_rdm_spin_trace_mo`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`test_dav`

 
.. c:function:: void routine_bielec_pqxx_no(void)


    File : :file:`test_chol.irp.f`

    .. code:: fortran

        subroutine routine_bielec_PQxx_no



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_no_array`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_num`
       * :c:data:`n_core_inact_act_orb`

 
.. c:function:: void routine_bielec_pxxq_no_array(void)


    File : :file:`test_chol.irp.f`

    .. code:: fortran

        subroutine routine_bielec_PxxQ_no_array



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pxxq_no_array`
       * :c:data:`list_core_inact_act`
       * :c:data:`mo_num`
       * :c:data:`n_core_inact_act_orb`

 
.. c:function:: void save_energy(void)


    File : :file:`save_energy.irp.f`

    .. code:: fortran

        subroutine save_energy(E,pt2)


    Saves the energy in |EZFIO|.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run_cipsi`
       * :c:func:`run_stochastic_cipsi`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_casscf_cipsi_energy`
       * :c:func:`ezfio_set_casscf_cipsi_energy_pt2`

 
.. c:function:: void test_bielec_pqxx(void)


    File : :file:`test_chol.irp.f`

    .. code:: fortran

        subroutine test_bielec_PQxx(i_mo, j_mo, i_ca, j_ca)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pqxx_array`
       * :c:data:`mo_num`
       * :c:data:`n_core_inact_act_orb`

 
.. c:function:: void test_bielec_pxxq_chol(void)


    File : :file:`test_chol.irp.f`

    .. code:: fortran

        subroutine test_bielec_PxxQ_chol(i_mo, i_ca, j_ca, j_mo)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`bielec_pxxq_array`
       * :c:data:`mo_num`
       * :c:data:`n_core_inact_act_orb`

