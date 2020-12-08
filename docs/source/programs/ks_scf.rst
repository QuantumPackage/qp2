.. _ks_scf: 
 
.. program:: ks_scf 
 
====== 
ks_scf 
====== 
 
 
 
 
 Produce `Kohn_Sham` MO orbital 
 output: mo_basis.mo_num mo_basis.mo_label mo_basis.ao_md5 mo_basis.mo_coef mo_basis.mo_occ 
 output: kohn_sham.energy 
 optional: mo_basis.mo_coef 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`density_for_dft` 
    * :c:data:`io_ao_one_e_integrals` 
    * :c:data:`io_mo_one_e_integrals` 
    * :c:data:`mu_erf_dft` 
 
 Calls: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:func:`check_coherence_functional` 
    * :c:func:`create_guess` 
    * :c:func:`orthonormalize_mos` 
    * :c:func:`run` 
 
 Touches: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`fock_matrix_ao_alpha` 
    * :c:data:`fock_matrix_ao_alpha` 
    * :c:data:`fock_matrix_ao_alpha` 
    * :c:data:`fock_matrix_ao_alpha` 
    * :c:data:`mo_coef` 
    * :c:data:`density_for_dft` 
    * :c:data:`io_ao_one_e_integrals` 
    * :c:data:`io_mo_one_e_integrals` 
    * :c:data:`level_shift` 
    * :c:data:`mo_coef` 
    * :c:data:`mo_label` 
