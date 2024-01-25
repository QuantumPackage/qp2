===========
bi_ort_ints
===========

This module contains all necessary integrals for the TC Hamiltonian in a bi-orthonormal (BO) MO Basis.
See in bi_ortho_basis for more information. 
The main providers are : 

One-electron integrals 
----------------------
+) ao_one_e_integrals_tc_tot : total one-electron Hamiltonian which might include non hermitian part coming from one-e correlation factor. 
+) mo_bi_ortho_tc_one_e : one-electron Hamiltonian (h_core+one-J terms) on the BO-MO basis. 
+) mo_bi_orth_bipole_x  : x-component of the dipole operator on the BO-MO basis. (Same for y,z) 

Two-electron integrals 
----------------------
+) ao_two_e_tc_tot : Total two-electron operator (including the non-hermitian term of the TC Hamiltonian) on the AO basis
+) mo_bi_ortho_tc_two_e : Total two-electron operator on the BO-MO basis

Three-electron integrals 
------------------------
+) three_body_ints_bi_ort : 6-indices three-electron tensor (-L) on the BO-MO basis. WARNING :: N^6 storage !
+) three_e_3_idx_direct_bi_ort : DIRECT term with 3 different indices of the -L operator. These terms appear in the DIAGONAL matrix element of the -L operator. The 5 other permutations needed to compute matrix elements can be found in three_body_ijm.irp.f 
+) three_e_4_idx_direct_bi_ort : DIRECT term with 4 different indices of the -L operator. These terms appear in the OFF-DIAGONAL matrix element of the -L operator including SINGLE EXCITATIONS. The 5 other permutations needed to compute matrix elements can be found in three_body_ijmk.irp.f 
+) three_e_5_idx_direct_bi_ort : DIRECT term with 5 different indices of the -L operator. These terms appear in the OFF-DIAGONAL matrix element of the -L operator including DOUBLE EXCITATIONS. The 5 other permutations needed to compute matrix elements can be found in three_body_ijmkl.irp.f 
