=============
non_h_ints_mu
=============

Computes the non hermitian potential of the mu-TC Hamiltonian on the AO and BI-ORTHO MO basis.
The operator is defined in Eq. 33 of JCP 154, 084119 (2021)

The two providers are :
+) ao_non_hermit_term_chemist which returns the non hermitian part of the two-electron TC Hamiltonian on the MO basis. 
+) mo_non_hermit_term_chemist which returns the non hermitian part of the two-electron TC Hamiltonian on the BI-ORTHO MO basis. 


!\sum_mm = 1,3 \sum_R phi_i(R) \phi_k(R) grad_1_u_ij_mu(j,l,R,mm)  grad_1_u_ij_mu(m,n,R,mm)
!\sum_mm+= 1,3 \sum_R phi_j(R) \phi_l(R) grad_1_u_ij_mu(i,k,R,mm)  grad_1_u_ij_mu(m,n,R,mm)
!\sum_mm+= 1,3 \sum_R phi_m(R) \phi_n(R) grad_1_u_ij_mu(i,k,R,mm)  grad_1_u_ij_mu(j,l,R,mm)
