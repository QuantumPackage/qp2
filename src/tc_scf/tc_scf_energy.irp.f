
 BEGIN_PROVIDER [ double precision, TC_HF_energy]
&BEGIN_PROVIDER [ double precision, TC_HF_one_e_energy]
&BEGIN_PROVIDER [ double precision, TC_HF_two_e_energy]

  BEGIN_DOC
  ! TC Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.
  END_DOC

  implicit none
  integer :: i, j

  PROVIDE mo_l_coef mo_r_coef

  TC_HF_energy = nuclear_repulsion
  TC_HF_one_e_energy = 0.d0
  TC_HF_two_e_energy = 0.d0

  do j = 1, ao_num
    do i = 1, ao_num
      TC_HF_two_e_energy += 0.5d0 * ( two_e_tc_non_hermit_integral_alpha(i,j) * TCSCF_density_matrix_ao_alpha(i,j) &
                                    + two_e_tc_non_hermit_integral_beta (i,j) * TCSCF_density_matrix_ao_beta (i,j) )
      TC_HF_one_e_energy += ao_one_e_integrals_tc_tot(i,j) &
                          * (TCSCF_density_matrix_ao_alpha(i,j) + TCSCF_density_matrix_ao_beta (i,j) )
    enddo
  enddo

  TC_HF_energy += TC_HF_one_e_energy + TC_HF_two_e_energy
  TC_HF_energy += diag_three_elem_hf

END_PROVIDER

! ---

