BEGIN_PROVIDER [double precision, extra_e_contrib_density]
 implicit none
 BEGIN_DOC
! Extra contribution to the SCF energy coming from the density.
!
! For a Hartree-Fock calculation: extra_e_contrib_density = 0
!
! For a Kohn-Sham or Range-separated Kohn-Sham: the exchange/correlation - trace of the V_xc potential
 END_DOC
 extra_e_contrib_density = 0.D0

END_PROVIDER

 BEGIN_PROVIDER [ double precision, HF_energy]
&BEGIN_PROVIDER [ double precision, HF_two_electron_energy]
&BEGIN_PROVIDER [ double precision, HF_one_electron_energy]
 implicit none
 BEGIN_DOC
 ! Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.
 END_DOC
 integer :: i,j
 HF_energy = nuclear_repulsion
 HF_two_electron_energy = 0.d0
 HF_one_electron_energy = 0.d0
 do j=1,ao_num
   do i=1,ao_num
    HF_two_electron_energy += 0.5d0 * ( ao_two_e_integral_alpha(i,j) * SCF_density_matrix_ao_alpha(i,j) &
                                       +ao_two_e_integral_beta(i,j)  * SCF_density_matrix_ao_beta(i,j) )
    HF_one_electron_energy += ao_one_e_integrals(i,j) * (SCF_density_matrix_ao_alpha(i,j) + SCF_density_matrix_ao_beta (i,j) )
   enddo
 enddo
 HF_energy += HF_two_electron_energy + HF_one_electron_energy
END_PROVIDER

