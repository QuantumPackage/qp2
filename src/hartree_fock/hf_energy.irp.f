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
 double precision :: tmp1, tmp2
 HF_energy = 0.d0
 HF_two_electron_energy = 0.d0
 HF_one_electron_energy = 0.d0
 do j=1,ao_num
   do i=1,ao_num
    tmp1 = 0.5d0 * ( ao_two_e_integral_alpha(i,j) * SCF_density_matrix_ao_alpha(i,j) &
                    +ao_two_e_integral_beta (i,j) * SCF_density_matrix_ao_beta (i,j) )
    tmp2 = ao_one_e_integrals(i,j) * (SCF_density_matrix_ao_alpha(i,j) + SCF_density_matrix_ao_beta (i,j) )
    HF_two_electron_energy += tmp1
    HF_one_electron_energy += tmp2
    HF_energy += tmp1 + tmp2
   enddo
 enddo
 HF_energy += nuclear_repulsion
END_PROVIDER


 BEGIN_PROVIDER [ double precision, HF_kinetic_energy]
&BEGIN_PROVIDER [ double precision, HF_n_e_energy]
 implicit none
 BEGIN_DOC
 ! Hartree-Fock energy containing the nuclear repulsion, and its one- and two-body components.
 END_DOC
 integer :: i,j
 double precision :: tmp1, tmp2
 HF_n_e_energy = 0.d0
 HF_kinetic_energy = 0.d0
 do j=1,ao_num
   do i=1,ao_num
    tmp1 = ao_integrals_n_e(i,j) * (SCF_density_matrix_ao_alpha(i,j) + SCF_density_matrix_ao_beta (i,j) )
    tmp2 = ao_kinetic_integrals(i,j) * (SCF_density_matrix_ao_alpha(i,j) + SCF_density_matrix_ao_beta (i,j) )
    HF_n_e_energy += tmp1
    HF_kinetic_energy += tmp2
   enddo
 enddo
END_PROVIDER

