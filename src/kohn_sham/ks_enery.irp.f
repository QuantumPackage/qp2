 BEGIN_PROVIDER [ double precision, KS_energy]
&BEGIN_PROVIDER [ double precision, two_e_energy]
&BEGIN_PROVIDER [ double precision, one_e_energy]
&BEGIN_PROVIDER [ double precision, Fock_matrix_energy]
&BEGIN_PROVIDER [ double precision, trace_potential_xc ]
 implicit none
 BEGIN_DOC
 ! Kohn-Sham energy containing the nuclear repulsion energy, and the various components of this quantity.
 END_DOC

 integer                        :: i,j
 double precision :: accu_mono,accu_fock
 KS_energy = nuclear_repulsion
 one_e_energy = 0.d0
 two_e_energy = 0.d0
 Fock_matrix_energy = 0.d0
 trace_potential_xc = 0.d0
 do j=1,ao_num
   do i=1,ao_num
    Fock_matrix_energy +=   Fock_matrix_ao_alpha(i,j) * SCF_density_matrix_ao_alpha(i,j) + &
                            Fock_matrix_ao_beta(i,j) * SCF_density_matrix_ao_beta(i,j)
    two_e_energy += 0.5d0 * ( ao_two_e_integral_alpha(i,j) * SCF_density_matrix_ao_alpha(i,j) &
                +ao_two_e_integral_beta(i,j) * SCF_density_matrix_ao_beta(i,j) )
    one_e_energy += ao_one_e_integrals(i,j) * (SCF_density_matrix_ao_alpha(i,j) + SCF_density_matrix_ao_beta (i,j) )
    trace_potential_xc += ao_potential_alpha_xc(i,j) * SCF_density_matrix_ao_alpha(i,j) + ao_potential_beta_xc(i,j) *  SCF_density_matrix_ao_beta (i,j)
   enddo
 enddo

 KS_energy +=  e_exchange_dft + e_correlation_dft + one_e_energy + two_e_energy
END_PROVIDER

BEGIN_PROVIDER [double precision, extra_e_contrib_density]
 implicit none
 BEGIN_DOC
! Extra contribution to the SCF energy coming from the density.
!
! For a Hartree-Fock calculation: extra_e_contrib_density = 0
!
! For a Kohn-Sham or Range-separated Kohn-Sham: the exchange/correlation - 1/2 trace of the V_xc potential
 END_DOC
 extra_e_contrib_density = e_exchange_dft + e_correlation_dft - 0.5d0 * trace_potential_xc
END_PROVIDER

