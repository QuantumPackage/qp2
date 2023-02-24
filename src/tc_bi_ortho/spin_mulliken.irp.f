
BEGIN_PROVIDER [double precision, tc_spin_population, (ao_num,ao_num,N_states)]
 implicit none
 integer :: i,j,istate
 BEGIN_DOC
! spin population on the ao basis :
! tc_spin_population(i,j) = rho_AO(alpha)(i,j) - rho_AO(beta)(i,j) * <AO_i|AO_j>
 END_DOC
 tc_spin_population = 0.d0
 if(only_spin_tc_right)then
  do i = 1, ao_num
   do j = 1, ao_num
    tc_spin_population(j,i,1) = tc_spin_dens_right_only(j,i) * ao_overlap(j,i)
   enddo
  enddo
 else
  do istate = 1, N_states
   do i = 1, ao_num
    do j = 1, ao_num
     tc_spin_population(j,i,istate) = tc_spin_transition_matrix_ao(j,i,istate,istate) * ao_overlap(j,i)
    enddo
   enddo
  enddo
 endif
END_PROVIDER

 BEGIN_PROVIDER [double precision, tc_spin_population_angular_momentum, (0:ao_l_max,N_states)]
&BEGIN_PROVIDER [double precision, tc_spin_population_angular_momentum_per_atom, (0:ao_l_max,nucl_num,N_states)]
 implicit none
 integer :: i,istate
 double precision :: accu
 tc_spin_population_angular_momentum = 0.d0
 tc_spin_population_angular_momentum_per_atom = 0.d0
 do istate = 1, N_states
  do i = 1,  ao_num
   tc_spin_population_angular_momentum(ao_l(i),istate) += tc_spin_gross_orbital_product(i,istate)
   tc_spin_population_angular_momentum_per_atom(ao_l(i),ao_nucl(i),istate) += tc_spin_gross_orbital_product(i,istate)
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [double precision, tc_spin_gross_orbital_product, (ao_num,N_states)]
 implicit none
 tc_spin_gross_orbital_product = 0.d0
 integer :: i,j,istate
 BEGIN_DOC
! gross orbital product for the spin population
 END_DOC
 do istate = 1, N_states
  do i = 1, ao_num
   do j = 1, ao_num
    tc_spin_gross_orbital_product(i,istate) += tc_spin_population(j,i,istate)
   enddo
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, tc_mulliken_spin_densities, (nucl_num,N_states)]
 implicit none
 integer :: i,j,istate
 BEGIN_DOC
!ATOMIC SPIN POPULATION (ALPHA MINUS BETA)
 END_DOC
 tc_mulliken_spin_densities = 0.d0
 do istate = 1, N_states
  do i = 1, ao_num
   tc_mulliken_spin_densities(ao_nucl(i),istate) += tc_spin_gross_orbital_product(i,istate)
  enddo
 enddo

END_PROVIDER

subroutine tc_print_mulliken_sd
 implicit none
 double precision :: accu
 integer :: i
 integer :: j
 print*,'Mulliken spin densities'
 accu= 0.d0
 do i = 1, nucl_num
  print*,i,nucl_charge(i),tc_mulliken_spin_densities(i,1)
  accu += tc_mulliken_spin_densities(i,1)
 enddo
 print*,'Sum of Mulliken SD = ',accu
 print*,'AO SPIN POPULATIONS'
 accu = 0.d0
 do i = 1, ao_num
  accu += tc_spin_gross_orbital_product(i,1)
  write(*,'(1X,I3,1X,A4,1X,I2,1X,A4,1X,F10.7)')i,trim(element_name(int(nucl_charge(ao_nucl(i))))),ao_nucl(i),trim(l_to_character(ao_l(i))),tc_spin_gross_orbital_product(i,1)
 enddo
 print*,'sum = ',accu
 accu = 0.d0
 print*,'Angular momentum analysis'
 do i = 0,  ao_l_max
  accu += tc_spin_population_angular_momentum(i,1)
  print*,' ',trim(l_to_character(i)),tc_spin_population_angular_momentum(i,1)
 print*,'sum = ',accu
 enddo
 print*,'Angular momentum analysis per atom'
 print*,'Angular momentum analysis'
 do j = 1,nucl_num
  accu = 0.d0
  do i = 0,  ao_l_max
   accu += tc_spin_population_angular_momentum_per_atom(i,j,1)
   write(*,'(1X,I3,1X,A4,1X,A4,1X,F10.7)')j,trim(element_name(int(nucl_charge(j)))),trim(l_to_character(i)),tc_spin_population_angular_momentum_per_atom(i,j,1)
   print*,'sum = ',accu
  enddo
 enddo

end

