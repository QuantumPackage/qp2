program pouet
 implicit none
! call ref_overlap
 call ref_pot

end

subroutine ref_overlap
 implicit none
 integer :: i,j
 do i = 1, ao_num
  do j = 1, ao_num
   write(34,*)ao_overlap(j,i)
  enddo
 enddo

end

subroutine ref_pot
 implicit none
 integer :: i,j
 double precision :: integral, C_center(3), mu_in
 double precision :: NAI_pol_mult_erf_ao
 C_center(1) = 0.1d0
 C_center(2) = -0.3d0
 C_center(3) =  0.8d0
 mu_in = 1.d10
 do i = 1, ao_num
  do j = 1, ao_num
   integral = NAI_pol_mult_erf_ao(i, j, mu_in, C_center)
   write(34,*)j,i,integral
  enddo
 enddo

end
