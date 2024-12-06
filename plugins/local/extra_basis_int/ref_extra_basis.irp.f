program pouet
 implicit none
! call ref_overlap
! call ref_pot
! call ref_pot_mixed
! call routine_pot_ne_extra
! call ref_pot_ne_mixed
! call ref_pot_ne
 call ref_pot_ne_extra_mixed

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

subroutine ref_pot_mixed
 implicit none
 integer ::i,j
 double precision :: integral, C_center(3), mu_in
 double precision :: NAI_pol_mult_erf_ao
 C_center(1) = 0.1d0
 C_center(2) = -0.3d0
 C_center(3) =  0.8d0
 mu_in = 1.d10
 do i=1, 5
  do j = 6, ao_num
   integral = NAI_pol_mult_erf_ao(i, j, mu_in, C_center)
   write(34,*)integral
  enddo
 enddo

end

subroutine routine_pot_ne_extra
 implicit none
 integer :: i,j
 double precision :: v_extra_nucl_extra_ao
 do i = 1, ao_num
  do j = 1, ao_num
   write(34,*)ao_integrals_n_e(i,j)
  enddo
 enddo
end


subroutine ref_pot_ne_mixed
 implicit none
 integer ::i,j,k
 double precision :: integral
 do i=1, 5
  do j = 6, ao_num
   integral = 0.d0
   do k = 2,2
    integral += ao_integrals_n_e_per_atom(i,j,k) * nucl_charge(k)
   enddo
   write(34,*)integral
  enddo
 enddo
end

subroutine ref_pot_ne
 implicit none
 integer ::i,j,k
 double precision :: integral
 do i=6,ao_num
  do j = 6, ao_num
   integral = 0.d0
   do k = 1,1
    integral += ao_integrals_n_e_per_atom(i,j,k) * nucl_charge(k)
   enddo
   write(34,*)integral
  enddo
 enddo

end

subroutine ref_pot_ne_extra_mixed
 implicit none
 integer ::i,j,k
 double precision :: integral
 do i=1, 5
  do j = 6, ao_num
   integral = 0.d0
   do k = 1,1
    integral += ao_integrals_n_e_per_atom(i,j,k) * nucl_charge(k)
   enddo
   write(34,*)integral
  enddo
 enddo
end
