program extra_basis_int
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
!  call test_overlap
!  call routine_test_pot_ne
! call routine_test_pot_ne_mixed
! call routine_pot_ne_extra
! call routine_test_pot_ne_mixed
! call routine_pot_ne
 call routine_test_pot_ne_extra_mixed
 
end

subroutine test_overlap
 implicit none
  integer :: i,j
  do i = 1, ao_extra_num
   do j = 1, ao_extra_num 
    write(33,*)ao_extra_overlap(j,i)
   enddo
  enddo
end

subroutine test_overlap_mixed
 implicit none
  integer :: i,j
  double precision, allocatable :: ao_mixed_overlap(:,:)
  allocate(ao_mixed_overlap(ao_extra_num,ao_num))
  call get_ao_mixed_overlap(extra_nucl_coord,ao_mixed_overlap)
  do i = 1, ao_extra_num
   do j = 1, ao_num 
    write(33,*)dabs(ao_extra_overlap_mixed(j,i)-ao_mixed_overlap(i,j))
    write(*,*)ao_extra_overlap_mixed(j,i),ao_mixed_overlap(i,j),dabs(ao_extra_overlap_mixed(j,i)-ao_mixed_overlap(i,j))
   enddo
  enddo
end

subroutine  routine_test_pot_ne
 implicit none
 integer :: i,j
 double precision :: integral, C_center(3), mu_in
 double precision :: NAI_pol_mult_erf_ao_extra
 C_center(1) = 0.1d0
 C_center(2) = -0.3d0
 C_center(3) =  0.8d0
 mu_in = 1.d10
 do i = 1, ao_extra_num
  do j = 1, ao_extra_num
   integral = NAI_pol_mult_erf_ao_extra(i, j, mu_in, C_center)
   write(33,*)j,i,integral
  enddo
 enddo

end

subroutine  routine_test_pot_mixed
 implicit none
 integer :: i,j
 double precision :: integral, C_center(3), mu_in
 double precision :: NAI_pol_mult_erf_ao_extra_mixed
 C_center(1) = 0.1d0
 C_center(2) = -0.3d0
 C_center(3) =  0.8d0
 mu_in = 1.d10
 do j = 1, ao_num
  do i = 1, ao_extra_num
   integral = NAI_pol_mult_erf_ao_extra_mixed(i, j, mu_in, C_center)
   write(33,*)integral
  enddo
 enddo

end

subroutine routine_pot_ne_extra
 implicit none
 integer :: i,j
 double precision :: v_extra_nucl_extra_ao
 do i = 1, ao_extra_num
  do j = 1, ao_extra_num
   write(33,*)v_extra_nucl_extra_ao(i,j)
  enddo
 enddo
end


subroutine routine_pot_ne
 implicit none
 integer :: i,j
 double precision :: v_nucl_extra_ao
 do i = 1, ao_extra_num
  do j = 1, ao_extra_num
   write(33,*)v_nucl_extra_ao(i,j)
  enddo
 enddo
end


subroutine  routine_test_pot_ne_mixed
 implicit none
 integer :: i,j
 double precision :: integral,v_extra_nucl_mixed_ao
 do j = 1, ao_num
  do i = 1, ao_extra_num
   integral = v_extra_nucl_mixed_ao(i,j)
   write(33,*)integral
  enddo
 enddo

end

subroutine  routine_test_pot_ne_extra_mixed
 implicit none
 integer :: i,j
 double precision :: integral,v_nucl_mixed_ao
 do j = 1, ao_num
  do i = 1, ao_extra_num
   integral = v_nucl_mixed_ao(i,j)
   write(33,*)integral
  enddo
 enddo

end
