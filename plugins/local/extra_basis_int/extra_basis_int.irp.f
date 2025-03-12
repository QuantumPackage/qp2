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
! call routine_test_pot_ne_extra_mixed
! call routine_test_coul_1s
! call print_v_ne_extra_basis
! call print_v_ne_basis
! call test_v_ne_a_extra_basis
 call print_v_ee_mixed_direct
 
end

subroutine test_v_ne_a_extra_basis
 implicit none
 integer :: i,j
 do i = 1, ao_extra_num
  write(*,'(100(F16.10,X))')pot_vne_A_extra_basis(1:ao_extra_num,i)
 enddo
end


subroutine test_overlap
 implicit none
  integer :: i,j
  do i = 1, ao_num
!   do j = 1, ao_num 
    write(33,'(100(F16.10,X))')ao_extra_overlap_mixed(i,1:ao_extra_num)
!   enddo
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

subroutine routine_test_coul_1s
 implicit none
 integer :: i,j
 double precision :: r(3) ,mu_in,NAI_pol_mult_erf_ao_extra
 double precision :: ref,new, accu,coul_full_ao_pq_r_1s,v_nucl_extra_ao
 r(1) = 0.d0
 r(2) = 0.5d0
 r(3) = -1.5d0
 r=nucl_coord(1,1:3)
 mu_in = 1.d+10
 accu = 0.d0
 do i = 1, ao_extra_num
  do j = 1, ao_extra_num
! do i = 1, 1
!  do j = 1, 1
   ref = NAI_pol_mult_erf_ao_extra(i, j, mu_in, r)
   new = coul_full_ao_pq_r_1s(i,j,r,ao_extra_center_1s(1,i),ao_extra_center_1s(1,j))
!   new = v_nucl_extra_ao(i,j)
   accu += dabs(new-ref)
  enddo
 enddo
 print*,'accu = ',accu
end

subroutine print_v_ne_extra_basis
 implicit none
 integer :: i_ao,j_ao,i
 double precision :: integral, accu, charge, coord(3), coul_pq_r_1s

 accu = 0.d0
 do i_ao = 1, ao_extra_num
  do j_ao = 1, ao_extra_num
   do i = 1, nucl_num
    charge = nucl_charge(i)
    coord(1:3) = nucl_coord_transp(1:3,i)
    integral = coul_pq_r_1s(i_ao,j_ao,coord,ao_extra_center_1s(1,i_ao),ao_extra_center_1s(1,j_ao))
    accu += -charge * integral * effective_ao_extra_dm(j_ao,i_ao) 
   enddo
  enddo
 enddo
 print*,'accu = ',accu

end

subroutine print_v_ne_basis
 implicit none
 integer :: i_ao,j_ao,i
 double precision :: integral, accu, charge, coord(3), coul_full_pq_r_1s,NAI_pol_mult_erf_ao

 accu = 0.d0
 do i_ao = 1, ao_num
  do j_ao = 1, ao_num
   do i = 1, nucl_num
    charge = nucl_charge(i)
    coord(1:3) = nucl_coord_transp(1:3,i)
    integral = NAI_pol_mult_erf_ao(i_ao, j_ao, 1d+10, coord)
    accu += -charge * integral * one_e_dm_ao(j_ao,i_ao) 
   enddo
  enddo
 enddo
 print*,'accu = ',accu

end

subroutine print_v_ee_mixed_direct
 implicit none
 integer :: i,j,k,l
 double precision :: ao_two_e_integral_mixed_direct
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_extra_num
    do l = 1, ao_extra_num
     write(34,*)ao_two_e_integral_mixed_direct(i, j, k, l)
    enddo
   enddo
  enddo
 enddo

end
