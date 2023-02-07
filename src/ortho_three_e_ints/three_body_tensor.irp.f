BEGIN_PROVIDER [ double precision, three_body_ints, (mo_num, mo_num, mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! matrix element of the -L  three-body operator 
!
! notice the -1 sign: in this way three_body_ints can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 three_body_ints = 0.d0
 print*,'Providing the three_body_ints ...'
 call wall_time(wall0)
 name_file = 'six_index_tensor'
 if(read_three_body_ints)then
  call read_fcidump_3_tc(three_body_ints)
 else
  if(read_three_body_ints)then
   print*,'Reading three_body_ints from disk ...'
   call read_array_6_index_tensor(mo_num,three_body_ints,name_file)
  else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,k,l,m,n,integral) & 
  !$OMP SHARED (mo_num,three_body_ints)
  !$OMP DO SCHEDULE (dynamic)
   do n = 1, mo_num
    do l = 1, mo_num
     do k = 1, mo_num
      do m = n, mo_num
       do j = l, mo_num
        do i = k, mo_num
!!         if(i>=j)then
           integral = 0.d0
           call give_integrals_3_body(i,j,m,k,l,n,integral)
 
           three_body_ints(i,j,m,k,l,n) = -1.d0 * integral 
   
           ! permutation with k,i
           three_body_ints(k,j,m,i,l,n) = -1.d0 * integral ! i,k
           ! two permutations with k,i
           three_body_ints(k,l,m,i,j,n) = -1.d0 * integral 
           three_body_ints(k,j,n,i,l,m) = -1.d0 * integral 
           ! three permutations with k,i
           three_body_ints(k,l,n,i,j,m) = -1.d0 * integral 
   
           ! permutation with l,j
           three_body_ints(i,l,m,k,j,n) = -1.d0 * integral ! j,l
           ! two permutations with l,j
           three_body_ints(k,l,m,i,j,n) = -1.d0 * integral 
           three_body_ints(i,l,n,k,j,m) = -1.d0 * integral 
           ! two permutations with l,j
!!!!        three_body_ints(k,l,n,i,j,m) = -1.d0 * integral 
   
           ! permutation with m,n
           three_body_ints(i,j,n,k,l,m) = -1.d0 * integral ! m,n
           ! two permutations with m,n
           three_body_ints(k,j,n,i,l,m) = -1.d0 * integral ! m,n
           three_body_ints(i,l,n,k,j,m) = -1.d0 * integral ! m,n
           ! three permutations with k,i
!!!!        three_body_ints(k,l,n,i,j,m) = -1.d0 * integral ! m,n
 
!!         endif
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
  endif
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_ints',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_ints on disk ...'
  call write_array_6_index_tensor(mo_num,three_body_ints,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 



subroutine give_integrals_3_body(i,j,m,k,l,n,integral)
 implicit none
 double precision, intent(out) :: integral
 integer, intent(in) :: i,j,m,k,l,n
 double precision :: weight
 BEGIN_DOC
! <ijm|L|kln>
 END_DOC
 integer :: ipoint,mm
 integral = 0.d0
 do mm = 1, 3
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)                                                                          
   integral += weight * mos_in_r_array_transp(ipoint,i) * mos_in_r_array_transp(ipoint,k) * x_W_ij_erf_rk(ipoint,mm,m,n) * x_W_ij_erf_rk(ipoint,mm,j,l) 
   integral += weight * mos_in_r_array_transp(ipoint,j) * mos_in_r_array_transp(ipoint,l) * x_W_ij_erf_rk(ipoint,mm,m,n) * x_W_ij_erf_rk(ipoint,mm,i,k) 
   integral += weight * mos_in_r_array_transp(ipoint,m) * mos_in_r_array_transp(ipoint,n) * x_W_ij_erf_rk(ipoint,mm,j,l) * x_W_ij_erf_rk(ipoint,mm,i,k) 
  enddo
 enddo
end

