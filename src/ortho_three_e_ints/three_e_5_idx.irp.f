
BEGIN_PROVIDER [ double precision, three_body_5_index, (mo_num, mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 5 index matrix element of the -L  three-body operator 
!
! three_body_5_index(i,j,m,l,n) = < phi_i phi_j phi_m | phi_i phi_l phi_n >
!
! notice the -1 sign: in this way three_body_5_index can be directly used to compute Slater rules :)
 END_DOC
 integer :: j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 three_body_5_index(1:mo_num, 1:mo_num, 1:mo_num, 1:mo_num, 1:mo_num) = 0.d0
 print*,'Providing the three_body_5_index ...'
 name_file = 'three_body_5_index'
 call wall_time(wall0)
 if(read_three_body_ints)then
  print*,'Reading three_body_5_index from disk ...'
  call read_array_5_index_tensor(mo_num,three_body_5_index,name_file)
 else
  provide x_W_ij_erf_rk
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (j,k,l,m,n,integral) & 
  !$OMP SHARED (mo_num,three_body_5_index)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do n = 1, mo_num
    do l = 1, mo_num
     do k = 1, mo_num
!      do m = n, mo_num
!       do j = l, mo_num
      do m = 1, mo_num
       do j = 1, mo_num
         integral = 0.d0
         
         call give_integrals_3_body(j,m,k,l,n,k,integral)
 
         three_body_5_index(k,j,m,l,n) = -1.d0 * integral 
   
       enddo
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_5_index',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_5_index on disk ...' 
  call write_array_5_index_tensor(mo_num,three_body_5_index,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif
! do n = 1, mo_num
!  do l = 1, mo_num
!   do k = 1, mo_num
!    do m = 1, n-1
!     do j = 1, l-1
!      three_body_5_index(k,j,m,l,n) = three_body_5_index(k,l,n,j,m)
!      three_body_5_index(k,j,m,l,n)
!     enddo
!    enddo
!   enddo
!  enddo
! enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_5_index_exch_13, (mo_num, mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 5 index matrix element of the -L  three-body operator 
!
! three_body_5_index_exch_13(k,j,m,l,n) = < phi_j phi_m phi_k | phi_k phi_n phi_l >
!
! notice the -1 sign: in this way three_body_5_index_exch_13 can be directly used to compute Slater rules :)
 END_DOC
 integer :: j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 
 three_body_5_index_exch_13 = 0.d0

 name_file = 'three_body_5_index_exch_13'
 print*,'Providing the three_body_5_index_exch_13 ...'
 call wall_time(wall0)
 if(read_three_body_ints)then
  print*,'Reading three_body_5_index_exch_13 from disk ...'
  call read_array_5_index_tensor(mo_num,three_body_5_index_exch_13,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (j,k,l,m,n,integral) & 
  !$OMP SHARED (mo_num,three_body_5_index_exch_13)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do n = 1, mo_num
    do l = 1, mo_num
     do k = 1, mo_num
!      do m = n, mo_num
!       do j = l, mo_num
      do m = 1, mo_num
       do j = 1, mo_num
         integral = 0.d0
!!                                  j,m,k,l,n,k : direct (case 2)
         call give_integrals_3_body(j,m,k,k,n,l,integral)
!!                                  j,m,k,k,n,l : exchange 1 3
 
         three_body_5_index_exch_13(k,j,m,l,n) = -1.d0 * integral 
       enddo
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_5_index_exch_13',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_5_index_exch_13 on disk ...' 
  call write_array_5_index_tensor(mo_num,three_body_5_index_exch_13,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif
! do n = 1, mo_num
!  do l = 1, mo_num
!   do k = 1, mo_num
!    do m = n, mo_num
!     do j = l, mo_num
!      three_body_5_index_exch_13(k,l,n,j,m) = three_body_5_index_exch_13(k,j,m,l,n)
!     enddo
!    enddo
!   enddo
!  enddo
! enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_5_index_exch_32, (mo_num, mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 5 index matrix element of the -L  three-body operator 
!
! three_body_5_index_exch_32(i,j,m,l,n) = < phi_i phi_j phi_m | phi_i phi_l phi_n >
!
! notice the -1 sign: in this way three_body_5_index_exch_32 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(328) :: name_file 
 
 three_body_5_index_exch_32 = 0.d0
 name_file = 'three_body_5_index_exch_32'
 print*,'Providing the three_body_5_index_exch_32 ...'
 call wall_time(wall0)

 if(read_three_body_ints)then
  print*,'Reading three_body_5_index_exch_32 from disk ...'
  call read_array_5_index_tensor(mo_num,three_body_5_index_exch_32,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (j,k,l,m,n,integral) & 
  !$OMP SHARED (mo_num,three_body_5_index_exch_32)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do n = 1, mo_num
    do l = 1, mo_num
     do k = 1, mo_num
!      do m = n, mo_num
!       do j = l, mo_num
      do m = 1, mo_num
       do j = 1, mo_num
         integral = 0.d0
!!                                  j,m,k,l,n,k : direct (case 3)
         call give_integrals_3_body(j,m,k,l,k,n,integral)
!!                                  j,m,k,l,k,n : exchange 2 3
 
         three_body_5_index_exch_32(k,j,m,l,n) = -1.d0 * integral 
   
       enddo
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_5_index_exch_32',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_5_index_exch_32 on disk ...' 
  call write_array_5_index_tensor(mo_num,three_body_5_index_exch_32,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif
! do n = 1, mo_num
!  do l = 1, mo_num
!   do k = 1, mo_num
!    do m = n, mo_num
!     do j = l, mo_num
!      three_body_5_index_exch_32(k,l,n,j,m) = three_body_5_index_exch_32(k,j,m,l,n)
!     enddo
!    enddo
!   enddo
!  enddo
! enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_5_index_exch_12, (mo_num, mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 5 index matrix element of the -L  three-body operator 
!
! three_body_5_index_exch_12(i,j,m,l,n) = < phi_i phi_j phi_m | phi_i phi_l phi_n >
!
! notice the -1 sign: in this way three_body_5_index_exch_12 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(328) :: name_file 

 three_body_5_index_exch_12 = 0.d0
 name_file = 'three_body_5_index_exch_12'
 print*,'Providing the three_body_5_index_exch_12 ...'
 call wall_time(wall0)

 if(read_three_body_ints)then
  print*,'Reading three_body_5_index_exch_12 from disk ...'
  call read_array_5_index_tensor(mo_num,three_body_5_index_exch_12,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (j,k,l,m,n,integral) & 
  !$OMP SHARED (mo_num,three_body_5_index_exch_12)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do n = 1, mo_num
    do l = 1, mo_num
     do k = 1, mo_num
!      do m = n, mo_num
!       do j = l, mo_num
      do m = 1, mo_num
       do j = 1, mo_num
         integral = 0.d0
!!                                  j,m,k,l,n,k : direct (case 1)
         call give_integrals_3_body(j,m,k,n,l,k,integral)
!!                                  j,m,k,l,k,n : exchange 2 3
 
         three_body_5_index_exch_12(k,j,m,l,n) = -1.d0 * integral 
   
       enddo
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_5_index_exch_12',wall1 - wall0
! do n = 1, mo_num
!  do l = 1, mo_num
!   do k = 1, mo_num
!    do m = n, mo_num
!     do j = l, mo_num
!      three_body_5_index_exch_12(k,l,n,j,m) = three_body_5_index_exch_12(k,j,m,l,n)
!     enddo
!    enddo
!   enddo
!  enddo
! enddo
 if(write_three_body_ints)then
  print*,'Writing three_body_5_index_exch_12 on disk ...' 
  call write_array_5_index_tensor(mo_num,three_body_5_index_exch_12,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 


BEGIN_PROVIDER [ double precision, three_body_5_index_312, (mo_num, mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 5 index matrix element of the -L  three-body operator 
!
! three_body_5_index_312(i,j,m,l,n) = < phi_i phi_j phi_m | phi_i phi_l phi_n >
!
! notice the -1 sign: in this way three_body_5_index_312 can be directly used to compute Slater rules :)
 END_DOC
 integer :: j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 
 three_body_5_index_312 = 0.d0
 name_file = 'three_body_5_index_312'
 print*,'Providing the three_body_5_index_312 ...'
 call wall_time(wall0)

 if(read_three_body_ints)then
  print*,'Reading three_body_5_index_312 from disk ...'
  call read_array_5_index_tensor(mo_num,three_body_5_index_312,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (j,k,l,m,n,integral) & 
  !$OMP SHARED (mo_num,three_body_5_index_312)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do n = 1, mo_num
    do l = 1, mo_num
     do k = 1, mo_num
!      do m = n, mo_num
!       do j = l, mo_num
      do m = 1, mo_num
       do j = 1, mo_num
         integral = 0.d0
         
         !                         <j m k|l n k> - > <j m k|n k l>
         call give_integrals_3_body(j,m,k,n,k,l,integral)
 
         three_body_5_index_312(k,j,m,l,n) = -1.d0 * integral 
   
       enddo
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_5_index_312',wall1 - wall0
! do n = 1, mo_num
!  do l = 1, mo_num
!   do k = 1, mo_num
!    do m = n, mo_num
!     do j = l, mo_num
!      three_body_5_index_312(k,l,n,j,m) = three_body_5_index_312(k,j,m,l,n)
!     enddo
!    enddo
!   enddo
!  enddo
! enddo
 if(write_three_body_ints)then
  print*,'Writing three_body_5_index_312 on disk ...' 
  call write_array_5_index_tensor(mo_num,three_body_5_index_312,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_5_index_132, (mo_num, mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 5 index matrix element of the -L  three-body operator 
!
! three_body_5_index_132(i,j,m,l,n) = < phi_i phi_j phi_m | phi_i phi_l phi_n >
!
! notice the -1 sign: in this way three_body_5_index_132 can be directly used to compute Slater rules :)
 END_DOC
 integer :: j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 three_body_5_index_132 = 0.d0
 name_file = 'three_body_5_index_132'
 print*,'Providing the three_body_5_index_132 ...'
 call wall_time(wall0)

 if(read_three_body_ints)then
  print*,'Reading three_body_5_index_132 from disk ...'
  call read_array_5_index_tensor(mo_num,three_body_5_index_132,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (j,k,l,m,n,integral) & 
  !$OMP SHARED (mo_num,three_body_5_index_132)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do n = 1, mo_num
    do l = 1, mo_num
     do k = 1, mo_num
!      do m = n, mo_num
!       do j = l, mo_num
      do m = 1, mo_num
       do j = 1, mo_num
         integral = 0.d0
         
         !                         <j m k|l n k> - > <j m k|k l n>
         call give_integrals_3_body(j,m,k,k,l,n,integral)
 
         three_body_5_index_132(k,j,m,l,n) = -1.d0 * integral 
   
       enddo
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_5_index_132',wall1 - wall0
! do n = 1, mo_num
!  do l = 1, mo_num
!   do k = 1, mo_num
!    do m = n, mo_num
!     do j = l, mo_num
!      three_body_5_index_132(k,l,n,j,m) = three_body_5_index_132(k,j,m,l,n)
!     enddo
!    enddo
!   enddo
!  enddo
! enddo
 if(write_three_body_ints)then
  print*,'Writing three_body_5_index_132 on disk ...' 
  call write_array_5_index_tensor(mo_num,three_body_5_index_132,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 

subroutine write_array_5_index_tensor(n_orb,array_tmp,name_file)
 implicit none
 integer, intent(in) :: n_orb
 character*(128),  intent(in) :: name_file 
 double precision, intent(in) :: array_tmp(n_orb,n_orb,n_orb,n_orb,n_orb)

 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 PROVIDE ezfio_filename                                                                                                  
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'W')
 write(i_unit_output)array_tmp
 close(unit=i_unit_output)
end

subroutine read_array_5_index_tensor(n_orb,array_tmp,name_file)
 implicit none
 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 integer, intent(in) :: n_orb
 character*(128),  intent(in)  :: name_file 
 double precision, intent(out) :: array_tmp(n_orb,n_orb,n_orb,n_orb,n_orb)
 PROVIDE ezfio_filename                                                                                                  
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'R')
 read(i_unit_output)array_tmp
 close(unit=i_unit_output)
end
