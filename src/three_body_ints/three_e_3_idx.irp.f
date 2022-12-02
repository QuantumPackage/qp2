
BEGIN_PROVIDER [ double precision, three_body_3_index, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix element of the -L  three-body operator 
!
! three_body_3_index(k,l,n) = < phi_k phi_l phi_n | phi_k phi_l phi_n >
!
! notice the -1 sign: in this way three_body_3_index can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 print*,'Providing the three_body_3_index ...'
 name_file = 'three_body_3_index'
 call wall_time(wall0)
 if(read_three_body_ints)then
  print*,'Reading three_body_ints from disk ...'
  call read_array_3_index_tensor(mo_num,three_body_3_index,name_file)
 else
  provide x_W_ij_erf_rk
  three_body_3_index = 0.d0
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,integral) & 
  !$OMP SHARED (mo_num,three_body_3_index)
  !$OMP DO SCHEDULE (guided) COLLAPSE(3)
   do m = 1, mo_num ! 3
    do j = 1, mo_num ! 2 
     do i = 1, mo_num ! 1 
      integral = 0.d0
      !                          1 2 3 1 2 3
      call give_integrals_3_body(i,j,m,i,j,m,integral)
 
      three_body_3_index(i,j,m) = -1.d0 * integral 
   
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_3_index',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_3_index on disk ...' 
  call write_array_3_index_tensor(mo_num,three_body_3_index,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif


END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_3_index_exch_12, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix EXCHANGE element of the -L  three-body operator 
!
! three_body_3_index_exch_12(k,l,n) = < phi_k phi_l phi_n | phi_l phi_k phi_n >
!
! notice the -1 sign: in this way three_body_3_index_exch_12 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 name_file = 'three_body_3_index_exch_12'
 print*,'Providing the three_body_3_index_exch_12 ...'
 call wall_time(wall0)
 if(read_three_body_ints)then
  print*,'Reading three_body_ints from disk ...'
  call read_array_3_index_tensor(mo_num,three_body_3_index_exch_12,name_file)
 else
  provide x_W_ij_erf_rk
  three_body_3_index_exch_12 = 0.d0
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,integral) & 
  !$OMP SHARED (mo_num,three_body_3_index_exch_12)
  !$OMP DO SCHEDULE (guided) COLLAPSE(3)
   do m = 1, mo_num ! 3
    do j = 1, mo_num ! 2 
     do i = 1, mo_num ! 1 
      integral = 0.d0
      !                          1 2 3 1 2 3
      call give_integrals_3_body(i,j,m,j,i,m,integral)
 
      three_body_3_index_exch_12(i,j,m) = -1.d0 * integral 
   
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_3_index_exch_12',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_3_index_exch_12 on disk ...' 
  call write_array_3_index_tensor(mo_num,three_body_3_index_exch_12,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif
END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_3_index_exch_23, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix EXCHANGE element of the -L  three-body operator 
!
! three_body_3_index_exch_12(k,l,n) = < phi_k phi_l phi_n | phi_k phi_n phi_l >
!
! notice the -1 sign: in this way three_body_3_index_exch_12 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 print*,'Providing the three_body_3_index_exch_23 ...'
 call wall_time(wall0)
 name_file = 'three_body_3_index_exch_23'
 if(read_three_body_ints)then
  print*,'Reading three_body_ints from disk ...'
  call read_array_3_index_tensor(mo_num,three_body_3_index_exch_23,name_file)
 else
  provide x_W_ij_erf_rk
  three_body_3_index_exch_23 = 0.d0
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,integral) & 
  !$OMP SHARED (mo_num,three_body_3_index_exch_23)
  !$OMP DO SCHEDULE (guided) COLLAPSE(3)
   do m = 1, mo_num ! 3
    do j = 1, mo_num ! 2 
     do i = 1, mo_num ! 1 
      integral = 0.d0
      !                          1 2 3 1 2 3
      call give_integrals_3_body(i,j,m,i,m,j,integral)
 
      three_body_3_index_exch_23(i,j,m) = -1.d0 * integral 
   
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
  call wall_time(wall1)
 endif
 print*,'wall time for three_body_3_index_exch_23',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_3_index_exch_23 on disk ...' 
  call write_array_3_index_tensor(mo_num,three_body_3_index_exch_23,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_3_index_exch_13, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix EXCHANGE element of the -L  three-body operator 
!
! three_body_3_index_exch_12(k,l,n) = < phi_k phi_l phi_n | phi_k phi_n phi_l >
!
! notice the -1 sign: in this way three_body_3_index_exch_12 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 print*,'Providing the three_body_3_index_exch_13 ...'
 call wall_time(wall0)
 name_file = 'three_body_3_index_exch_13'
 if(read_three_body_ints)then
  print*,'Reading three_body_ints from disk ...'
  call read_array_3_index_tensor(mo_num,three_body_3_index_exch_13,name_file)
 else
  provide x_W_ij_erf_rk
  three_body_3_index_exch_13 = 0.d0
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,integral) & 
  !$OMP SHARED (mo_num,three_body_3_index_exch_13)
  !$OMP DO SCHEDULE (guided)
   do m = 1, mo_num ! 3
    do j = 1, mo_num ! 2 
     do i = 1, mo_num ! 1 
      integral = 0.d0
      !                          1 2 3 1 2 3
      call give_integrals_3_body(i,j,m,m,j,i,integral)
 
      three_body_3_index_exch_13(i,j,m) = -1.d0 * integral 
   
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif

 call wall_time(wall1)
 print*,'wall time for three_body_3_index_exch_13',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_3_index_exch_13 on disk ...' 
  call write_array_3_index_tensor(mo_num,three_body_3_index_exch_13,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif
END_PROVIDER 


BEGIN_PROVIDER [ double precision, three_body_3_index_exch_231, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix element of the -L  three-body operator 
!
! three_body_3_index_exch_231(k,l,n) = < phi_k phi_l phi_n | phi_l phi_n phi_k >
!
! notice the -1 sign: in this way three_body_3_index can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 print*,'Providing the three_body_3_index_231 ...'
 call wall_time(wall0)
 name_file = 'three_body_3_index_exch_231'
 if(read_three_body_ints)then
  print*,'Reading three_body_ints from disk ...'
  call read_array_3_index_tensor(mo_num,three_body_3_index_exch_231,name_file)
 else
  provide x_W_ij_erf_rk
  three_body_3_index_exch_231 = 0.d0
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,integral) & 
  !$OMP SHARED (mo_num,three_body_3_index_exch_231)
  !$OMP DO SCHEDULE (guided) COLLAPSE(3)
   do m = 1, mo_num ! 3
    do j = 1, mo_num ! 2 
     do i = 1, mo_num ! 1 
      integral = 0.d0
      !                          1 2 3 1 2 3
      call give_integrals_3_body(i,j,m,j,m,i,integral)
 
      three_body_3_index_exch_231(i,j,m) = -1.d0 * integral 
   
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_3_index_exch_231 ',wall1 - wall0

 if(write_three_body_ints)then
  print*,'Writing three_body_3_index_exch_231 on disk ...' 
  call write_array_3_index_tensor(mo_num,three_body_3_index_exch_231,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_3_index_exch_312, (mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 3 index matrix element of the -L  three-body operator 
!
! three_body_3_index(k,l,n) = < phi_k phi_l phi_n | phi_l phi_n phi_k >
!
! notice the -1 sign: in this way three_body_3_index can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,m
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 print*,'Providing the three_body_3_index_312 ...'
 call wall_time(wall0)
 name_file = 'three_body_3_index_exch_312'
 if(read_three_body_ints)then
  print*,'Reading three_body_ints from disk ...'
  call read_array_3_index_tensor(mo_num,three_body_3_index_exch_312,name_file)
 else
  provide x_W_ij_erf_rk
  three_body_3_index_exch_312 = 0.d0
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,integral) & 
  !$OMP SHARED (mo_num,three_body_3_index_exch_312)
  !$OMP DO SCHEDULE (guided) COLLAPSE(3)
   do m = 1, mo_num ! 3
    do j = 1, mo_num ! 2 
     do i = 1, mo_num ! 1 
      integral = 0.d0
      !                          1 2 3 1 2 3
      call give_integrals_3_body(i,j,m,m,i,j,integral)
 
      three_body_3_index_exch_312(i,j,m) = -1.d0 * integral 
   
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_3_index_312',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_3_index_exch_312 on disk ...' 
  call write_array_3_index_tensor(mo_num,three_body_3_index_exch_312,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 

subroutine write_array_3_index_tensor(n_orb,array_tmp,name_file)
 implicit none
 integer, intent(in) :: n_orb
 character*(128),  intent(in) :: name_file 
 double precision, intent(in) :: array_tmp(n_orb,n_orb,n_orb)

 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 PROVIDE ezfio_filename                                                                                                  
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'W')
 write(i_unit_output)array_tmp
 close(unit=i_unit_output)
end

subroutine read_array_3_index_tensor(n_orb,array_tmp,name_file)
 implicit none
 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 integer, intent(in) :: n_orb
 character*(128),  intent(in)  :: name_file 
 double precision, intent(out) :: array_tmp(n_orb,n_orb,n_orb)
 PROVIDE ezfio_filename                                                                                                  
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'R')
 read(i_unit_output)array_tmp
 close(unit=i_unit_output)
end
