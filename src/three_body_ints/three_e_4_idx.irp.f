
BEGIN_PROVIDER [ double precision, three_body_4_index, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 4 index matrix direct element of the -L  three-body operator 
!
! three_body_4_index(j,m,k,i) = < phi_j phi_m phi_k | phi_j phi_m phi_i >
!
! notice the -1 sign: in this way three_body_4_index can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 three_body_4_index = 0.d0
 print*,'Providing the three_body_4_index ...'
 call wall_time(wall0)

 name_file = 'three_body_4_index'
 if(read_three_body_ints)then
  print*,'Reading three_body_4_index from disk ...'
  call read_array_4_index_tensor(mo_num,three_body_4_index,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,k,integral) & 
  !$OMP SHARED (mo_num,three_body_4_index)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do i = 1, mo_num
    do k = 1, mo_num
     do m = 1, mo_num
      do j = 1, mo_num
       integral = 0.d0
       call give_integrals_3_body(i,j,m,k,j,m,integral)
 
       three_body_4_index(j,m,k,i) = -1.d0 * integral 
   
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_4_index',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_4_index on disk ...' 
  call write_array_4_index_tensor(mo_num,three_body_4_index,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_4_index_exch_12, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 4 index matrix EXCHANGE element of the -L  three-body operator 
!                                         
! three_body_4_index_exch_12(j,m,k,i) = < phi_m phi_j phi_i | phi_j phi_m phi_k >
!
! notice the -1 sign: in this way three_body_3_index_exch_12 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 three_body_4_index_exch_12 = 0.d0
 print*,'Providing the three_body_4_index_exch_12 ...'
 call wall_time(wall0)

 name_file = 'three_body_4_index_exch_12'
 if(read_three_body_ints)then
  print*,'Reading three_body_4_index_exch_12 from disk ...'
  call read_array_4_index_tensor(mo_num,three_body_4_index_exch_12,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,k,integral) & 
  !$OMP SHARED (mo_num,three_body_4_index_exch_12)
  !$OMP DO SCHEDULE (guided) COLLAPSE(4)
   do i = 1, mo_num
    do k = 1, mo_num
     do m = 1, mo_num
      do j = 1, mo_num
       integral = 0.d0
       call give_integrals_3_body(i,m,j,k,j,m,integral)
 
       three_body_4_index_exch_12(j,m,k,i) = -1.d0 * integral 
   
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_4_index_exch_12',wall1 - wall0

 if(write_three_body_ints)then
  print*,'Writing three_body_4_index_exch_12 on disk ...' 
  call write_array_4_index_tensor(mo_num,three_body_4_index_exch_12,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_4_index_exch_12_part, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 4 index matrix EXCHANGE element of the -L  three-body operator 
!
! three_body_4_index_exch_12_part(j,m,k,i) = < phi_m phi_j phi_i | phi_m phi_k phi_j >
!
! notice the -1 sign: in this way three_body_3_index_exch_12 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 three_body_4_index_exch_12_part = 0.d0
 print*,'Providing the three_body_4_index_exch_12_part ...'
 call wall_time(wall0)

 name_file = 'three_body_4_index_exch_12_part'
 if(read_three_body_ints)then
  print*,'Reading three_body_4_index_exch_12_part from disk ...'
  call read_array_4_index_tensor(mo_num,three_body_4_index_exch_12_part,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,k,integral) & 
  !$OMP SHARED (mo_num,three_body_4_index_exch_12_part)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do i = 1, mo_num
    do k = 1, mo_num
     do m = 1, mo_num
      do j = 1, mo_num
       integral = 0.d0
       !                          
       call give_integrals_3_body(i,j,m,j,k,m,integral)
       three_body_4_index_exch_12_part(j,m,k,i) = -1.d0 * integral 
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
  call wall_time(wall1)
 endif
 print*,'wall time for three_body_4_index_exch_12_part',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_4_index_exch_12_part on disk ...' 
  call write_array_4_index_tensor(mo_num,three_body_4_index_exch_12_part,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_4_index_exch_12_part_bis, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 4 index matrix EXCHANGE element of the -L  three-body operator 
!
! three_body_4_index_exch_12_part_bis(j,m,k,i) = < phi_m phi_j phi_i | phi_m phi_k phi_j >
!
! notice the -1 sign: in this way three_body_3_index_exch_12 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 three_body_4_index_exch_12_part_bis = 0.d0
 print*,'Providing the three_body_4_index_exch_12_part_bis ...'
 call wall_time(wall0)

 name_file = 'three_body_4_index_exch_12_part_bis'
 if(read_three_body_ints)then
  print*,'Reading three_body_4_index_exch_12_part_bisfrom disk ...'
  call read_array_4_index_tensor(mo_num,three_body_4_index_exch_12_part_bis,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,k,integral) & 
  !$OMP SHARED (mo_num,three_body_4_index_exch_12_part_bis)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do i = 1, mo_num
    do k = 1, mo_num
     do m = 1, mo_num
      do j = 1, mo_num
       integral = 0.d0
       !                          
       call give_integrals_3_body(i,j,m,m,j,k,integral)
 
       three_body_4_index_exch_12_part_bis(j,m,k,i) = -1.d0 * integral 
   
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_4_index_exch_12_part_bis',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_4_index_exch_12_part_bis on disk ...' 
  call write_array_4_index_tensor(mo_num,three_body_4_index_exch_12_part_bis,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 


BEGIN_PROVIDER [ double precision, three_body_4_index_exch_231, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 4 index matrix direct element of the -L  three-body operator 
!
! three_body_4_index_exch_231(j,m,k,i) = < phi_j phi_m phi_k | phi_j phi_m phi_i >
!
! notice the -1 sign: in this way three_body_4_index_exch_231 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 three_body_4_index_exch_231 = 0.d0
 print*,'Providing the three_body_4_index_exch_231 ...'
 call wall_time(wall0)
 name_file = 'three_body_4_index_exch_231'
 if(read_three_body_ints)then
  print*,'Reading three_body_4_index_exch_231 from disk ...'
  call read_array_4_index_tensor(mo_num,three_body_4_index_exch_231,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,k,integral) & 
  !$OMP SHARED (mo_num,three_body_4_index_exch_231)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do i = 1, mo_num
    do k = 1, mo_num
     do m = 1, mo_num
      do j = 1, mo_num
       integral = 0.d0
       call give_integrals_3_body(i,j,m,j,m,k,integral)
 
       three_body_4_index_exch_231(j,m,k,i) = -1.d0 * integral 
   
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_4_index_exch_231',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_4_index_exch_231 on disk ...' 
  call write_array_4_index_tensor(mo_num,three_body_4_index_exch_231,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_body_4_index_exch_312, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! 4 index matrix direct element of the -L  three-body operator 
!
! three_body_4_index_exch_312(j,m,k,i) = < phi_j phi_m phi_k | phi_j phi_m phi_i >
!
! notice the -1 sign: in this way three_body_4_index_exch_312 can be directly used to compute Slater rules :)
 END_DOC
 integer :: i,j,k,l,m,n
 double precision :: integral, wall1, wall0
 character*(128) :: name_file 
 three_body_4_index_exch_312 = 0.d0
 print*,'Providing the three_body_4_index_exch_312 ...'
 call wall_time(wall0)
 name_file = 'three_body_4_index_exch_312'
 if(read_three_body_ints)then
  print*,'Reading three_body_4_index_exch_312 from disk ...'
  call read_array_4_index_tensor(mo_num,three_body_4_index_exch_312,name_file)
 else
  provide x_W_ij_erf_rk
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,m,k,integral) & 
  !$OMP SHARED (mo_num,three_body_4_index_exch_312)
  !$OMP DO SCHEDULE (guided) COLLAPSE(2)
   do i = 1, mo_num
    do k = 1, mo_num
     do m = 1, mo_num
      do j = 1, mo_num
       integral = 0.d0
       call give_integrals_3_body(i,j,m,m,k,j,integral)
 
       three_body_4_index_exch_312(j,m,k,i) = -1.d0 * integral 
   
      enddo
     enddo
    enddo
   enddo
  !$OMP END DO
  !$OMP END PARALLEL
 endif
 call wall_time(wall1)
 print*,'wall time for three_body_4_index_exch_312',wall1 - wall0
 if(write_three_body_ints)then
  print*,'Writing three_body_4_index_exch_312 on disk ...' 
  call write_array_4_index_tensor(mo_num,three_body_4_index_exch_312,name_file)
  call ezfio_set_three_body_ints_io_three_body_ints("Read")
 endif

END_PROVIDER 

subroutine write_array_4_index_tensor(n_orb,array_tmp,name_file)
 implicit none
 integer, intent(in) :: n_orb
 character*(128),  intent(in) :: name_file 
 double precision, intent(in) :: array_tmp(n_orb,n_orb,n_orb,n_orb)

 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 PROVIDE ezfio_filename                                                                                                  
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'W')
 write(i_unit_output)array_tmp
 close(unit=i_unit_output)
end

subroutine read_array_4_index_tensor(n_orb,array_tmp,name_file)
 implicit none
 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 integer, intent(in) :: n_orb
 character*(128),  intent(in)  :: name_file 
 double precision, intent(out) :: array_tmp(n_orb,n_orb,n_orb,n_orb)
 PROVIDE ezfio_filename                                                                                                  
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'R')
 read(i_unit_output)array_tmp
 close(unit=i_unit_output)
end
