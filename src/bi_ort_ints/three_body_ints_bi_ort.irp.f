
! ---

BEGIN_PROVIDER [ double precision, three_body_ints_bi_ort, (mo_num, mo_num, mo_num, mo_num, mo_num, mo_num)]

 BEGIN_DOC
! matrix element of the -L  three-body operator 
!
! notice the -1 sign: in this way three_body_ints_bi_ort can be directly used to compute Slater rules :)
 END_DOC

 implicit none
 integer          :: i, j, k, l, m, n
 double precision :: integral, wall1, wall0
 character*(128)  :: name_file 

  three_body_ints_bi_ort = 0.d0
  print *, ' Providing the three_body_ints_bi_ort ...'
  call wall_time(wall0)
  name_file = 'six_index_tensor'

! if(read_three_body_ints_bi_ort)then
!  call read_fcidump_3_tc(three_body_ints_bi_ort)
! else
!  if(read_three_body_ints_bi_ort)then
!   print*,'Reading three_body_ints_bi_ort from disk ...'
!   call read_array_6_index_tensor(mo_num,three_body_ints_bi_ort,name_file)
!  else

  !provide x_W_ki_bi_ortho_erf_rk 
  provide mos_r_in_r_array_transp mos_l_in_r_array_transp

 !$OMP PARALLEL                       &
 !$OMP DEFAULT (NONE)                 &
 !$OMP PRIVATE (i,j,k,l,m,n,integral) & 
 !$OMP SHARED (mo_num,three_body_ints_bi_ort)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do j = 1, mo_num
      do m = 1, mo_num
        do k = 1, mo_num
          do l = 1, mo_num
            do n = 1, mo_num
              call give_integrals_3_body_bi_ort(n, l, k, m, j, i, integral)

              three_body_ints_bi_ort(n,l,k,m,j,i) = -1.d0 * integral 
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL
!  endif
! endif

  call wall_time(wall1)
  print *, ' wall time for three_body_ints_bi_ort', wall1 - wall0
! if(write_three_body_ints_bi_ort)then
!  print*,'Writing three_body_ints_bi_ort on disk ...'
!  call write_array_6_index_tensor(mo_num,three_body_ints_bi_ort,name_file)
!  call ezfio_set_three_body_ints_bi_ort_io_three_body_ints_bi_ort("Read")
! endif

END_PROVIDER 

! ---

subroutine give_integrals_3_body_bi_ort(n, l, k, m, j, i, integral)

  BEGIN_DOC
  !
  ! < n l k | -L | m j i > with a BI-ORTHONORMAL MOLECULAR ORBITALS 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n, l, k, m, j, i
  double precision, intent(out) :: integral
  integer                       :: ipoint
  double precision              :: weight

  integral = 0.d0
  do ipoint = 1, n_points_final_grid
    weight = final_weight_at_r_vector(ipoint)                                                                          

    integral += weight * mos_l_in_r_array_transp(ipoint,k) * mos_r_in_r_array_transp(ipoint,i)        & 
              * ( int2_grad1_u12_bimo_t(ipoint,1,n,m) * int2_grad1_u12_bimo_t(ipoint,1,l,j) &
                + int2_grad1_u12_bimo_t(ipoint,2,n,m) * int2_grad1_u12_bimo_t(ipoint,2,l,j) &
                + int2_grad1_u12_bimo_t(ipoint,3,n,m) * int2_grad1_u12_bimo_t(ipoint,3,l,j) )
    integral += weight * mos_l_in_r_array_transp(ipoint,l) * mos_r_in_r_array_transp(ipoint,j)        & 
              * ( int2_grad1_u12_bimo_t(ipoint,1,n,m) * int2_grad1_u12_bimo_t(ipoint,1,k,i) &
                + int2_grad1_u12_bimo_t(ipoint,2,n,m) * int2_grad1_u12_bimo_t(ipoint,2,k,i) &
                + int2_grad1_u12_bimo_t(ipoint,3,n,m) * int2_grad1_u12_bimo_t(ipoint,3,k,i) )
    integral += weight * mos_l_in_r_array_transp(ipoint,n) * mos_r_in_r_array_transp(ipoint,m)        &
              * ( int2_grad1_u12_bimo_t(ipoint,1,l,j) * int2_grad1_u12_bimo_t(ipoint,1,k,i) &
                + int2_grad1_u12_bimo_t(ipoint,2,l,j) * int2_grad1_u12_bimo_t(ipoint,2,k,i) &
                + int2_grad1_u12_bimo_t(ipoint,3,l,j) * int2_grad1_u12_bimo_t(ipoint,3,k,i) )

  enddo

end subroutine give_integrals_3_body_bi_ort

! ---

subroutine give_integrals_3_body_bi_ort_old(n, l, k, m, j, i, integral)

  BEGIN_DOC
  !
  ! < n l k | -L | m j i > with a BI-ORTHONORMAL MOLECULAR ORBITALS 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n, l, k, m, j, i
  double precision, intent(out) :: integral
  integer                       :: ipoint
  double precision              :: weight

  integral = 0.d0
  do ipoint = 1, n_points_final_grid
    weight = final_weight_at_r_vector(ipoint)                                                                          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    integral += weight * mos_l_in_r_array_transp(ipoint,k) * mos_r_in_r_array_transp(ipoint,i) & 
!              * ( x_W_ki_bi_ortho_erf_rk(ipoint,1,n,m) * x_W_ki_bi_ortho_erf_rk(ipoint,1,l,j)  &
!                + x_W_ki_bi_ortho_erf_rk(ipoint,2,n,m) * x_W_ki_bi_ortho_erf_rk(ipoint,2,l,j)  &
!                + x_W_ki_bi_ortho_erf_rk(ipoint,3,n,m) * x_W_ki_bi_ortho_erf_rk(ipoint,3,l,j)  )
!    integral += weight * mos_l_in_r_array_transp(ipoint,l) * mos_r_in_r_array_transp(ipoint,j) & 
!              * ( x_W_ki_bi_ortho_erf_rk(ipoint,1,n,m) * x_W_ki_bi_ortho_erf_rk(ipoint,1,k,i)  &
!                + x_W_ki_bi_ortho_erf_rk(ipoint,2,n,m) * x_W_ki_bi_ortho_erf_rk(ipoint,2,k,i)  &
!                + x_W_ki_bi_ortho_erf_rk(ipoint,3,n,m) * x_W_ki_bi_ortho_erf_rk(ipoint,3,k,i)  )
!    integral += weight * mos_l_in_r_array_transp(ipoint,n) * mos_r_in_r_array_transp(ipoint,m) &
!              * ( x_W_ki_bi_ortho_erf_rk(ipoint,1,l,j) * x_W_ki_bi_ortho_erf_rk(ipoint,1,k,i)  &
!                + x_W_ki_bi_ortho_erf_rk(ipoint,2,l,j) * x_W_ki_bi_ortho_erf_rk(ipoint,2,k,i)  &
!                + x_W_ki_bi_ortho_erf_rk(ipoint,3,l,j) * x_W_ki_bi_ortho_erf_rk(ipoint,3,k,i)  )

!    integral += weight * mos_l_in_r_array_transp(ipoint,k) * mos_r_in_r_array_transp(ipoint,i) & 
!              * ( int2_grad1_u12_bimo(1,n,m,ipoint) * int2_grad1_u12_bimo(1,l,j,ipoint)        &
!                + int2_grad1_u12_bimo(2,n,m,ipoint) * int2_grad1_u12_bimo(2,l,j,ipoint)        &
!                + int2_grad1_u12_bimo(3,n,m,ipoint) * int2_grad1_u12_bimo(3,l,j,ipoint)        )
!    integral += weight * mos_l_in_r_array_transp(ipoint,l) * mos_r_in_r_array_transp(ipoint,j) & 
!              * ( int2_grad1_u12_bimo(1,n,m,ipoint) * int2_grad1_u12_bimo(1,k,i,ipoint)        &
!                + int2_grad1_u12_bimo(2,n,m,ipoint) * int2_grad1_u12_bimo(2,k,i,ipoint)        &
!                + int2_grad1_u12_bimo(3,n,m,ipoint) * int2_grad1_u12_bimo(3,k,i,ipoint)        )
!    integral += weight * mos_l_in_r_array_transp(ipoint,n) * mos_r_in_r_array_transp(ipoint,m) &
!              * ( int2_grad1_u12_bimo(1,l,j,ipoint) * int2_grad1_u12_bimo(1,k,i,ipoint)        &
!                + int2_grad1_u12_bimo(2,l,j,ipoint) * int2_grad1_u12_bimo(2,k,i,ipoint)        &
!                + int2_grad1_u12_bimo(3,l,j,ipoint) * int2_grad1_u12_bimo(3,k,i,ipoint)        )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    integral += weight * mos_l_in_r_array_transp(ipoint,k) * mos_r_in_r_array_transp(ipoint,i)        & 
              * ( int2_grad1_u12_bimo_transp(n,m,1,ipoint) * int2_grad1_u12_bimo_transp(l,j,1,ipoint) &
                + int2_grad1_u12_bimo_transp(n,m,2,ipoint) * int2_grad1_u12_bimo_transp(l,j,2,ipoint) &
                + int2_grad1_u12_bimo_transp(n,m,3,ipoint) * int2_grad1_u12_bimo_transp(l,j,3,ipoint) )
    integral += weight * mos_l_in_r_array_transp(ipoint,l) * mos_r_in_r_array_transp(ipoint,j)        & 
              * ( int2_grad1_u12_bimo_transp(n,m,1,ipoint) * int2_grad1_u12_bimo_transp(k,i,1,ipoint) &
                + int2_grad1_u12_bimo_transp(n,m,2,ipoint) * int2_grad1_u12_bimo_transp(k,i,2,ipoint) &
                + int2_grad1_u12_bimo_transp(n,m,3,ipoint) * int2_grad1_u12_bimo_transp(k,i,3,ipoint) )
    integral += weight * mos_l_in_r_array_transp(ipoint,n) * mos_r_in_r_array_transp(ipoint,m)        &
              * ( int2_grad1_u12_bimo_transp(l,j,1,ipoint) * int2_grad1_u12_bimo_transp(k,i,1,ipoint) &
                + int2_grad1_u12_bimo_transp(l,j,2,ipoint) * int2_grad1_u12_bimo_transp(k,i,2,ipoint) &
                + int2_grad1_u12_bimo_transp(l,j,3,ipoint) * int2_grad1_u12_bimo_transp(k,i,3,ipoint) )

  enddo

end subroutine give_integrals_3_body_bi_ort_old

! ---

subroutine give_integrals_3_body_bi_ort_ao(n, l, k, m, j, i, integral)

  BEGIN_DOC
  !
  ! < n l k | -L | m j i > with a BI-ORTHONORMAL ATOMIC ORBITALS 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n, l, k, m, j, i
  double precision, intent(out) :: integral
  integer                       :: ipoint
  double precision              :: weight

  integral = 0.d0
  do ipoint = 1, n_points_final_grid
    weight = final_weight_at_r_vector(ipoint)                                                                          

    integral += weight * aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) & 
              * ( int2_grad1_u12_ao_t(ipoint,1,n,m) * int2_grad1_u12_ao_t(ipoint,1,l,j)    &
                + int2_grad1_u12_ao_t(ipoint,2,n,m) * int2_grad1_u12_ao_t(ipoint,2,l,j)    &
                + int2_grad1_u12_ao_t(ipoint,3,n,m) * int2_grad1_u12_ao_t(ipoint,3,l,j) )
    integral += weight * aos_in_r_array_transp(ipoint,l) * aos_in_r_array_transp(ipoint,j) & 
              * ( int2_grad1_u12_ao_t(ipoint,1,n,m) * int2_grad1_u12_ao_t(ipoint,1,k,i)    &
                + int2_grad1_u12_ao_t(ipoint,2,n,m) * int2_grad1_u12_ao_t(ipoint,2,k,i)    &
                + int2_grad1_u12_ao_t(ipoint,3,n,m) * int2_grad1_u12_ao_t(ipoint,3,k,i) )
    integral += weight * aos_in_r_array_transp(ipoint,n) * aos_in_r_array_transp(ipoint,m) &
              * ( int2_grad1_u12_ao_t(ipoint,1,l,j) * int2_grad1_u12_ao_t(ipoint,1,k,i)    &
                + int2_grad1_u12_ao_t(ipoint,2,l,j) * int2_grad1_u12_ao_t(ipoint,2,k,i)    &
                + int2_grad1_u12_ao_t(ipoint,3,l,j) * int2_grad1_u12_ao_t(ipoint,3,k,i) )

  enddo

end subroutine give_integrals_3_body_bi_ort_ao

! ---
