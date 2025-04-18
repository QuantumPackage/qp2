BEGIN_PROVIDER [ integer*8, qmckl_ctx_jastrow ]
  use qmckl
  use iso_c_binding
  implicit none
  BEGIN_DOC
  ! Context for the QMCKL library
  END_DOC
  integer(qmckl_exit_code) :: rc

  qmckl_ctx_jastrow = qmckl_context_create()

  if (.not.jast_qmckl_spin_independent) then
    print *, 'WARNING: In QMCkl Jastrow, jast_qmckl_spin_independent should to be set to True'
  endif
  rc =  qmckl_set_jastrow_champ_spin_independent(qmckl_ctx_jastrow, jast_qmckl_spin_independent)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc = qmckl_set_nucleus_num(qmckl_ctx_jastrow, nucl_num*1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc = qmckl_set_nucleus_charge(qmckl_ctx_jastrow, nucl_charge, nucl_num*1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc = qmckl_set_nucleus_coord(qmckl_ctx_jastrow, 'T', nucl_coord, nucl_num*3_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc = qmckl_set_electron_num(qmckl_ctx_jastrow, 1_8, 1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1


  ! Jastrow parameters
  rc =  qmckl_set_jastrow_champ_type_nucl_num(qmckl_ctx_jastrow, 1_8*jast_qmckl_type_nucl_num)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_type_nucl_vector(qmckl_ctx_jastrow, 1_8*jast_qmckl_type_nucl_vector-1_8, 1_8*nucl_num)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_rescale_factor_ee(qmckl_ctx_jastrow, jast_qmckl_rescale_ee)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_rescale_factor_en(qmckl_ctx_jastrow, jast_qmckl_rescale_en, 1_8*jast_qmckl_type_nucl_num)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_aord_num(qmckl_ctx_jastrow, jast_qmckl_aord_num*1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_a_vector(qmckl_ctx_jastrow, jast_qmckl_a_vector, 1_8*size(jast_qmckl_a_vector))
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_bord_num(qmckl_ctx_jastrow, jast_qmckl_bord_num*1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  rc =  qmckl_set_jastrow_champ_b_vector(qmckl_ctx_jastrow, jast_qmckl_b_vector, 1_8*size(jast_qmckl_b_vector))
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1


  rc =  qmckl_set_jastrow_champ_cord_num(qmckl_ctx_jastrow, jast_qmckl_cord_num*1_8)
  rc = qmckl_check(qmckl_ctx_jastrow, rc)
  if (rc /= QMCKL_SUCCESS) stop -1

  if (jast_qmckl_cord_num > 0) then
    rc =  qmckl_set_jastrow_champ_c_vector(qmckl_ctx_jastrow, jast_qmckl_c_vector, 1_8*jast_qmckl_c_vector_size)
    rc = qmckl_check(qmckl_ctx_jastrow, rc)
    if (rc /= QMCKL_SUCCESS) stop -1
  endif
!  print*,'jast_qmckl_cord_num = ',jast_qmckl_cord_num
!  integer :: i
!  do i = 1, jast_qmckl_c_vector_size
!   print*,jast_qmckl_c_vector(i) 
!  enddo

END_PROVIDER


 BEGIN_PROVIDER [ double precision, aos_in_r_array_qmckl, (ao_num,n_points_final_grid)]
&BEGIN_PROVIDER [ double precision, aos_grad_in_r_array_qmckl, (ao_num,n_points_final_grid,3)]
&BEGIN_PROVIDER [ double precision, aos_lapl_in_r_array_qmckl, (ao_num, n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! AOS computed with qmckl
 END_DOC
 use qmckl

 integer*8 :: qmckl_ctx
 integer(qmckl_exit_code) :: rc

 qmckl_ctx = qmckl_context_create()

 rc = qmckl_trexio_read(qmckl_ctx, trexio_file, 1_8*len(trim(trexio_filename)))
 if (rc /= QMCKL_SUCCESS) then
   print *, irp_here, 'qmckl error in read_trexio'
   rc = qmckl_check(qmckl_ctx, rc)
   stop -1
 endif

 rc = qmckl_set_point(qmckl_ctx, 'N', n_points_final_grid*1_8, final_grid_points, n_points_final_grid*3_8)
 if (rc /= QMCKL_SUCCESS) then
   print *, irp_here, 'qmckl error in set_electron_point'
   rc = qmckl_check(qmckl_ctx, rc)
   stop -1
 endif

 double precision, allocatable :: vgl(:,:,:)
 allocate( vgl(ao_num,5,n_points_final_grid))
 rc = qmckl_get_ao_basis_ao_vgl_inplace(qmckl_ctx, vgl, n_points_final_grid*ao_num*5_8)
 if (rc /= QMCKL_SUCCESS) then
   print *, irp_here, 'qmckl error in get_ao_vgl'
   rc = qmckl_check(qmckl_ctx, rc)
   stop -1
 endif

 integer :: i,k
 do k=1,n_points_final_grid
   do i=1,ao_num
     aos_in_r_array_qmckl(i,k) = vgl(i,1,k)
     aos_grad_in_r_array_qmckl(i,k,1) = vgl(i,2,k)
     aos_grad_in_r_array_qmckl(i,k,2) = vgl(i,3,k)
     aos_grad_in_r_array_qmckl(i,k,3) = vgl(i,4,k)
     aos_lapl_in_r_array_qmckl(i,k) = vgl(i,5,k)
   enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [ double precision, mos_in_r_array_qmckl, (mo_num,n_points_final_grid)]
&BEGIN_PROVIDER [ double precision, mos_grad_in_r_array_qmckl, (mo_num,n_points_final_grid,3)]
&BEGIN_PROVIDER [ double precision, mos_lapl_in_r_array_qmckl, (mo_num, n_points_final_grid)]
 implicit none
 BEGIN_DOC
 ! moS computed with qmckl
 END_DOC
 use qmckl

 integer*8 :: qmckl_ctx
 integer(qmckl_exit_code) :: rc

 qmckl_ctx = qmckl_context_create()

 rc = qmckl_trexio_read(qmckl_ctx, trexio_file, 1_8*len(trim(trexio_filename)))
 if (rc /= QMCKL_SUCCESS) then
   print *, irp_here, 'qmckl error in read_trexio'
   rc = qmckl_check(qmckl_ctx, rc)
   stop -1
 endif

 rc = qmckl_set_point(qmckl_ctx, 'N', n_points_final_grid*1_8, final_grid_points, n_points_final_grid*3_8)
 if (rc /= QMCKL_SUCCESS) then
   print *, irp_here, 'qmckl error in set_electron_point'
   rc = qmckl_check(qmckl_ctx, rc)
   stop -1
 endif

 double precision, allocatable :: vgl(:,:,:)
 allocate( vgl(mo_num,5,n_points_final_grid))
 rc = qmckl_get_mo_basis_mo_vgl(qmckl_ctx, vgl, n_points_final_grid*mo_num*5_8)
 if (rc /= QMCKL_SUCCESS) then
   print *, irp_here, 'qmckl error in get_mo_vgl'
   rc = qmckl_check(qmckl_ctx, rc)
   stop -1
 endif

 integer :: i,k
 do k=1,n_points_final_grid
   do i=1,mo_num
     mos_in_r_array_qmckl(i,k) = vgl(i,1,k)
     mos_grad_in_r_array_qmckl(i,k,1) = vgl(i,2,k)
     mos_grad_in_r_array_qmckl(i,k,2) = vgl(i,3,k)
     mos_grad_in_r_array_qmckl(i,k,3) = vgl(i,4,k)
     mos_lapl_in_r_array_qmckl(i,k) = vgl(i,5,k)
   enddo
 enddo

END_PROVIDER


