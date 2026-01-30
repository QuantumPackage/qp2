use gpu

 BEGIN_PROVIDER [ integer, cholesky_mo_num ]
&BEGIN_PROVIDER [ integer, cholesky_mo_num_split, (1:5)]
 implicit none
 BEGIN_DOC
 ! Number of Cholesky vectors in MO basis
 END_DOC
 integer, external              :: getUnitAndOpen
 integer                        :: iunit
 if (read_mo_cholesky) then
   iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_mo_transp', 'R')
   read(iunit) cholesky_mo_num
   close(iunit)
 else
   cholesky_mo_num = cholesky_ao_num
 endif
 cholesky_mo_num_split(1) = 0
 cholesky_mo_num_split(2) = cholesky_mo_num/4
 cholesky_mo_num_split(3) = 2*cholesky_mo_num_split(2)
 cholesky_mo_num_split(4) = 3*cholesky_mo_num_split(2)
 cholesky_mo_num_split(5) = cholesky_mo_num
 cholesky_mo_num_split += 1
END_PROVIDER

BEGIN_PROVIDER [ double precision, cholesky_mo, (mo_num, mo_num, cholesky_mo_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in MO basis
 END_DOC

 integer :: k, i, j

 !$OMP PARALLEL DO PRIVATE(k)
 do k=1,cholesky_mo_num
  do j=1,mo_num
    do i=1,mo_num
      cholesky_mo(i,j,k) = cholesky_mo_transp(k,i,j)
    enddo
  enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER


 BEGIN_PROVIDER [ double precision, cholesky_mo_transp, (cholesky_mo_num, mo_num, mo_num) ]
&BEGIN_PROVIDER [ type(gpu_double3), cholesky_mo_transp_d, (0:gpu_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in MO basis. Warning: it is transposed wrt cholesky_ao:
 !
 ! - cholesky_ao        is (ao_num^2 x cholesky_ao_num)
 !
 ! - cholesky_mo_transp is (cholesky_mo_num x mo_num^2)
 END_DOC

 type(gpu_double3) :: X_d, cholesky_ao_d
 type(gpu_double2) :: mo_coef_d
 double precision :: wall0, wall1
 integer, external              :: getUnitAndOpen
 integer                        :: iunit, ierr, rank
 integer :: igpu

 if (read_mo_cholesky) then
      print *,  'Reading Cholesky MO vectors from disk...'
     iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_mo_transp', 'R')
     read(iunit) rank
     if (cholesky_mo_num /= rank) then
        stop 'inconsistent rank'
     endif
     read(iunit) cholesky_mo_transp
     !$OMP PARALLEL DO PRIVATE(igpu)
     do igpu=0,gpu_num-1
       call gpu_allocate(cholesky_mo_transp_d(igpu), cholesky_mo_num, mo_num, mo_num)
       call gpu_set_device(igpu)
       call gpu_allocate(cholesky_mo_transp_d(igpu), cholesky_mo_num, mo_num, mo_num)
       call gpu_upload(cholesky_mo_transp, cholesky_mo_transp_d(igpu))
     enddo
     !$OMP END PARALLEL DO
     call gpu_set_device(0)
     close(iunit)
 else
     print *, ''
     print *, 'AO->MO Transformation of Cholesky vectors'
     print *, '-----------------------------------------'
     print *, ''

     call wall_time(wall0)

     if (size(cholesky_ao,3) /= cholesky_ao_num) then
       call qp_bug(irp_here, size(cholesky_ao,3), 'size(cholesky_ao,3) /= cholesky_ao_num')
     endif

     if (gpu_num > 0) then
       call gpu_allocate(cholesky_mo_transp_d(0), cholesky_mo_num, mo_num, mo_num)

       call gpu_allocate(mo_coef_d, ao_num,mo_num)
       call gpu_upload(mo_coef, mo_coef_d)

       call gpu_allocate(cholesky_ao_d, ao_num,ao_num,cholesky_ao_num)
       call gpu_upload(cholesky_ao, cholesky_ao_d)

       call gpu_allocate(X_d, mo_num,cholesky_mo_num,ao_num)

       call gpu_dgemm(blas_handle, 'T','N', ao_num*cholesky_mo_num, mo_num, ao_num, 1.d0, &
           cholesky_ao_d%f(1,1,1), ao_num, mo_coef_d%f(1,1), ao_num, 0.d0, X_d%f(1,1,1), ao_num*cholesky_mo_num)

       call gpu_dgemm(blas_handle, 'T','N', cholesky_mo_num*mo_num, mo_num, ao_num, 1.d0, &
           X_d%f(1,1,1), ao_num, mo_coef_d%f(1,1), ao_num, 0.d0, cholesky_mo_transp_d(0)%f(1,1,1), cholesky_mo_num*mo_num)

       call gpu_synchronize()
       call gpu_deallocate(cholesky_ao_d)
       call gpu_deallocate(X_d)
       call gpu_deallocate(mo_coef_d)

       call gpu_download(cholesky_mo_transp_d(0),cholesky_mo_transp)

       !$OMP PARALLEL DO PRIVATE(igpu)
       do igpu=1,gpu_num-1
         call gpu_set_device(igpu)
         call gpu_allocate(cholesky_mo_transp_d(igpu), cholesky_mo_num, mo_num, mo_num)
         call gpu_upload(cholesky_mo_transp, cholesky_mo_transp_d(igpu))
       enddo
       !$OMP END PARALLEL DO
       call gpu_set_device(0)
     else
       double precision, allocatable :: X(:,:,:)
       allocate(X(mo_num,cholesky_mo_num,ao_num))

       call dgemm('T','N', ao_num*cholesky_mo_num, mo_num, ao_num, 1.d0, &
           cholesky_ao(1,1,1), ao_num, mo_coef(1,1), ao_num, 0.d0, X(1,1,1), ao_num*cholesky_mo_num)

       call dgemm('T','N', cholesky_mo_num*mo_num, mo_num, ao_num, 1.d0, &
           X(1,1,1), ao_num, mo_coef(1,1), ao_num, 0.d0, cholesky_mo_transp(1,1,1), cholesky_mo_num*mo_num)

       deallocate(X)
     endif

     call wall_time(wall1)
     print*,'Time to provide MO cholesky vectors = ',(wall1-wall0)/60.d0, ' min'


     if (write_mo_cholesky) then
       print *,  'Writing Cholesky MO vectors to disk...'
       iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_mo_transp', 'W')
       write(iunit) cholesky_mo_num
       write(iunit) cholesky_mo_transp
       close(iunit)
       call ezfio_set_mo_two_e_ints_io_mo_cholesky('Read')
     endif
  endif

END_PROVIDER



BEGIN_PROVIDER [ double precision, cholesky_semi_mo_transp_simple, (cholesky_mo_num, ao_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in MO basis
 END_DOC

 double precision, allocatable :: X(:,:,:)
 double precision :: wall0, wall1
 integer :: ierr
 print *, 'Semi AO->MO Transformation of Cholesky vectors'
  call wall_time(wall0)

 allocate(X(mo_num,cholesky_mo_num,ao_num), stat=ierr)
 if (ierr /= 0) then
   print *, irp_here, ': Allocation failed'
 endif
 integer :: i_chol, i_mo, j_mo, i_ao 
 cholesky_semi_mo_transp_simple = 0.d0
 do i_mo = 1, mo_num
  do i_ao = 1, ao_num
   do j_mo = 1, mo_num
    do i_chol = 1, cholesky_mo_num
     cholesky_semi_mo_transp_simple(i_chol, i_ao,i_mo) += cholesky_mo_transp(i_chol,j_mo,i_mo) * mo_coef_transp(j_mo,i_ao)
    enddo
   enddo
  enddo
 enddo

END_PROVIDER




BEGIN_PROVIDER [ real, cholesky_mo_sp, (mo_num, mo_num, cholesky_mo_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in MO basis in stored in single precision
 END_DOC

 integer :: k, i, j

 !$OMP PARALLEL DO PRIVATE(k)
 do k=1,cholesky_mo_num
  do j=1,mo_num
    do i=1,mo_num
      cholesky_mo_sp(i,j,k) = cholesky_mo_transp_sp(k,i,j)
    enddo
  enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER

BEGIN_PROVIDER [ real, cholesky_mo_transp_sp, (cholesky_mo_num, mo_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in MO basis in s. Warning: it is transposed wrt cholesky_ao:
 !
 ! - cholesky_ao        is (ao_num^2 x cholesky_ao_num)
 !
 ! - cholesky_mo_transp is (cholesky_mo_num x mo_num^2)
 END_DOC

 integer :: i,j,k
 !$OMP PARALLEL DO PRIVATE(k)
 do j=1,mo_num
  do i=1,mo_num
   do k=1,cholesky_mo_num
      cholesky_mo_transp_sp(k,i,j) = real(cholesky_mo_transp(k,i,j),4)
    enddo
  enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER

BEGIN_PROVIDER [ type(gpu_real3), cholesky_mo_transp_sp_d, (0:gpu_num) ]
 BEGIN_DOC
 ! Cholesky vectors in MO basis, on GPU, single precision.
 END_DOC
 integer :: igpu
 if (gpu_num > 0) then
   !$OMP PARALLEL DO PRIVATE(igpu)
   do igpu=0,gpu_num-1
     call gpu_set_device(igpu)
     call gpu_allocate(cholesky_mo_transp_sp_d(igpu), cholesky_mo_num, mo_num, mo_num)
     call gpu_upload(cholesky_mo_transp_sp,cholesky_mo_transp_sp_d(igpu))
   enddo
   !$OMP END PARALLEL DO
   call gpu_set_device(0)
 endif
END_PROVIDER

