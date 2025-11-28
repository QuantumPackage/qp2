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

 call set_multiple_levels_omp(.False.)
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
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in MO basis. Warning: it is transposed wrt cholesky_ao:
 !
 ! - cholesky_ao        is (ao_num^2 x cholesky_ao_num)
 !
 ! - cholesky_mo_transp is (cholesky_mo_num x mo_num^2)
 END_DOC

 double precision, allocatable :: X(:,:,:)
 double precision :: wall0, wall1
 integer, external              :: getUnitAndOpen
 integer                        :: iunit, ierr, rank

  if (read_mo_cholesky) then
      print *,  'Reading Cholesky MO vectors from disk...'
     iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_mo_transp', 'R')
     read(iunit) rank
     if (cholesky_mo_num /= rank) then
        stop 'inconsistent rank'
     endif
     read(iunit) cholesky_mo_transp
     close(iunit)
  else
     print *, ''
     print *, 'AO->MO Transformation of Cholesky vectors'
     print *, '-----------------------------------------'
     print *, ''

     call wall_time(wall0)

     allocate(X(mo_num,cholesky_mo_num,ao_num), stat=ierr)
     if (ierr /= 0) then
       print *, irp_here, ': Allocation failed'
     endif
     call dgemm('T','N', ao_num*cholesky_mo_num, mo_num, ao_num, 1.d0, &
         cholesky_ao, ao_num, mo_coef, ao_num, 0.d0, X, ao_num*cholesky_mo_num)
     call dgemm('T','N', cholesky_mo_num*mo_num, mo_num, ao_num, 1.d0, &
         X, ao_num, mo_coef, ao_num, 0.d0, cholesky_mo_transp, cholesky_mo_num*mo_num)
     deallocate(X)
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

 call set_multiple_levels_omp(.False.)
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


