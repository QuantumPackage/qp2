BEGIN_TEMPLATE 

double precision function get_ao$_erf_integ_chol(i,j,k,l)
 implicit none
  BEGIN_DOC
  !  CHOLESKY representation of the integral of the AO$_erf basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC
 integer, intent(in) :: i,j,k,l
 double precision, external :: ddot
 get_ao$_erf_integ_chol = ddot(cholesky_ao$_erf_num, cholesky_ao$_erf_transp(1,i,j), 1, cholesky_ao$_erf_transp(1,k,l), 1)

end

 BEGIN_PROVIDER [ integer, cholesky_ao$_erf_num ]
 implicit none
 BEGIN_DOC
 ! Number of Cholesky vectors in MO basis
 END_DOC
 integer, external              :: getUnitAndOpen
 integer                        :: iunit
 if (read_ao$_erf_cholesky) then
   iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_ao$_erf', 'R')
   read(iunit) cholesky_ao$_erf_num
   close(iunit)
 else
   cholesky_ao$_erf_num = cholesky_ao$_erf_cart_num
 endif
END_PROVIDER

BEGIN_PROVIDER [ double precision, cholesky_ao$_erf, (ao_num, ao_num, cholesky_ao$_erf_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in ao basis
 END_DOC

 integer :: k, i, j

 call set_multiple_levels_omp(.False.)
 !$OMP PARALLEL DO PRIVATE(k)
 do k=1,cholesky_ao$_erf_num
  do j=1,ao_num
    do i=1,ao_num
      cholesky_ao$_erf(i,j,k) = cholesky_ao$_erf_transp(k,i,j)
    enddo
  enddo
 enddo
 !$OMP END PARALLEL DO
 free cholesky_ao$_erf_transp

END_PROVIDER

BEGIN_PROVIDER [ double precision, cholesky_ao$_erf_transp, (cholesky_ao$_erf_num, ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in ao basis. Warning: it is transposed wrt cholesky_ao$_erf_cart:
 !
 ! - cholesky_ao$_erf_cart        is (ao_cart_num^2 x cholesky_ao$_erf_cart_num)
 !
 ! - cholesky_ao$_erf_transp is (cholesky_ao$_erf_num x ao_num^2)
 END_DOC

 double precision, allocatable :: X(:,:,:)
 double precision :: wall0, wall1
 integer, external              :: getUnitAndOpen
 integer                        :: iunit, ierr, rank

  if (read_ao$_erf_cholesky) then
      print *,  'Reading Cholesky ao$_erf vectors from disk...'
     iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_ao$_erf_transp', 'R')
     read(iunit) rank
     if (cholesky_ao$_erf_num /= rank) then
        stop 'inconsistent rank'
     endif
     read(iunit) cholesky_ao$_erf_transp
     close(iunit)
  else
     print *, ''
     print *, 'ao$_erf_cart->ao$_erf Transformation of Cholesky vectors'
     print *, '-----------------------------------------'
     print *, ''

     call wall_time(wall0)

     allocate(X(ao_num,cholesky_ao$_erf_num,ao_cart_num), stat=ierr)
     if (ierr /= 0) then
       print *, irp_here, ': Allocation failed'
     endif
     call dgemm('T','N', ao_cart_num*cholesky_ao$_erf_num, ao_num, ao_cart_num, 1.d0, &
         cholesky_ao$_erf_cart, ao_cart_num, ao_cart_to_ao_basis_mat_transp, ao_cart_num, & 
         0.d0, X, ao_cart_num*cholesky_ao$_erf_num)
     call dgemm('T','N', cholesky_ao$_erf_num*ao_num, ao_num, ao_cart_num, 1.d0, &
         X, ao_cart_num, ao_cart_to_ao_basis_mat_transp, ao_cart_num, 0.d0, & 
         cholesky_ao$_erf_transp, cholesky_ao$_erf_num*ao_num)
     deallocate(X)
     call wall_time(wall1)
     print*,'Time to provide ao$_erf cholesky vectors = ',(wall1-wall0)/60.d0, ' min'


     if (write_ao$_erf_cholesky) then
       print *,  'Writing Cholesky ao$_erf vectors to disk...'
       iunit = getUnitAndOpen(trim(ezfio_work_dir)//'cholesky_ao$_erf_transp', 'W')
       write(iunit) cholesky_ao$_erf_num
       write(iunit) cholesky_ao$_erf_transp
       close(iunit)
       call ezfio_set_ao_two_e_ints_io_ao$_erf_cholesky('Read')
     endif
  endif

END_PROVIDER


BEGIN_PROVIDER [ real, cholesky_ao$_erf_sp, (ao_num, ao_num, cholesky_ao$_erf_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in ao$_erf basis in stored in single precision
 END_DOC

 integer :: k, i, j

 call set_multiple_levels_omp(.False.)
 !$OMP PARALLEL DO PRIVATE(k)
 do k=1,cholesky_ao$_erf_num
  do j=1,ao_num
    do i=1,ao_num
      cholesky_ao$_erf_sp(i,j,k) = cholesky_ao$_erf_transp_sp(k,i,j)
    enddo
  enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER

BEGIN_PROVIDER [ real, cholesky_ao$_erf_transp_sp, (cholesky_ao$_erf_num, ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in ao$_erf basis in s. Warning: it is transposed wrt cholesky_ao$_erf_cart:
 !
 ! - cholesky_ao$_erf_cart        is (ao_cart_num^2 x cholesky_ao$_erf_cart_num)
 !
 ! - cholesky_ao$_erf_transp is (cholesky_ao$_erf_num x ao_num^2)
 END_DOC

 integer :: i,j,k
 !$OMP PARALLEL DO PRIVATE(k)
 do j=1,ao_num
  do i=1,ao_num
   do k=1,cholesky_ao$_erf_num
      cholesky_ao$_erf_transp_sp(k,i,j) = cholesky_ao$_erf_transp(k,i,j)
    enddo
  enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER

SUBST [ _erf ]
  
;;
_erf;;
_cgtos;;

END_TEMPLATE 
