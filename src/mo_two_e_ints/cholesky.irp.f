BEGIN_PROVIDER [ integer, cholesky_mo_num ]
 implicit none
 BEGIN_DOC
 ! Number of Cholesky vectors in MO basis
 END_DOC
 cholesky_mo_num = cholesky_ao_num
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
 ! Cholesky vectors in MO basis
 END_DOC

 double precision, allocatable :: X(:,:,:)
 integer :: ierr
 print *, 'AO->MO Transformation of Cholesky vectors'

 allocate(X(mo_num,cholesky_mo_num,ao_num), stat=ierr)
 if (ierr /= 0) then
   print *, irp_here, ': Allocation failed'
 endif
 call dgemm('T','N', ao_num*cholesky_mo_num, mo_num, ao_num, 1.d0, &
     cholesky_ao, ao_num, mo_coef, ao_num, 0.d0, X, ao_num*cholesky_mo_num)
 call dgemm('T','N', cholesky_mo_num*mo_num, mo_num, ao_num, 1.d0, &
     X, ao_num, mo_coef, ao_num, 0.d0, cholesky_mo_transp, cholesky_mo_num*mo_num)
 deallocate(X)

END_PROVIDER

