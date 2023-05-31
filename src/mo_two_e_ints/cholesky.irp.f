BEGIN_PROVIDER [ double precision, cholesky_mo, (mo_num, mo_num, cholesky_ao_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in MO basis
 END_DOC

 integer :: k

 print *, 'AO->MO Transformation of Cholesky vectors'
 !$OMP PARALLEL DO PRIVATE(k)
 do k=1,cholesky_ao_num
  call ao_to_mo(cholesky_ao(1,1,k),ao_num,cholesky_mo(1,1,k),mo_num)
 enddo
 !$OMP END PARALLEL DO
 print *, ''

END_PROVIDER

BEGIN_PROVIDER [ double precision, cholesky_mo_transp, (cholesky_ao_num, mo_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in MO basis
 END_DOC

 integer :: i,j,k
 double precision, allocatable :: buffer(:,:)

 print *, 'AO->MO Transformation of Cholesky vectors  .'
 !$OMP PARALLEL PRIVATE(i,j,k,buffer)
 allocate(buffer(mo_num,mo_num))
 !$OMP DO SCHEDULE(static)
 do k=1,cholesky_ao_num
  call ao_to_mo(cholesky_ao(1,1,k),ao_num,buffer,mo_num)
  do j=1,mo_num
    do i=1,mo_num
      cholesky_mo_transp(k,i,j) = buffer(i,j)
    enddo
  enddo
 enddo
 !$OMP END DO
 deallocate(buffer)
 !$OMP END PARALLEL
 print *, ''

END_PROVIDER

