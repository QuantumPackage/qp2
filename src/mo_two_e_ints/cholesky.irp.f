BEGIN_PROVIDER [ double precision, cholesky_mo, (mo_num, mo_num, cholesky_ao_num) ]
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in MO basis
 END_DOC

 integer :: k

 !$OMP PARALLEL DO PRIVATE(k)
 do k=1,cholesky_ao_num
  call ao_to_mo(cholesky_ao(1,1,k),ao_num,cholesky_mo(1,1,k),mo_num)
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER

