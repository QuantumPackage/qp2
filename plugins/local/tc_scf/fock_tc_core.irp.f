 BEGIN_PROVIDER [ double precision, ao_two_e_tc_core_integral, (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Fock matrices in AO basis set for only the core orbitals 
 ! 
 ! WARNING : use the TCSCF_bi_ort_core_dm_ao matrix  
 END_DOC

 integer :: m,n,l,s,j
 double precision :: integral
 double precision, allocatable :: X(:), X2(:,:,:,:), X3(:,:,:,:)

 allocate (X(cholesky_ao_num))

 ! X(j) = \sum_{mn} TCSCF_bi_ort_core_dm_ao(m,n) * cholesky_ao(m,n,j)
 call dgemm('T','N',cholesky_ao_num,1,ao_num*ao_num,1.d0,            &
     cholesky_ao, ao_num*ao_num,                                     &
     TCSCF_bi_ort_core_dm_ao, ao_num*ao_num,0.d0,                      &
     X, cholesky_ao_num)
!

 ! ao_two_e_tc_core_integral(m,n) = \sum_{j} cholesky_ao(m,n,j) * X(j)
 call dgemm('N','N',ao_num*ao_num,1,cholesky_ao_num, 1.d0,           &
     cholesky_ao, ao_num*ao_num,                                     &
     X, cholesky_ao_num, 0.d0,                                       &
     ao_two_e_tc_core_integral_alpha, ao_num*ao_num)

 deallocate(X)

 double precision :: rss, mem0, mem
 double precision :: memory_of_double

 integer :: iblock
 integer :: block_size

 call resident_memory(mem0)

 block_size = 1024

 rss = memory_of_double(2.d0*ao_num*ao_num)
 do
   mem = mem0 + block_size*rss
   if ( (block_size < 2).or.(mem < qp_max_mem) ) exit
   block_size = block_size/2
 enddo

 call check_mem(block_size*rss, irp_here)

 allocate(X2(ao_num,ao_num,block_size,2))
 allocate(X3(ao_num,block_size,ao_num,2))

! ao_two_e_tc_core_integral (l,s) -= cholesky_ao(l,m,j) * TCSCF_bi_ort_core_dm_ao(m,n) * cholesky_ao(n,s,j)

 do iblock=1,cholesky_ao_num,block_size

   call dgemm('N','N',ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size),ao_num, 1.d0,      &
       TCSCF_bi_ort_core_dm_ao, ao_num,    &
       cholesky_ao(1,1,iblock), ao_num, 0.d0,  &
       X2(1,1,1,1), ao_num)

     do s=1,ao_num
      do j=1,min(cholesky_ao_num-iblock+1,block_size)
       do m=1,ao_num
        X3(m,j,s,1) = X2(m,s,j,1)
       enddo
      enddo
     enddo

   call dgemm('N','N',ao_num,ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size), -1.d0,     &
       cholesky_ao(1,1,iblock), ao_num,       &
       X3(1,1,1,1), ao_num*block_size, 1.d0,  &
       ao_two_e_tc_core_integral, ao_num)

 enddo

 deallocate(X2,X3)
 ao_two_e_tc_core_integral = 2.D0 * ao_two_e_tc_core_integral ! count for alpha + beta

END_PROVIDER

