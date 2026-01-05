BEGIN_PROVIDER [double precision, Fock_matrix_tc_mo_core_eri, (mo_num, mo_num)]

  BEGIN_DOC
  ! Two-electron part with BARE COULOMB interaction of the Fock matrix 
  ! 
  ! with contributions coming only from the CORE ELECTRONS
  END_DOC

  implicit none
  double precision              :: t0, t1, tt0, tt1
  double precision, allocatable :: tmp(:,:)

  PROVIDE mo_l_coef mo_r_coef

  call ao_to_mo_bi_ortho( ao_two_e_core_integral, size(ao_two_e_core_integral, 1) &
                        ,  Fock_matrix_tc_mo_core_eri, size( Fock_matrix_tc_mo_core_eri, 1) )


END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, Fock_matrix_alpha_tc_mo_valence_eri, (mo_num, mo_num)]

  BEGIN_DOC
  ! Two-electron part with BARE COULOMB interaction of the Fock matrix 
  ! 
  ! with contributions coming only from the VALENCE ELECTRONS
  END_DOC

  implicit none
  double precision              :: t0, t1, tt0, tt1
  double precision, allocatable :: tmp(:,:)

  PROVIDE mo_l_coef mo_r_coef

  call ao_to_mo_bi_ortho( ao_two_e_valence_integral_alpha, size(ao_two_e_valence_integral_alpha, 1) &
                        ,  Fock_matrix_alpha_tc_mo_valence_eri, size(Fock_matrix_alpha_tc_mo_valence_eri, 1) )


END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, Fock_matrix_beta_tc_mo_valence_eri, (mo_num, mo_num)]

  BEGIN_DOC
  ! Two-electron part with BARE COULOMB interaction of the Fock matrix 
  ! 
  ! with contributions coming only from the VALENCE ELECTRONS
  END_DOC

  implicit none
  double precision              :: t0, t1, tt0, tt1
  double precision, allocatable :: tmp(:,:)

  PROVIDE mo_l_coef mo_r_coef

  call ao_to_mo_bi_ortho( ao_two_e_valence_integral_beta, size(ao_two_e_valence_integral_beta, 1) &
                        ,  Fock_matrix_beta_tc_mo_valence_eri, size(Fock_matrix_beta_tc_mo_valence_eri, 1) )


END_PROVIDER

! ---


BEGIN_PROVIDER [ double precision, Fock_matrix_tc_eri_mo_valence, (mo_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Fock matrix containing only the 1/r12 interaction and only the VALENCE density matrix
 END_DOC
 integer :: i,j
 if (all_shells_closed) then
   Fock_matrix_tc_eri_mo_valence = Fock_matrix_alpha_tc_mo_valence_eri
 else
   ! Core
   do j = 1, elec_beta_num
     ! Core
     do i = 1, elec_beta_num
       Fock_matrix_tc_eri_mo_valence(i,j) = Fock_matrix_param(1,1) * Fock_matrix_alpha_tc_mo_valence_eri(i,j) + &
                             Fock_matrix_param(2,1) * Fock_matrix_beta_tc_mo_valence_eri(i,j)
     enddo
     ! Open
     do i = elec_beta_num+1, elec_alpha_num
       Fock_matrix_tc_eri_mo_valence(i,j) = Fock_matrix_beta_tc_mo_valence_eri(i,j)
     enddo
     ! Virtual
     do i = elec_alpha_num+1, mo_num
       Fock_matrix_tc_eri_mo_valence(i,j) = 0.5d0 * Fock_matrix_alpha_tc_mo_valence_eri(i,j) + &
                             0.5d0 * Fock_matrix_beta_tc_mo_valence_eri(i,j)
     enddo
   enddo
   ! Open
   do j = elec_beta_num+1, elec_alpha_num
     ! Core
     do i = 1, elec_beta_num
       Fock_matrix_tc_eri_mo_valence(i,j) = Fock_matrix_beta_tc_mo_valence_eri(i,j)
     enddo
     ! Open
     do i = elec_beta_num+1, elec_alpha_num
       Fock_matrix_tc_eri_mo_valence(i,j) = Fock_matrix_param(1,2) * Fock_matrix_alpha_tc_mo_valence_eri(i,j) + &
                             Fock_matrix_param(2,2) * Fock_matrix_beta_tc_mo_valence_eri(i,j)
     enddo
     ! Virtual
     do i = elec_alpha_num+1, mo_num
       Fock_matrix_tc_eri_mo_valence(i,j) =  Fock_matrix_alpha_tc_mo_valence_eri(i,j)
     enddo
   enddo
   ! Virtual
   do j = elec_alpha_num+1, mo_num
     ! Core
     do i = 1, elec_beta_num
       Fock_matrix_tc_eri_mo_valence(i,j) =   0.5d0 * Fock_matrix_alpha_tc_mo_valence_eri(i,j) &
                             + 0.5d0 * Fock_matrix_beta_tc_mo_valence_eri(i,j)
     enddo
     ! Open
     do i = elec_beta_num+1, elec_alpha_num
       Fock_matrix_tc_eri_mo_valence(i,j) = Fock_matrix_alpha_tc_mo_valence_eri(i,j)
     enddo
     ! Virtual
     do i = elec_alpha_num+1, mo_num
       Fock_matrix_tc_eri_mo_valence(i,j) = Fock_matrix_param(1,3) * Fock_matrix_alpha_tc_mo_valence_eri(i,j) + &
                             Fock_matrix_param(2,3) * Fock_matrix_beta_tc_mo_valence_eri(i,j)
     enddo
   enddo
 endif
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, ao_two_e_core_integral_chol, (ao_num, ao_num) ]
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

 ! ao_two_e_core_integral_chol(m,n) = \sum_{j} cholesky_ao(m,n,j) * X(j)
 call dgemm('N','N',ao_num*ao_num,1,cholesky_ao_num, 1.d0,           &
     cholesky_ao, ao_num*ao_num,                                     &
     X, cholesky_ao_num, 0.d0,                                       &
     ao_two_e_core_integral_chol, ao_num*ao_num)

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

! ao_two_e_core_integral_chol (l,s) -= cholesky_ao(l,m,j) * TCSCF_bi_ort_core_dm_ao(m,n) * cholesky_ao(n,s,j)

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
       ao_two_e_core_integral_chol, ao_num)

 enddo

 deallocate(X2,X3)
 ao_two_e_core_integral_chol = 2.D0 * ao_two_e_core_integral_chol ! count for alpha + beta

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_two_e_valence_integral_alpha_chol, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_two_e_valence_integral_beta_chol ,  (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Alpha and Beta Fock matrices in AO basis set
 END_DOC

 integer :: m,n,l,s,j
 double precision :: integral
 double precision, allocatable :: X(:), X2(:,:,:,:), X3(:,:,:,:)

 allocate (X(cholesky_ao_num))

 ! X(j) = \sum_{mn} SCF_density_matrix_ao(m,n) * cholesky_ao(m,n,j)
 call dgemm('T','N',cholesky_ao_num,1,ao_num*ao_num,1.d0,            &
     cholesky_ao, ao_num*ao_num,                                     &
     TCSCF_bi_ort_valence_dm_ao_full, ao_num*ao_num,0.d0,                      &
     X, cholesky_ao_num)
!

 ! ao_two_e_integral_alpha(m,n) = \sum_{j} cholesky_ao(m,n,j) * X(j)
 call dgemm('N','N',ao_num*ao_num,1,cholesky_ao_num, 1.d0,           &
     cholesky_ao, ao_num*ao_num,                                     &
     X, cholesky_ao_num, 0.d0,                                       &
     ao_two_e_valence_integral_alpha_chol, ao_num*ao_num)

 deallocate(X)

 if (elec_alpha_num > elec_beta_num) then
   ao_two_e_valence_integral_beta_chol = ao_two_e_valence_integral_alpha_chol
 endif


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

! ao_two_e_valence_integral_alpha_chol (l,s) -= cholesky_ao(l,m,j) * SCF_density_matrix_ao_beta (m,n) * cholesky_ao(n,s,j)

 do iblock=1,cholesky_ao_num,block_size

   call dgemm('N','N',ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size),ao_num, 1.d0,      &
       TCSCF_bi_ort_valence_dm_ao_alpha, ao_num,    &
       cholesky_ao(1,1,iblock), ao_num, 0.d0,  &
       X2(1,1,1,1), ao_num)

   if (elec_alpha_num > elec_beta_num) then
     call dgemm('N','N',ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size),ao_num, 1.d0,      &
       TCSCF_bi_ort_valence_dm_ao_beta, ao_num,     &
       cholesky_ao(1,1,iblock), ao_num, 0.d0,  &
       X2(1,1,1,2), ao_num)

     do s=1,ao_num
      do j=1,min(cholesky_ao_num-iblock+1,block_size)
       do m=1,ao_num
        X3(m,j,s,1) = X2(m,s,j,1)
        X3(m,j,s,2) = X2(m,s,j,2)
       enddo
      enddo
     enddo

   else

     do s=1,ao_num
      do j=1,min(cholesky_ao_num-iblock+1,block_size)
       do m=1,ao_num
        X3(m,j,s,1) = X2(m,s,j,1)
       enddo
      enddo
     enddo
   endif

   call dgemm('N','N',ao_num,ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size), -1.d0,     &
       cholesky_ao(1,1,iblock), ao_num,       &
       X3(1,1,1,1), ao_num*block_size, 1.d0,  &
       ao_two_e_valence_integral_alpha_chol, ao_num)

   if (elec_alpha_num > elec_beta_num) then
     call dgemm('N','N',ao_num,ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size), -1.d0,     &
       cholesky_ao(1,1,iblock), ao_num,        &
       X3(1,1,1,2), ao_num*block_size, 1.d0,   &
       ao_two_e_valence_integral_beta_chol, ao_num)
   endif

 enddo

 if (elec_alpha_num == elec_beta_num) then
   ao_two_e_valence_integral_beta_chol = ao_two_e_valence_integral_alpha_chol
 endif
 deallocate(X2,X3)

END_PROVIDER

!!!!!!!!!!!!!!!!
 BEGIN_PROVIDER [ double precision, ao_two_e_valence_integral_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_two_e_valence_integral_beta ,  (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Fock matrices in AO basis set for only the core orbitals 
 ! 
 ! WARNING : use the TCSCF_bi_ort_core_dm_ao matrix  
 END_DOC
 integer(omp_lock_kind)         :: lck(ao_num)
 integer(map_size_kind)         :: i8
 integer                        :: ii(8), jj(8), kk(8), ll(8), k2,i,j,k,l,k1
 integer(cache_map_size_kind)   :: n_elements_max, n_elements
 integer(key_kind), allocatable :: keys(:)
 double precision, allocatable  :: values(:)
 double precision, allocatable  :: ao_two_e_valence_integral_alpha_tmp(:,:)
 double precision, allocatable  :: ao_two_e_valence_integral_beta_tmp(:,:)
 double precision :: integral

 ao_two_e_valence_integral_alpha = 0.d0
 ao_two_e_valence_integral_beta = 0.d0
 PROVIDE ao_two_e_integrals_in_map

 !$OMP PARALLEL DEFAULT(NONE)                                    &
     !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max,&
     !$OMP  n_elements,ao_two_e_valence_integral_alpha_tmp,ao_two_e_valence_integral_beta_tmp)&
     !$OMP SHARED(ao_num,TCSCF_bi_ort_valence_dm_ao_beta,TCSCF_bi_ort_valence_dm_ao_alpha,&
     !$OMP  ao_integrals_map, ao_two_e_valence_integral_alpha, ao_two_e_valence_integral_beta)

 call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
 allocate(keys(n_elements_max), values(n_elements_max))
 allocate(ao_two_e_valence_integral_alpha_tmp(ao_num,ao_num),            &
     ao_two_e_valence_integral_beta_tmp(ao_num,ao_num))
 ao_two_e_valence_integral_alpha_tmp = 0.d0
 ao_two_e_valence_integral_beta_tmp  = 0.d0

 !$OMP DO SCHEDULE(static,1)
 do i8=0_8,ao_integrals_map%map_size
   n_elements = n_elements_max
   call get_cache_map(ao_integrals_map,i8,keys,values,n_elements)
   do k1=1,n_elements
     call two_e_integrals_index_reverse(kk,ii,ll,jj,keys(k1))

     do k2=1,8
       if (kk(k2)==0) then
         cycle
       endif
       i = ii(k2)
       j = jj(k2)
       k = kk(k2)
       l = ll(k2)
       integral = (TCSCF_bi_ort_valence_dm_ao_alpha(k,l)+TCSCF_bi_ort_valence_dm_ao_beta(k,l)) * values(k1)
       ao_two_e_valence_integral_alpha_tmp(i,j) += integral
       ao_two_e_valence_integral_beta_tmp (i,j) += integral
       integral = values(k1)
       ao_two_e_valence_integral_alpha_tmp(l,j) -= TCSCF_bi_ort_valence_dm_ao_alpha(k,i) * integral
       ao_two_e_valence_integral_beta_tmp (l,j) -= TCSCF_bi_ort_valence_dm_ao_beta(k,i)  * integral
     enddo
   enddo
 enddo
 !$OMP END DO NOWAIT
 !$OMP CRITICAL
 ao_two_e_valence_integral_alpha += ao_two_e_valence_integral_alpha_tmp
 ao_two_e_valence_integral_beta  += ao_two_e_valence_integral_beta_tmp
 !$OMP END CRITICAL
 deallocate(keys,values,ao_two_e_valence_integral_alpha_tmp,ao_two_e_valence_integral_beta_tmp)
 !$OMP END PARALLEL

END_PROVIDER
!!!!!!!!!!!!!!!! 
 BEGIN_PROVIDER [ double precision, ao_two_e_core_integral, (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Fock matrices in AO basis set for only the core orbitals 
 ! 
 ! WARNING : use the TCSCF_bi_ort_core_dm_ao matrix  
 END_DOC
 ao_two_e_core_integral = 0.d0
     PROVIDE ao_two_e_integrals_in_map

 integer                        :: i,j,k,l,k1,r,s
 integer                        :: i0,j0,k0,l0
 integer*8                      :: p,q
 double precision               :: integral, c0, c1, c2
 double precision               :: ao_two_e_integral, local_threshold
 double precision, allocatable  :: ao_two_e_integral_alpha_tmp(:,:)

 double precision               :: dm_core,exch_dm
 integer(omp_lock_kind)         :: lck(ao_num)
 integer(map_size_kind)         :: i8
 integer                        :: ii(8), jj(8), kk(8), ll(8), k2
 integer(cache_map_size_kind)   :: n_elements_max, n_elements
 integer(key_kind), allocatable :: keys(:)
 double precision, allocatable  :: values(:)

 !$OMP PARALLEL DEFAULT(NONE)                                    &
     !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max,&
     !$OMP  n_elements,ao_two_e_integral_alpha_tmp,dm_core,exch_dm)&
     !$OMP SHARED(ao_num,TCSCF_bi_ort_core_dm_ao,&
     !$OMP  ao_integrals_map, ao_two_e_core_integral)

 call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
 allocate(keys(n_elements_max), values(n_elements_max))
 allocate(ao_two_e_integral_alpha_tmp(ao_num,ao_num))
 ao_two_e_integral_alpha_tmp = 0.d0

 !$OMP DO SCHEDULE(static,1)
 do i8=0_8,ao_integrals_map%map_size
   n_elements = n_elements_max
   call get_cache_map(ao_integrals_map,i8,keys,values,n_elements)
   do k1=1,n_elements
     call two_e_integrals_index_reverse(kk,ii,ll,jj,keys(k1))

     do k2=1,8
       if (kk(k2)==0) then
         cycle
       endif
       i = ii(k2)
       j = jj(k2)
       k = kk(k2)
       l = ll(k2)
       dm_core = TCSCF_bi_ort_core_dm_ao(k,l)
       exch_dm = 0.5d0 * TCSCF_bi_ort_core_dm_ao(k,i)
       integral = dm_core * values(k1)
       ao_two_e_integral_alpha_tmp(i,j) += integral
       integral = values(k1)
       ao_two_e_integral_alpha_tmp(l,j) -= exch_dm * integral
     enddo
   enddo
 enddo
 !$OMP END DO NOWAIT
 !$OMP CRITICAL
 ao_two_e_core_integral += ao_two_e_integral_alpha_tmp
 !$OMP END CRITICAL
 deallocate(keys,values,ao_two_e_integral_alpha_tmp)
 !$OMP END PARALLEL
 END_PROVIDER 
