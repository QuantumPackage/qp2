
 BEGIN_PROVIDER [ double precision, ao_two_e_integral_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_two_e_integral_beta ,  (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Alpha and Beta Fock matrices in AO basis set
 END_DOC

 integer                        :: i,j,k,l,k1,r,s
 integer                        :: i0,j0,k0,l0
 integer*8                      :: p,q
 double precision               :: integral, c0, c1, c2
 double precision               :: ao_two_e_integral, local_threshold
 double precision, allocatable  :: ao_two_e_integral_alpha_tmp(:,:)
 double precision, allocatable  :: ao_two_e_integral_beta_tmp(:,:)

 if (do_ao_cholesky) then   ! Use Cholesky-decomposed integrals
   ao_two_e_integral_alpha(:,:) = ao_two_e_integral_alpha_chol(:,:)
   ao_two_e_integral_beta (:,:) = ao_two_e_integral_beta_chol (:,:)

 else ! Use integrals in AO basis set

   ao_two_e_integral_alpha = 0.d0
   ao_two_e_integral_beta  = 0.d0
   if (do_direct_integrals) then

     !$OMP PARALLEL DEFAULT(NONE)                                    &
         !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,keys,values,p,q,r,s,i0,j0,k0,l0,&
         !$OMP ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp, c0, c1, c2,&
         !$OMP local_threshold)                                      &
         !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
         !$OMP ao_integrals_map,ao_integrals_threshold, ao_two_e_integral_schwartz,&
         !$OMP ao_two_e_integral_alpha, ao_two_e_integral_beta)

     allocate(keys(1), values(1))
     allocate(ao_two_e_integral_alpha_tmp(ao_num,ao_num),            &
         ao_two_e_integral_beta_tmp(ao_num,ao_num))
     ao_two_e_integral_alpha_tmp = 0.d0
     ao_two_e_integral_beta_tmp  = 0.d0

     q = ao_num*ao_num*ao_num*ao_num
     !$OMP DO SCHEDULE(static,64)
     do p=1_8,q
       call two_e_integrals_index_reverse(kk,ii,ll,jj,p)
       if ( (kk(1)>ao_num).or.                                       &
             (ii(1)>ao_num).or.                                      &
             (jj(1)>ao_num).or.                                      &
             (ll(1)>ao_num) ) then
         cycle
       endif
       k = kk(1)
       i = ii(1)
       l = ll(1)
       j = jj(1)

       logical, external              :: ao_two_e_integral_zero
       if (ao_two_e_integral_zero(i,k,j,l)) then
         cycle
       endif
       local_threshold = ao_two_e_integral_schwartz(k,l)*ao_two_e_integral_schwartz(i,j)
       if (local_threshold < ao_integrals_threshold) then
         cycle
       endif
       i0 = i
       j0 = j
       k0 = k
       l0 = l
       values(1) = 0.d0
       local_threshold = ao_integrals_threshold/local_threshold
       do k2=1,8
         if (kk(k2)==0) then
           cycle
         endif
         i = ii(k2)
         j = jj(k2)
         k = kk(k2)
         l = ll(k2)
         c0 = SCF_density_matrix_ao_alpha(k,l)+SCF_density_matrix_ao_beta(k,l)
         c1 = SCF_density_matrix_ao_alpha(k,i)
         c2 = SCF_density_matrix_ao_beta(k,i)
         if ( dabs(c0)+dabs(c1)+dabs(c2) < local_threshold) then
           cycle
         endif
         if (values(1) == 0.d0) then
           values(1) = ao_two_e_integral(k0,l0,i0,j0)
         endif
         integral = c0 * values(1)
         ao_two_e_integral_alpha_tmp(i,j) += integral
         ao_two_e_integral_beta_tmp (i,j) += integral
         integral = values(1)
         ao_two_e_integral_alpha_tmp(l,j) -= c1 * integral
         ao_two_e_integral_beta_tmp (l,j) -= c2 * integral
       enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP CRITICAL
     ao_two_e_integral_alpha += ao_two_e_integral_alpha_tmp
     ao_two_e_integral_beta  += ao_two_e_integral_beta_tmp
     !$OMP END CRITICAL
     deallocate(keys,values,ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp)
     !$OMP END PARALLEL
   else
     PROVIDE ao_two_e_integrals_in_map

     integer(omp_lock_kind)         :: lck(ao_num)
     integer(map_size_kind)         :: i8
     integer                        :: ii(8), jj(8), kk(8), ll(8), k2
     integer(cache_map_size_kind)   :: n_elements_max, n_elements
     integer(key_kind), allocatable :: keys(:)
     double precision, allocatable  :: values(:)

     !$OMP PARALLEL DEFAULT(NONE)                                    &
         !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max,&
         !$OMP  n_elements,ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp)&
         !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
         !$OMP  ao_integrals_map, ao_two_e_integral_alpha, ao_two_e_integral_beta)

     call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
     allocate(keys(n_elements_max), values(n_elements_max))
     allocate(ao_two_e_integral_alpha_tmp(ao_num,ao_num),            &
         ao_two_e_integral_beta_tmp(ao_num,ao_num))
     ao_two_e_integral_alpha_tmp = 0.d0
     ao_two_e_integral_beta_tmp  = 0.d0

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
           integral = (SCF_density_matrix_ao_alpha(k,l)+SCF_density_matrix_ao_beta(k,l)) * values(k1)
           ao_two_e_integral_alpha_tmp(i,j) += integral
           ao_two_e_integral_beta_tmp (i,j) += integral
           integral = values(k1)
           ao_two_e_integral_alpha_tmp(l,j) -= SCF_density_matrix_ao_alpha(k,i) * integral
           ao_two_e_integral_beta_tmp (l,j) -= SCF_density_matrix_ao_beta (k,i) * integral
         enddo
       enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP CRITICAL
     ao_two_e_integral_alpha += ao_two_e_integral_alpha_tmp
     ao_two_e_integral_beta  += ao_two_e_integral_beta_tmp
     !$OMP END CRITICAL
     deallocate(keys,values,ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp)
     !$OMP END PARALLEL

   endif
 endif

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_two_e_integral_alpha_chol, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_two_e_integral_beta_chol ,  (ao_num, ao_num) ]
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
     SCF_density_matrix_ao, ao_num*ao_num,0.d0,                      &
     X, cholesky_ao_num)
!

 ! ao_two_e_integral_alpha(m,n) = \sum_{j} cholesky_ao(m,n,j) * X(j)
 call dgemm('N','N',ao_num*ao_num,1,cholesky_ao_num, 1.d0,           &
     cholesky_ao, ao_num*ao_num,                                     &
     X, cholesky_ao_num, 0.d0,                                       &
     ao_two_e_integral_alpha_chol, ao_num*ao_num)

 deallocate(X)

 if (elec_alpha_num > elec_beta_num) then
   ao_two_e_integral_beta_chol = ao_two_e_integral_alpha_chol
 endif


 double precision :: rss
 double precision :: memory_of_double

 integer :: iblock
 integer, parameter :: block_size = 32

 rss = memory_of_double(ao_num*ao_num)
 call check_mem(2.d0*block_size*rss, irp_here)
 allocate(X2(ao_num,ao_num,block_size,2))
 allocate(X3(ao_num,block_size,ao_num,2))
    
! ao_two_e_integral_alpha_chol (l,s) -= cholesky_ao(l,m,j) * SCF_density_matrix_ao_beta (m,n) * cholesky_ao(n,s,j)

 do iblock=1,cholesky_ao_num,block_size

   call dgemm('N','N',ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size),ao_num, 1.d0,      &
       SCF_density_matrix_ao_alpha, ao_num,    &
       cholesky_ao(1,1,iblock), ao_num, 0.d0,  &
       X2(1,1,1,1), ao_num)

   if (elec_alpha_num > elec_beta_num) then
     call dgemm('N','N',ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size),ao_num, 1.d0,      &
       SCF_density_matrix_ao_beta, ao_num,     &
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
       ao_two_e_integral_alpha_chol, ao_num)

   if (elec_alpha_num > elec_beta_num) then
     call dgemm('N','N',ao_num,ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size), -1.d0,     &
       cholesky_ao(1,1,iblock), ao_num,        &
       X3(1,1,1,2), ao_num*block_size, 1.d0,   &
       ao_two_e_integral_beta_chol, ao_num)
   endif

 enddo

 if (elec_alpha_num == elec_beta_num) then
   ao_two_e_integral_beta_chol = ao_two_e_integral_alpha_chol
 endif
 deallocate(X2,X3)

END_PROVIDER

 BEGIN_PROVIDER [ double precision, Fock_matrix_ao_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_ao_beta,  (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Alpha Fock matrix in AO basis set
 END_DOC

 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     Fock_matrix_ao_alpha(i,j) = ao_one_e_integrals(i,j) + ao_two_e_integral_alpha(i,j)
     Fock_matrix_ao_beta (i,j) = ao_one_e_integrals(i,j) + ao_two_e_integral_beta (i,j)
   enddo
 enddo

END_PROVIDER
