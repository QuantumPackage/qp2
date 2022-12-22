! ---

 BEGIN_PROVIDER [ double precision, ao_two_e_integral_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_two_e_integral_beta , (ao_num, ao_num) ]

 BEGIN_DOC
 !
 ! 2-e part of alpha and beta Fock matrices (F^{a} & F^{b}) in AO basis set
 ! 
 ! F^{a} = h + G^{a}
 ! F^{b} = h + G^{b}
 ! 
 ! where : 
 !         F^{a} = J^{a} + J^{b} - K^{a} ==> G_{ij}^{a} = \sum_{k,l} P_{kl} (kl|ij) - P_{kl}^{a} (ki|lj)
 !         F^{b} = J^{a} + J^{b} - K^{b} ==> G_{ij}^{b} = \sum_{k,l} P_{kl} (kl|ij) - P_{kl}^{b} (ki|lj)
 ! 
 ! and P_{kl} = P_{kl}^{a} + P_{kl}^{b} 
 !
 END_DOC

 use map_module

 implicit none

 integer                        :: i,j,k,l,k1,r,s
 integer                        :: i0,j0,k0,l0
 integer*8                      :: p,q
 double precision               :: integral, c0, c1, c2
 double precision               :: ao_two_e_integral, local_threshold
 double precision, allocatable  :: ao_two_e_integral_alpha_tmp(:,:)
 double precision, allocatable  :: ao_two_e_integral_beta_tmp(:,:)

 ao_two_e_integral_alpha = 0.d0
 ao_two_e_integral_beta  = 0.d0
 if (do_direct_integrals) then

   !$OMP PARALLEL DEFAULT(NONE)                                      &
       !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,keys,values,p,q,r,s,i0,j0,k0,l0, &
       !$OMP ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp, c0, c1, c2, &
       !$OMP local_threshold)&
       !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
       !$OMP ao_integrals_map,ao_integrals_threshold, ao_two_e_integral_schwartz, &
       !$OMP ao_two_e_integral_alpha, ao_two_e_integral_beta)

   allocate(keys(1), values(1))
   allocate(ao_two_e_integral_alpha_tmp(ao_num,ao_num), &
            ao_two_e_integral_beta_tmp(ao_num,ao_num))
   ao_two_e_integral_alpha_tmp = 0.d0
   ao_two_e_integral_beta_tmp  = 0.d0

   q = ao_num*ao_num*ao_num*ao_num
   !$OMP DO SCHEDULE(static,64)
   do p=1_8,q
           call two_e_integrals_index_reverse(kk,ii,ll,jj,p)
           if ( (kk(1)>ao_num).or. &
                (ii(1)>ao_num).or. &
                (jj(1)>ao_num).or. &
                (ll(1)>ao_num) ) then
                cycle
           endif
           k = kk(1)
           i = ii(1)
           l = ll(1)
           j = jj(1)

           logical, external :: ao_two_e_integral_zero
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

   integer(omp_lock_kind) :: lck(ao_num)
   integer(map_size_kind)     :: i8
   integer                        :: ii(8), jj(8), kk(8), ll(8), k2
   integer(cache_map_size_kind)   :: n_elements_max, n_elements
   integer(key_kind), allocatable :: keys(:)
   double precision, allocatable  :: values(:)

   !$OMP PARALLEL DEFAULT(NONE)                                      &
       !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
       !$OMP  n_elements,ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp)&
       !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
       !$OMP  ao_integrals_map, ao_two_e_integral_alpha, ao_two_e_integral_beta)

   call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
   allocate(keys(n_elements_max), values(n_elements_max))
   allocate(ao_two_e_integral_alpha_tmp(ao_num,ao_num), &
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

END_PROVIDER

! ---

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
