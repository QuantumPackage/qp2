 BEGIN_PROVIDER [ double precision, ao_two_e_integral_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_two_e_integral_beta ,  (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Alpha Fock matrix in ao basis set
 END_DOC

 integer                        :: i,j,k,l,k1,r,s
 integer                        :: i0,j0,k0,l0
 integer*8                      :: p,q
 double precision               :: integral, c0, c1, c2
 double precision               :: ao_two_e_integral, local_threshold
 double precision, allocatable  :: ao_two_e_integral_alpha_tmp(:,:)
 double precision, allocatable  :: ao_two_e_integral_beta_tmp(:,:)
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: ao_two_e_integral_beta_tmp
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: ao_two_e_integral_alpha_tmp

 ao_two_e_integral_alpha = 0.d0
 ao_two_e_integral_beta  = 0.d0
   PROVIDE ao_two_e_integrals_in_map
   PROVIDE ao_two_e_integrals_erf_in_map

   integer(omp_lock_kind) :: lck(ao_num)
   integer*8                      :: i8
   integer                        :: ii(8), jj(8), kk(8), ll(8), k2
   integer(cache_map_size_kind)   :: n_elements_max, n_elements
   integer(key_kind), allocatable :: keys(:)
   double precision, allocatable  :: values(:)
   integer(cache_map_size_kind)   :: n_elements_max_erf, n_elements_erf
   integer(key_kind), allocatable :: keys_erf(:)
   double precision, allocatable  :: values_erf(:)

   !$OMP PARALLEL DEFAULT(NONE) if (ao_num > 100) &
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
   !DIR$ NOVECTOR
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

   !$OMP PARALLEL DEFAULT(NONE) if (ao_num > 100) &
       !$OMP PRIVATE(i,j,l,k1,k,integral_erf,ii,jj,kk,ll,i8,keys_erf,values_erf,n_elements_max_erf, &
       !$OMP  n_elements_erf,ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp)&
       !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
       !$OMP  ao_integrals_erf_map, ao_two_e_integral_alpha, ao_two_e_integral_beta)


   call get_cache_map_n_elements_max(ao_integrals_erf_map,n_elements_max_erf)
   allocate(ao_two_e_integral_alpha_tmp(ao_num,ao_num), &
            ao_two_e_integral_beta_tmp(ao_num,ao_num))
   allocate(keys_Erf(n_elements_max_erf), values_erf(n_elements_max_erf))

   ao_two_e_integral_alpha_tmp = 0.d0
   ao_two_e_integral_beta_tmp  = 0.d0
   !$OMP DO SCHEDULE(static,1)
   !DIR$ NOVECTOR
   do i8=0_8,ao_integrals_erf_map%map_size
     n_elements_erf = n_elements_max_erf
     call get_cache_map(ao_integrals_erf_map,i8,keys_erf,values_erf,n_elements_erf)
     do k1=1,n_elements_erf
       call two_e_integrals_index_reverse(kk,ii,ll,jj,keys_erf(k1))

       do k2=1,8
         if (kk(k2)==0) then
           cycle
         endif
         i = ii(k2)
         j = jj(k2)
         k = kk(k2)
         l = ll(k2)
         double precision :: integral_erf
         integral_erf = values_erf(k1)
         ao_two_e_integral_alpha_tmp(l,j) -= (SCF_density_matrix_ao_alpha(k,i) * integral_erf)
         ao_two_e_integral_beta_tmp (l,j) -= (SCF_density_matrix_ao_beta (k,i) * integral_erf)
       enddo
     enddo
   enddo

   !$OMP END DO NOWAIT
   !$OMP CRITICAL
   ao_two_e_integral_alpha  =  ao_two_e_integral_alpha + ao_two_e_integral_alpha_tmp
   ao_two_e_integral_beta   =  ao_two_e_integral_beta  + ao_two_e_integral_beta_tmp
   !$OMP END CRITICAL
   deallocate(ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp)
   deallocate(keys_erf,values_erf)
   !$OMP END PARALLEL


END_PROVIDER

 BEGIN_PROVIDER [ double precision, Fock_matrix_ao_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_ao_beta,  (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Alpha Fock matrix in ao basis set
 END_DOC

 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     Fock_matrix_ao_alpha(i,j) = Fock_matrix_alpha_no_xc_ao(i,j) + ao_potential_alpha_xc(i,j)
     Fock_matrix_ao_beta(i,j) = Fock_matrix_beta_no_xc_ao(i,j)  + ao_potential_beta_xc(i,j)
   enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [ double precision, Fock_matrix_alpha_no_xc_ao, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_beta_no_xc_ao,  (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Mono electronic an Coulomb matrix in ao basis set
 END_DOC

 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     Fock_matrix_alpha_no_xc_ao(i,j) = ao_one_e_integrals(i,j) + ao_two_e_integral_alpha(i,j)
     Fock_matrix_beta_no_xc_ao(i,j) = ao_one_e_integrals(i,j) + ao_two_e_integral_beta (i,j)
   enddo
 enddo

END_PROVIDER

