
 BEGIN_PROVIDER [ complex*16, ao_two_e_integral_alpha_complex, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ complex*16, ao_two_e_integral_beta_complex ,  (ao_num, ao_num) ]
  use map_module
  implicit none
  BEGIN_DOC
  ! Alpha and Beta Fock matrices in AO basis set
  END_DOC
  !TODO: finish implementing this: see complex qp1 (different mapping)
 
  integer                        :: i,j,k,l,k1,r,s
  integer                        :: i0,j0,k0,l0
  integer*8                      :: p,q
  complex*16               :: integral, c0
  double precision               :: ao_two_e_integral, local_threshold
  double precision, allocatable  :: ao_two_e_integral_alpha_tmp(:,:)
  double precision, allocatable  :: ao_two_e_integral_beta_tmp(:,:)
 
  ao_two_e_integral_alpha_complex = (0.d0,0.d0)
  ao_two_e_integral_beta_complex  = (0.d0,0.d0)
  PROVIDE ao_two_e_integrals_in_map
 
  integer(omp_lock_kind) :: lck(ao_num)
  integer(map_size_kind)     :: i8
  integer                        :: ii(4), jj(4), kk(4), ll(4), k2
  integer(cache_map_size_kind)   :: n_elements_max, n_elements
  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  complex*16, parameter    :: i_sign(4) = (/(0.d0,1.d0),(0.d0,1.d0),(0.d0,-1.d0),(0.d0,-1.d0)/)
  integer(key_kind) :: key1
 
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
      !$OMP  n_elements,ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp, &
      !$OMP  c0,key1)&
      !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha_complex, &
      !$OMP  SCF_density_matrix_ao_beta_complex,i_sign, &
      !$OMP  ao_integrals_map, ao_two_e_integral_alpha_complex, ao_two_e_integral_beta_complex)
 
  call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
  allocate(ao_two_e_integral_alpha_tmp(ao_num,ao_num), &
           ao_two_e_integral_beta_tmp(ao_num,ao_num))
  ao_two_e_integral_alpha_tmp = (0.d0,0.d0)
  ao_two_e_integral_beta_tmp  = (0.d0,0.d0)
 
  !$OMP DO SCHEDULE(static,1)
  do i8=0_8,ao_integrals_map%map_size
    n_elements = n_elements_max
    call get_cache_map(ao_integrals_map,i8,keys,values,n_elements)
    do k1=1,n_elements
      ! get original key
      ! reverse of 2*key (imag part) and 2*key-1 (real part)
      key1 = shiftr(keys(k1)+1,1)

      call two_e_integrals_index_reverse_complex_1(ii,jj,kk,ll,key1)
      ! i<=k, j<=l, ik<=jl
      ! ijkl, jilk, klij*, lkji*

      if (shiftl(key1,1)==keys(k1)) then !imaginary part (even)
        do k2=1,4
          if (ii(k2)==0) then
            cycle
          endif
          i = ii(k2)
          j = jj(k2)
          k = kk(k2)
          l = ll(k2)
          integral = i_sign(k2)*values(k1) !for klij and lkji, take complex conjugate

          !G_a(i,k) += D_{ab}(l,j)*(<ij|kl>)
          !G_b(i,k) += D_{ab}(l,j)*(<ij|kl>)
          !G_a(i,l) -= D_a   (k,j)*(<ij|kl>)
          !G_b(i,l) -= D_b   (k,j)*(<ij|kl>)

          c0 = (scf_density_matrix_ao_alpha_complex(l,j)+scf_density_matrix_ao_beta_complex(l,j)) * integral

          ao_two_e_integral_alpha_tmp(i,k) += c0
          ao_two_e_integral_beta_tmp (i,k) += c0

          ao_two_e_integral_alpha_tmp(i,l) -= SCF_density_matrix_ao_alpha_complex(k,j) * integral
          ao_two_e_integral_beta_tmp (i,l) -= scf_density_matrix_ao_beta_complex (k,j) * integral
        enddo
      else ! real part
        do k2=1,4
          if (ii(k2)==0) then
            cycle
          endif
          i = ii(k2)
          j = jj(k2)
          k = kk(k2)
          l = ll(k2)
          integral = values(k1)

          c0 = (scf_density_matrix_ao_alpha_complex(l,j)+scf_density_matrix_ao_beta_complex(l,j)) * integral

          ao_two_e_integral_alpha_tmp(i,k) += c0
          ao_two_e_integral_beta_tmp (i,k) += c0

          ao_two_e_integral_alpha_tmp(i,l) -= SCF_density_matrix_ao_alpha_complex(k,j) * integral
          ao_two_e_integral_beta_tmp (i,l) -= scf_density_matrix_ao_beta_complex (k,j) * integral
        enddo
      endif
    enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  ao_two_e_integral_alpha_complex += ao_two_e_integral_alpha_tmp
  ao_two_e_integral_beta_complex  += ao_two_e_integral_beta_tmp
  !$OMP END CRITICAL
  deallocate(keys,values,ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp)
  !$OMP END PARALLEL

 
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
      !$OMP  n_elements,ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp, &
      !$OMP  c0,key1)&
      !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha_complex, &
      !$OMP   SCF_density_matrix_ao_beta_complex,i_sign, &
      !$OMP  ao_integrals_map_2, ao_two_e_integral_alpha_complex, ao_two_e_integral_beta_complex)
 
  call get_cache_map_n_elements_max(ao_integrals_map_2,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
  allocate(ao_two_e_integral_alpha_tmp(ao_num,ao_num), &
           ao_two_e_integral_beta_tmp(ao_num,ao_num))
  ao_two_e_integral_alpha_tmp = (0.d0,0.d0)
  ao_two_e_integral_beta_tmp  = (0.d0,0.d0)
 
  !$OMP DO SCHEDULE(static,1)
  do i8=0_8,ao_integrals_map_2%map_size
    n_elements = n_elements_max
    call get_cache_map(ao_integrals_map_2,i8,keys,values,n_elements)
    do k1=1,n_elements
      ! get original key
      ! reverse of 2*key (imag part) and 2*key-1 (real part)
      key1 = shiftr(keys(k1)+1,1)

      call two_e_integrals_index_reverse_complex_2(ii,jj,kk,ll,key1)
      ! i>=k, j<=l, ik<=jl
      ! ijkl, jilk, klij*, lkji*
      if (shiftl(key1,1)==keys(k1)) then !imaginary part
        do k2=1,4
          if (ii(k2)==0) then
            cycle
          endif
          i = ii(k2)
          j = jj(k2)
          k = kk(k2)
          l = ll(k2)
          integral = i_sign(k2)*values(k1) ! for klij and lkji, take conjugate

          !G_a(i,k) += D_{ab}(l,j)*(<ij|kl>)
          !G_b(i,k) += D_{ab}(l,j)*(<ij|kl>)
          !G_a(i,l) -= D_a   (k,j)*(<ij|kl>)
          !G_b(i,l) -= D_b   (k,j)*(<ij|kl>)

          c0 = (scf_density_matrix_ao_alpha_complex(l,j)+scf_density_matrix_ao_beta_complex(l,j)) * integral

          ao_two_e_integral_alpha_tmp(i,k) += c0
          ao_two_e_integral_beta_tmp (i,k) += c0

          ao_two_e_integral_alpha_tmp(i,l) -= SCF_density_matrix_ao_alpha_complex(k,j) * integral
          ao_two_e_integral_beta_tmp (i,l) -= scf_density_matrix_ao_beta_complex (k,j) * integral
        enddo
      else ! real part
        do k2=1,4
          if (ii(k2)==0) then
            cycle
          endif
          i = ii(k2)
          j = jj(k2)
          k = kk(k2)
          l = ll(k2)
          integral = values(k1)

          c0 = (scf_density_matrix_ao_alpha_complex(l,j)+scf_density_matrix_ao_beta_complex(l,j)) * integral

          ao_two_e_integral_alpha_tmp(i,k) += c0
          ao_two_e_integral_beta_tmp (i,k) += c0

          ao_two_e_integral_alpha_tmp(i,l) -= SCF_density_matrix_ao_alpha_complex(k,j) * integral
          ao_two_e_integral_beta_tmp (i,l) -= scf_density_matrix_ao_beta_complex (k,j) * integral
        enddo
      endif
    enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  ao_two_e_integral_alpha_complex += ao_two_e_integral_alpha_tmp
  ao_two_e_integral_beta_complex  += ao_two_e_integral_beta_tmp
  !$OMP END CRITICAL
  deallocate(keys,values,ao_two_e_integral_alpha_tmp,ao_two_e_integral_beta_tmp)
  !$OMP END PARALLEL


END_PROVIDER

 BEGIN_PROVIDER [ complex*16, Fock_matrix_ao_alpha_complex, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ complex*16, Fock_matrix_ao_beta_complex,  (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Alpha Fock matrix in AO basis set
 END_DOC

 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     Fock_matrix_ao_alpha_complex(i,j) = ao_one_e_integrals_complex(i,j) + ao_two_e_integral_alpha_complex(i,j)
     Fock_matrix_ao_beta_complex (i,j) = ao_one_e_integrals_complex(i,j) + ao_two_e_integral_beta_complex (i,j)
   enddo
 enddo

END_PROVIDER
