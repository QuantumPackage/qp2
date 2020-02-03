program print_2e_integrals_from_map
  call run
end

subroutine run
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
 
  PROVIDE ao_two_e_integrals_in_map
 
  integer(omp_lock_kind) :: lck(ao_num)
  integer(map_size_kind)     :: i8
  integer                        :: ii(4), jj(4), kk(4), ll(4), k2
  integer(cache_map_size_kind)   :: n_elements_max, n_elements
  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  complex*16, parameter    :: i_sign(4) = (/(0.d0,1.d0),(0.d0,1.d0),(0.d0,-1.d0),(0.d0,-1.d0)/)
  integer(key_kind) :: key1
 
  call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
 
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
          print'((A),4(I4),1(E15.7),2(I),2(E9.1))','imag1  ',i,j,k,l,values(k1),k1,k2,i_sign(k2)

          !G_a(i,k) += D_{ab}(l,j)*(<ij|kl>)
          !G_b(i,k) += D_{ab}(l,j)*(<ij|kl>)
          !G_a(i,l) -= D_a   (k,j)*(<ij|kl>)
          !G_b(i,l) -= D_b   (k,j)*(<ij|kl>)

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
          print'((A),4(I4),1(E15.7),2(I))','real1  ',i,j,k,l,values(k1),k1,k2
        enddo
      endif
    enddo
  enddo
  deallocate(keys,values)

 
  call get_cache_map_n_elements_max(ao_integrals_map_2,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
 
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
          print'((A),4(I4),1(E15.7),2(I),2(E9.1))','imag2  ',i,j,k,l,values(k1),k1,k2,i_sign(k2)
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
          print'((A),4(I4),1(E15.7),2(I))','real2  ',i,j,k,l,values(k1),k1,k2
        enddo
      endif
    enddo
  enddo
  deallocate(keys,values)
end
