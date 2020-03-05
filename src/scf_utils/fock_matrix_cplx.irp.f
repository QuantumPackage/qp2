 BEGIN_PROVIDER [ complex*16, Fock_matrix_mo_complex, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_diag_mo_complex, (mo_num)]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis.
   ! For open shells, the ROHF Fock Matrix is ::
   !
   !       |   F-K    |  F + K/2  |    F     |
   !       |---------------------------------|
   !       | F + K/2  |     F     |  F - K/2 |
   !       |---------------------------------|
   !       |    F     |  F - K/2  |  F + K   |
   !
   !
   ! F = 1/2 (Fa + Fb)
   !
   ! K = Fb - Fa
   !
   END_DOC
   integer                        :: i,j,n
   if (elec_alpha_num == elec_beta_num) then
     Fock_matrix_mo_complex = Fock_matrix_mo_alpha_complex
   else

     do j=1,elec_beta_num
       ! F-K
       do i=1,elec_beta_num !CC
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))&
             - (Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
       ! F+K/2
       do i=elec_beta_num+1,elec_alpha_num  !CA
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))&
             + 0.5d0*(Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
       ! F
       do i=elec_alpha_num+1, mo_num !CV
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))
       enddo
     enddo

     do j=elec_beta_num+1,elec_alpha_num
       ! F+K/2
       do i=1,elec_beta_num !AC
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))&
             + 0.5d0*(Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
       ! F
       do i=elec_beta_num+1,elec_alpha_num !AA
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))
       enddo
       ! F-K/2
       do i=elec_alpha_num+1, mo_num !AV
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))&
             - 0.5d0*(Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
     enddo

     do j=elec_alpha_num+1, mo_num
       ! F
       do i=1,elec_beta_num !VC
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))
       enddo
       ! F-K/2
       do i=elec_beta_num+1,elec_alpha_num !VA
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j))&
             - 0.5d0*(Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
       ! F+K
       do i=elec_alpha_num+1,mo_num !VV
         Fock_matrix_mo_complex(i,j) = 0.5d0*(Fock_matrix_mo_alpha_complex(i,j)+Fock_matrix_mo_beta_complex(i,j)) &
             + (Fock_matrix_mo_beta_complex(i,j) - Fock_matrix_mo_alpha_complex(i,j))
       enddo
     enddo

   endif

   do i = 1, mo_num
     Fock_matrix_diag_mo_complex(i) = dble(Fock_matrix_mo_complex(i,i))
     if (dabs(dimag(Fock_matrix_mo_complex(i,i))) .gt. 1.0d-12) then
       !stop 'diagonal elements of Fock matrix should be real'
       print *, 'diagonal elements of Fock matrix should be real',i,Fock_matrix_mo_complex(i,i)
       stop -1
     endif
   enddo


   if(frozen_orb_scf)then
     integer                        :: iorb,jorb
     do i = 1, n_core_orb
      iorb = list_core(i)
      do j = 1, n_act_orb
       jorb = list_act(j)
       Fock_matrix_mo_complex(iorb,jorb) = (0.d0,0.d0)
       Fock_matrix_mo_complex(jorb,iorb) = (0.d0,0.d0)
      enddo
     enddo
   endif

END_PROVIDER



BEGIN_PROVIDER [ complex*16, Fock_matrix_mo_alpha_complex, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis
   END_DOC
   call ao_to_mo_complex(Fock_matrix_ao_alpha_complex,size(Fock_matrix_ao_alpha_complex,1), &
                 Fock_matrix_mo_alpha_complex,size(Fock_matrix_mo_alpha_complex,1))
END_PROVIDER

BEGIN_PROVIDER [ complex*16, Fock_matrix_mo_beta_complex, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis
   END_DOC
   call ao_to_mo_complex(Fock_matrix_ao_beta_complex,size(Fock_matrix_ao_beta_complex,1), &
                 Fock_matrix_mo_beta_complex,size(Fock_matrix_mo_beta_complex,1))
END_PROVIDER


BEGIN_PROVIDER [ complex*16, Fock_matrix_ao_complex, (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Fock matrix in AO basis set
 END_DOC

 if(frozen_orb_scf)then
   call mo_to_ao_complex(Fock_matrix_mo_complex,size(Fock_matrix_mo_complex,1),              &
       Fock_matrix_ao_complex,size(Fock_matrix_ao_complex,1))
 else
   if ( (elec_alpha_num == elec_beta_num).and.                       &
         (level_shift == 0.) )                                       &
         then
     integer                        :: i,j
     do j=1,ao_num
       do i=1,ao_num
         Fock_matrix_ao_complex(i,j) = Fock_matrix_ao_alpha_complex(i,j)
       enddo
     enddo
   else
     call mo_to_ao_complex(Fock_matrix_mo_complex,size(Fock_matrix_mo_complex,1),            &
         Fock_matrix_ao_complex,size(Fock_matrix_ao_complex,1))
   endif
 endif
END_PROVIDER


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
  complex*16, allocatable  :: ao_two_e_integral_alpha_tmp(:,:)
  complex*16, allocatable  :: ao_two_e_integral_beta_tmp(:,:)
 
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
      !$OMP  SCF_density_matrix_ao_beta_complex, &
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
      !$OMP   SCF_density_matrix_ao_beta_complex, &
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
