 BEGIN_PROVIDER [ complex*16, one_e_dm_mo_alpha_average_complex, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ complex*16, one_e_dm_mo_beta_average_complex, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! $\alpha$ and $\beta$ one-body density matrix for each state
   END_DOC
   integer                        :: i
   one_e_dm_mo_alpha_average_complex = (0.d0,0.d0)
   one_e_dm_mo_beta_average_complex  = (0.d0,0.d0)
   do i = 1,N_states
     one_e_dm_mo_alpha_average_complex(:,:) += one_e_dm_mo_alpha_complex(:,:,i) * state_average_weight(i)
     one_e_dm_mo_beta_average_complex(:,:) += one_e_dm_mo_beta_complex(:,:,i) * state_average_weight(i)
   enddo
END_PROVIDER

BEGIN_PROVIDER [ complex*16, one_e_dm_mo_diff_complex, (mo_num,mo_num,2:N_states) ]
  implicit none
  BEGIN_DOC
  ! Difference of the one-body density matrix with respect to the ground state
  END_DOC
  integer                        :: i,j, istate

  do istate=2,N_states
    do j=1,mo_num
      do i=1,mo_num
        one_e_dm_mo_diff_complex(i,j,istate) =                            &
            one_e_dm_mo_alpha_complex(i,j,istate) - one_e_dm_mo_alpha_complex(i,j,1) +&
            one_e_dm_mo_beta_complex (i,j,istate) - one_e_dm_mo_beta_complex (i,j,1)
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ complex*16, one_e_dm_mo_spin_index_complex, (mo_num,mo_num,N_states,2) ]
  implicit none
  integer                        :: i,j,ispin,istate
  ispin = 1
  do istate = 1, N_states
    do j = 1, mo_num
      do i = 1, mo_num
        one_e_dm_mo_spin_index_complex(i,j,istate,ispin) = one_e_dm_mo_alpha_complex(i,j,istate)
      enddo
    enddo
  enddo

  ispin = 2
  do istate = 1, N_states
    do j = 1, mo_num
      do i = 1, mo_num
        one_e_dm_mo_spin_index_complex(i,j,istate,ispin) = one_e_dm_mo_beta_complex(i,j,istate)
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ complex*16, one_e_dm_dagger_mo_spin_index_complex, (mo_num,mo_num,N_states,2) ]
  print*,irp_here,' not implemented for complex'
  stop -1
!   implicit none
!   integer                        :: i,j,ispin,istate
!   ispin = 1
!   do istate = 1, N_states
!     do j = 1, mo_num
!       one_e_dm_dagger_mo_spin_index(j,j,istate,ispin) = 1 - one_e_dm_mo_alpha(j,j,istate)
!       do i = j+1, mo_num
!         one_e_dm_dagger_mo_spin_index(i,j,istate,ispin) = -one_e_dm_mo_alpha(i,j,istate)
!         one_e_dm_dagger_mo_spin_index(j,i,istate,ispin) = -one_e_dm_mo_alpha(i,j,istate)
!       enddo
!     enddo
!   enddo
!
!   ispin = 2
!   do istate = 1, N_states
!     do j = 1, mo_num
!       one_e_dm_dagger_mo_spin_index(j,j,istate,ispin) = 1 - one_e_dm_mo_beta(j,j,istate)
!       do i = j+1, mo_num
!         one_e_dm_dagger_mo_spin_index(i,j,istate,ispin) = -one_e_dm_mo_beta(i,j,istate)
!         one_e_dm_dagger_mo_spin_index(j,i,istate,ispin) = -one_e_dm_mo_beta(i,j,istate)
!       enddo
!     enddo
!   enddo
!
END_PROVIDER

 BEGIN_PROVIDER [ complex*16, one_e_dm_mo_alpha_complex, (mo_num,mo_num,N_states) ]
&BEGIN_PROVIDER [ complex*16, one_e_dm_mo_beta_complex, (mo_num,mo_num,N_states) ]
  implicit none
  BEGIN_DOC
  ! $\alpha$ and $\beta$ one-body density matrix for each state
  ! $\gamma_{\mu\nu} = \langle\Psi|a_{\nu}^{\dagger}a_{\mu}|\Psi\rangle$
  ! $\gamma_{\mu\nu} = \langle a_{\nu} \Psi|a_{\mu} \Psi\rangle$
  ! $\gamma_{\mu\nu} = \sum_{IJ} c^*_J c_I \langle a_{\nu} I|a_{\mu} J\rangle$
  END_DOC

  integer                        :: j,k,l,m,k_a,k_b
  integer                        :: occ(N_int*bit_kind_size,2)
  complex*16               :: ck, cl, ckl
  double precision               :: phase
  integer                        :: h1,h2,p1,p2,s1,s2, degree
  integer(bit_kind)              :: tmp_det(N_int,2), tmp_det2(N_int)
  integer                        :: exc(0:2,2),n_occ(2)
  complex*16, allocatable  :: tmp_a(:,:,:), tmp_b(:,:,:)
  integer                        :: krow, kcol, lrow, lcol

  PROVIDE psi_det psi_coef_complex

  one_e_dm_mo_alpha_complex = (0.d0,0.d0)
  one_e_dm_mo_beta_complex  = (0.d0,0.d0)
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(j,k,k_a,k_b,l,m,occ,ck, cl, ckl,phase,h1,h2,p1,p2,s1,s2, degree,exc,&
      !$OMP  tmp_a, tmp_b, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2)&
      !$OMP SHARED(psi_det,psi_coef_complex,N_int,N_states,elec_alpha_num,  &
      !$OMP  elec_beta_num,one_e_dm_mo_alpha_complex,one_e_dm_mo_beta_complex,N_det,&
      !$OMP  mo_num,psi_bilinear_matrix_rows,psi_bilinear_matrix_columns,&
      !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns,&
      !$OMP  psi_bilinear_matrix_order_reverse, psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP  psi_bilinear_matrix_values_complex, psi_bilinear_matrix_transp_values_complex,&
      !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here)
  allocate(tmp_a(mo_num,mo_num,N_states), tmp_b(mo_num,mo_num,N_states) )
  tmp_a = (0.d0,0.d0)
  !$OMP DO SCHEDULE(dynamic,64)
  do k_a=1,N_det
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

    ! Diagonal part
    ! -------------

    call bitstring_to_list_ab(tmp_det, occ, n_occ, N_int)
    do m=1,N_states
      ck = cdabs(psi_bilinear_matrix_values_complex(k_a,m)*psi_bilinear_matrix_values_complex(k_a,m))
      do l=1,elec_alpha_num
        j = occ(l,1)
        tmp_a(j,j,m) += ck
      enddo
    enddo

    if (k_a == N_det) cycle
    l = k_a+1
    lrow = psi_bilinear_matrix_rows(l)
    lcol = psi_bilinear_matrix_columns(l)
    ! Fix beta determinant, loop over alphas
    do while ( lcol == kcol )
      tmp_det2(:) = psi_det_alpha_unique(:, lrow)
      call get_excitation_degree_spin(tmp_det(1,1),tmp_det2,degree,N_int)
      if (degree == 1) then
        exc = 0
        call get_single_excitation_spin(tmp_det(1,1),tmp_det2,exc,phase,N_int)
        call decode_exc_spin(exc,h1,p1,h2,p2)
        ! h1 occ in k
        ! p1 occ in l
        do m=1,N_states
          ckl = dconjg(psi_bilinear_matrix_values_complex(k_a,m))*psi_bilinear_matrix_values_complex(l,m) * phase
          tmp_a(h1,p1,m) += dconjg(ckl)
          tmp_a(p1,h1,m) += ckl
        enddo
      endif
      l = l+1
      if (l>N_det) exit
      lrow = psi_bilinear_matrix_rows(l)
      lcol = psi_bilinear_matrix_columns(l)
    enddo

  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  one_e_dm_mo_alpha_complex(:,:,:) = one_e_dm_mo_alpha_complex(:,:,:) + tmp_a(:,:,:)
  !$OMP END CRITICAL
  deallocate(tmp_a)

  tmp_b = (0.d0,0.d0)
  !$OMP DO SCHEDULE(dynamic,64)
  do k_b=1,N_det
    krow = psi_bilinear_matrix_transp_rows(k_b)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_transp_columns(k_b)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

    ! Diagonal part
    ! -------------

    call bitstring_to_list_ab(tmp_det, occ, n_occ, N_int)
    do m=1,N_states
      ck = cdabs(psi_bilinear_matrix_transp_values_complex(k_b,m)*psi_bilinear_matrix_transp_values_complex(k_b,m))
      do l=1,elec_beta_num
        j = occ(l,2)
        tmp_b(j,j,m) += ck
      enddo
    enddo

    if (k_b == N_det) cycle
    l = k_b+1
    lrow = psi_bilinear_matrix_transp_rows(l)
    lcol = psi_bilinear_matrix_transp_columns(l)
    ! Fix beta determinant, loop over alphas
    do while ( lrow == krow )
      tmp_det2(:) = psi_det_beta_unique(:, lcol)
      call get_excitation_degree_spin(tmp_det(1,2),tmp_det2,degree,N_int)
      if (degree == 1) then
        exc = 0
        call get_single_excitation_spin(tmp_det(1,2),tmp_det2,exc,phase,N_int)
        call decode_exc_spin(exc,h1,p1,h2,p2)
        do m=1,N_states
          ckl = dconjg(psi_bilinear_matrix_transp_values_complex(k_b,m))*psi_bilinear_matrix_transp_values_complex(l,m) * phase
          tmp_b(h1,p1,m) += dconjg(ckl)
          tmp_b(p1,h1,m) += ckl
        enddo
      endif
      l = l+1
      if (l>N_det) exit
      lrow = psi_bilinear_matrix_transp_rows(l)
      lcol = psi_bilinear_matrix_transp_columns(l)
    enddo

  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  one_e_dm_mo_beta_complex(:,:,:)  = one_e_dm_mo_beta_complex(:,:,:)  + tmp_b(:,:,:)
  !$OMP END CRITICAL

  deallocate(tmp_b)
  !$OMP END PARALLEL

END_PROVIDER

BEGIN_PROVIDER [ complex*16, one_e_dm_mo_complex, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! One-body density matrix
   END_DOC
   one_e_dm_mo_complex = one_e_dm_mo_alpha_average_complex + one_e_dm_mo_beta_average_complex
END_PROVIDER

BEGIN_PROVIDER [ complex*16, one_e_spin_density_mo_complex, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! $\rho(\alpha) - \rho(\beta)$
   END_DOC
   one_e_spin_density_mo_complex = one_e_dm_mo_alpha_average_complex - one_e_dm_mo_beta_average_complex
END_PROVIDER


BEGIN_PROVIDER [ complex*16, one_e_spin_density_ao_complex, (ao_num,ao_num) ]
   BEGIN_DOC
   ! One body spin density matrix on the |AO| basis : $\rho_{AO}(\alpha) - \rho_{AO}(\beta)$
   ! todo: verify that this is correct for complex
   ! equivalent to using mo_to_ao_no_overlap?
   END_DOC
   implicit none
   integer                        :: i,j,k,l
   complex*16               :: dm_mo

   one_e_spin_density_ao_complex = (0.d0,0.d0)
   do k = 1, ao_num
     do l = 1, ao_num
       do i = 1, mo_num
         do j = 1, mo_num
           dm_mo = one_e_spin_density_mo_complex(j,i)
           !    if(dabs(dm_mo).le.1.d-10)cycle
           one_e_spin_density_ao_complex(l,k) += dconjg(mo_coef_complex(k,i)) * mo_coef_complex(l,j) * dm_mo

         enddo
       enddo
     enddo
   enddo

END_PROVIDER

 BEGIN_PROVIDER [ complex*16, one_e_dm_ao_alpha_complex, (ao_num,ao_num) ]
&BEGIN_PROVIDER [ complex*16, one_e_dm_ao_beta_complex, (ao_num,ao_num) ]
   BEGIN_DOC
   ! One body density matrix on the |AO| basis : $\rho_{AO}(\alpha), \rho_{AO}(\beta)$.
   END_DOC
   implicit none
   integer                        :: i,j,k,l
   complex*16               :: mo_alpha,mo_beta

   one_e_dm_ao_alpha_complex = (0.d0,0.d0)
   one_e_dm_ao_beta_complex = (0.d0,0.d0)
   do k = 1, ao_num
     do l = 1, ao_num
       do i = 1, mo_num
         do j = 1, mo_num
           mo_alpha = one_e_dm_mo_alpha_average_complex(j,i)
           mo_beta = one_e_dm_mo_beta_average_complex(j,i)
           !    if(dabs(dm_mo).le.1.d-10)cycle
           one_e_dm_ao_alpha_complex(l,k) += dconjg(mo_coef_complex(k,i)) * mo_coef_complex(l,j) *  mo_alpha
           one_e_dm_ao_beta_complex(l,k) += dconjg(mo_coef_complex(k,i)) * mo_coef_complex(l,j)  *  mo_beta
         enddo
       enddo
     enddo
   enddo

END_PROVIDER


!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!

 BEGIN_PROVIDER [ complex*16, one_e_dm_mo_alpha_average_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num) ]
&BEGIN_PROVIDER [ complex*16, one_e_dm_mo_beta_average_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! $\alpha$ and $\beta$ one-body density matrix for each state
   END_DOC
   integer                        :: i,k
   one_e_dm_mo_alpha_average_kpts = (0.d0,0.d0)
   one_e_dm_mo_beta_average_kpts  = (0.d0,0.d0)
   do i = 1,N_states
     do k=1,kpt_num
       one_e_dm_mo_alpha_average_kpts(:,:,k) += one_e_dm_mo_alpha_kpts(:,:,k,i) * state_average_weight(i)
       one_e_dm_mo_beta_average_kpts(:,:,k) += one_e_dm_mo_beta_kpts(:,:,k,i) * state_average_weight(i)
     enddo
   enddo
END_PROVIDER

BEGIN_PROVIDER [ complex*16, one_e_dm_mo_diff_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num,2:N_states) ]
  implicit none
  BEGIN_DOC
  ! Difference of the one-body density matrix with respect to the ground state
  END_DOC
  integer                        :: i,j, istate,k

  do istate=2,N_states
    do k=1,kpt_num
      do j=1,mo_num_per_kpt
        do i=1,mo_num_per_kpt
          one_e_dm_mo_diff_kpts(i,j,k,istate) =                            &
              one_e_dm_mo_alpha_kpts(i,j,k,istate) - one_e_dm_mo_alpha_kpts(i,j,k,1) +&
              one_e_dm_mo_beta_kpts (i,j,k,istate) - one_e_dm_mo_beta_kpts (i,j,k,1)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ complex*16, one_e_dm_mo_spin_index_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num,N_states,2) ]
  implicit none
  integer                        :: i,j,k,ispin,istate
  ispin = 1
  do istate = 1, N_states
    do k=1,kpt_num
      do j = 1, mo_num_per_kpt
        do i = 1, mo_num_per_kpt
          one_e_dm_mo_spin_index_kpts(i,j,k,istate,ispin) = one_e_dm_mo_alpha_kpts(i,j,k,istate)
        enddo
      enddo
    enddo
  enddo

  ispin = 2
  do istate = 1, N_states
    do k=1,kpt_num
      do j = 1, mo_num_per_kpt
        do i = 1, mo_num_per_kpt
          one_e_dm_mo_spin_index_kpts(i,j,k,istate,ispin) = one_e_dm_mo_beta_kpts(i,j,k,istate)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ complex*16, one_e_dm_dagger_mo_spin_index_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num,N_states,2) ]
  print*,irp_here,' not implemented for kpts'
  stop -1
!   implicit none
!   integer                        :: i,j,ispin,istate
!   ispin = 1
!   do istate = 1, N_states
!     do j = 1, mo_num
!       one_e_dm_dagger_mo_spin_index(j,j,istate,ispin) = 1 - one_e_dm_mo_alpha(j,j,istate)
!       do i = j+1, mo_num
!         one_e_dm_dagger_mo_spin_index(i,j,istate,ispin) = -one_e_dm_mo_alpha(i,j,istate)
!         one_e_dm_dagger_mo_spin_index(j,i,istate,ispin) = -one_e_dm_mo_alpha(i,j,istate)
!       enddo
!     enddo
!   enddo
!
!   ispin = 2
!   do istate = 1, N_states
!     do j = 1, mo_num
!       one_e_dm_dagger_mo_spin_index(j,j,istate,ispin) = 1 - one_e_dm_mo_beta(j,j,istate)
!       do i = j+1, mo_num
!         one_e_dm_dagger_mo_spin_index(i,j,istate,ispin) = -one_e_dm_mo_beta(i,j,istate)
!         one_e_dm_dagger_mo_spin_index(j,i,istate,ispin) = -one_e_dm_mo_beta(i,j,istate)
!       enddo
!     enddo
!   enddo
!
END_PROVIDER

 BEGIN_PROVIDER [ complex*16, one_e_dm_mo_alpha_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num,N_states) ]
&BEGIN_PROVIDER [ complex*16, one_e_dm_mo_beta_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num,N_states) ]
  implicit none
  BEGIN_DOC
  ! $\alpha$ and $\beta$ one-body density matrix for each state
  ! $\gamma_{\mu\nu} = \langle\Psi|a_{\nu}^{\dagger}a_{\mu}|\Psi\rangle$
  ! $\gamma_{\mu\nu} = \langle a_{\nu} \Psi|a_{\mu} \Psi\rangle$
  ! $\gamma_{\mu\nu} = \sum_{IJ} c^*_J c_I \langle a_{\nu} I|a_{\mu} J\rangle$
  END_DOC
  !todo: implement for kpts
  integer                        :: j,k,l,m,k_a,k_b
  integer                        :: occ(N_int*bit_kind_size,2)
  complex*16               :: ck, cl, ckl
  double precision               :: phase
  integer                        :: h1,h2,p1,p2,s1,s2, degree
  integer                        :: ih1,ip1,kh1,kp1,kk,k_shft,ii
  integer(bit_kind)              :: tmp_det(N_int,2), tmp_det2(N_int)
  integer(bit_kind)              :: tmp_det_kpts(N_int,2)
  integer                        :: exc(0:2,2),n_occ(2)
  complex*16, allocatable  :: tmp_a(:,:,:,:), tmp_b(:,:,:,:)
  integer                        :: krow, kcol, lrow, lcol

  PROVIDE psi_det psi_coef_complex

  one_e_dm_mo_alpha_kpts = (0.d0,0.d0)
  one_e_dm_mo_beta_kpts  = (0.d0,0.d0)
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(j,k,k_a,k_b,l,m,occ,ck, cl, ckl,phase,h1,h2,p1,p2,s1,s2, degree,exc,&
      !$OMP  tmp_a, tmp_b, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2,ih1,ip1,kh1,kp1,kk,&
      !$OMP  tmp_det_kpts,k_shft,ii)&
      !$OMP SHARED(psi_det,psi_coef_complex,N_int,N_states,elec_alpha_num_kpts,  &
      !$OMP  elec_beta_num_kpts,one_e_dm_mo_alpha_kpts,one_e_dm_mo_beta_kpts,N_det,&
      !$OMP  mo_num_per_kpt,psi_bilinear_matrix_rows,psi_bilinear_matrix_columns,&
      !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns,&
      !$OMP  psi_bilinear_matrix_order_reverse, psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP  psi_bilinear_matrix_values_complex, psi_bilinear_matrix_transp_values_complex,&
      !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here,kpt_num,kpts_bitmask)
  allocate(tmp_a(mo_num_per_kpt,mo_num_per_kpt,kpt_num,N_states), tmp_b(mo_num_per_kpt,mo_num_per_kpt,kpt_num,N_states) )
  tmp_a = (0.d0,0.d0)
  !$OMP DO SCHEDULE(dynamic,64)
  do k_a=1,N_det
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

    ! Diagonal part
    ! -------------

    do kk=1,kpt_num
      k_shft = (kk-1)*mo_num_per_kpt
      do ii=1,N_int
        tmp_det_kpts(ii,1) = iand(tmp_det(ii,1),kpts_bitmask(ii,kk))
        tmp_det_kpts(ii,2) = iand(tmp_det(ii,2),kpts_bitmask(ii,kk))
      enddo
      call bitstring_to_list_ab(tmp_det_kpts, occ, n_occ, N_int)
      do m=1,N_states
        ck = cdabs(psi_bilinear_matrix_values_complex(k_a,m)*psi_bilinear_matrix_values_complex(k_a,m))
        !do l=1,elec_alpha_num_kpts(kk)
        do l=1,n_occ(1)
          j = occ(l,1) - k_shft
          tmp_a(j,j,kk,m) += ck
        enddo
      enddo
    enddo

    if (k_a == N_det) cycle
    l = k_a+1
    lrow = psi_bilinear_matrix_rows(l)
    lcol = psi_bilinear_matrix_columns(l)
    ! Fix beta determinant, loop over alphas
    do while ( lcol == kcol )
      tmp_det2(:) = psi_det_alpha_unique(:, lrow)
      call get_excitation_degree_spin(tmp_det(1,1),tmp_det2,degree,N_int)
      if (degree == 1) then
        exc = 0
        call get_single_excitation_spin(tmp_det(1,1),tmp_det2,exc,phase,N_int)
        call decode_exc_spin(exc,h1,p1,h2,p2)
        ! h1 occ in k
        ! p1 occ in l
        ih1 = mod(h1-1,mo_num_per_kpt)+1
        ip1 = mod(p1-1,mo_num_per_kpt)+1
        kh1 = (h1-1)/mo_num_per_kpt + 1
        kp1 = (p1-1)/mo_num_per_kpt + 1
        if (kh1.ne.kp1) then
          print *,'problem in: ',irp_here,'a'
          print *,' h1 = ',h1
          print *,' p1 = ',p1
          print *,'ih1 = ',ih1
          print *,'ip1 = ',ip1
          print *,'kh1 = ',kh1
          print *,'kp1 = ',kp1
          !call debug_det(tmp_det,N_int)
          !call debug_spindet(tmp_det2,N_int)
          !call print_spindet(tmp_det2,N_int)
          !stop -2
        endif
        do m=1,N_states
          ckl = dconjg(psi_bilinear_matrix_values_complex(k_a,m))*psi_bilinear_matrix_values_complex(l,m) * phase
          tmp_a(ih1,ip1,kh1,m) += dconjg(ckl)
          tmp_a(ip1,ih1,kh1,m) += ckl
        enddo
      endif
      l = l+1
      if (l>N_det) exit
      lrow = psi_bilinear_matrix_rows(l)
      lcol = psi_bilinear_matrix_columns(l)
    enddo

  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  one_e_dm_mo_alpha_kpts(:,:,:,:) = one_e_dm_mo_alpha_kpts(:,:,:,:) + tmp_a(:,:,:,:)
  !$OMP END CRITICAL
  deallocate(tmp_a)

  tmp_b = (0.d0,0.d0)
  !$OMP DO SCHEDULE(dynamic,64)
  do k_b=1,N_det
    krow = psi_bilinear_matrix_transp_rows(k_b)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_transp_columns(k_b)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

    ! Diagonal part
    ! -------------

    do kk=1,kpt_num
      k_shft = (kk-1)*mo_num_per_kpt
      do ii=1,N_int
        tmp_det_kpts(ii,1) = iand(tmp_det(ii,1),kpts_bitmask(ii,kk))
        tmp_det_kpts(ii,2) = iand(tmp_det(ii,2),kpts_bitmask(ii,kk))
      enddo
      call bitstring_to_list_ab(tmp_det_kpts, occ, n_occ, N_int)
      do m=1,N_states
        ck = cdabs(psi_bilinear_matrix_transp_values_complex(k_b,m)*psi_bilinear_matrix_transp_values_complex(k_b,m))
        do l=1,n_occ(2)
          j = occ(l,2) - k_shft
          tmp_b(j,j,kk,m) += ck
        enddo
      enddo
    enddo

    if (k_b == N_det) cycle
    l = k_b+1
    lrow = psi_bilinear_matrix_transp_rows(l)
    lcol = psi_bilinear_matrix_transp_columns(l)
    ! Fix beta determinant, loop over alphas
    do while ( lrow == krow )
      tmp_det2(:) = psi_det_beta_unique(:, lcol)
      call get_excitation_degree_spin(tmp_det(1,2),tmp_det2,degree,N_int)
      if (degree == 1) then
        exc = 0
        call get_single_excitation_spin(tmp_det(1,2),tmp_det2,exc,phase,N_int)
        call decode_exc_spin(exc,h1,p1,h2,p2)
        ih1 = mod(h1-1,mo_num_per_kpt)+1
        ip1 = mod(p1-1,mo_num_per_kpt)+1
        kh1 = (h1-1)/mo_num_per_kpt + 1
        kp1 = (p1-1)/mo_num_per_kpt + 1
        if (kh1.ne.kp1) then
          print *,'problem in: ',irp_here,'b'
          print *,' h1 = ',h1
          print *,' p1 = ',p1
          print *,'ih1 = ',ih1
          print *,'ip1 = ',ip1
          print *,'kh1 = ',kh1
          print *,'kp1 = ',kp1
          !stop -3
        endif
        do m=1,N_states
          ckl = dconjg(psi_bilinear_matrix_transp_values_complex(k_b,m))*psi_bilinear_matrix_transp_values_complex(l,m) * phase
          tmp_b(ih1,ip1,kh1,m) += dconjg(ckl)
          tmp_b(ip1,ih1,kh1,m) += ckl
        enddo
      endif
      l = l+1
      if (l>N_det) exit
      lrow = psi_bilinear_matrix_transp_rows(l)
      lcol = psi_bilinear_matrix_transp_columns(l)
    enddo

  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  one_e_dm_mo_beta_kpts(:,:,:,:)  = one_e_dm_mo_beta_kpts(:,:,:,:)  + tmp_b(:,:,:,:)
  !$OMP END CRITICAL

  deallocate(tmp_b)
  !$OMP END PARALLEL

END_PROVIDER

BEGIN_PROVIDER [ complex*16, one_e_dm_mo_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! One-body density matrix
   END_DOC
   one_e_dm_mo_kpts = one_e_dm_mo_alpha_average_kpts + one_e_dm_mo_beta_average_kpts
END_PROVIDER

BEGIN_PROVIDER [ complex*16, one_e_spin_density_mo_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num) ]
   implicit none
   BEGIN_DOC
   ! $\rho(\alpha) - \rho(\beta)$
   END_DOC
   one_e_spin_density_mo_kpts = one_e_dm_mo_alpha_average_kpts - one_e_dm_mo_beta_average_kpts
END_PROVIDER


BEGIN_PROVIDER [ complex*16, one_e_spin_density_ao_kpts, (ao_num_per_kpt,ao_num_per_kpt,kpt_num) ]
   BEGIN_DOC
   ! One body spin density matrix on the |AO| basis : $\rho_{AO}(\alpha) - \rho_{AO}(\beta)$
   ! todo: verify that this is correct for complex
   ! equivalent to using mo_to_ao_no_overlap?
   END_DOC
   implicit none
   integer                        :: i,j,k,l,kk
   complex*16               :: dm_mo

   one_e_spin_density_ao_kpts = (0.d0,0.d0)
   do kk=1,kpt_num
     do k = 1, ao_num_per_kpt
       do l = 1, ao_num_per_kpt
         do i = 1, mo_num_per_kpt
           do j = 1, mo_num_per_kpt
             dm_mo = one_e_spin_density_mo_kpts(j,i,kk)
             !    if(dabs(dm_mo).le.1.d-10)cycle
             one_e_spin_density_ao_kpts(l,k,kk) += dconjg(mo_coef_kpts(k,i,kk)) * mo_coef_kpts(l,j,kk) * dm_mo

           enddo
         enddo
       enddo
     enddo
   enddo

END_PROVIDER

 BEGIN_PROVIDER [ complex*16, one_e_dm_ao_alpha_kpts, (ao_num_per_kpt,ao_num_per_kpt,kpt_num) ]
&BEGIN_PROVIDER [ complex*16, one_e_dm_ao_beta_kpts, (ao_num_per_kpt,ao_num_per_kpt,kpt_num) ]
   BEGIN_DOC
   ! One body density matrix on the |AO| basis : $\rho_{AO}(\alpha), \rho_{AO}(\beta)$.
   END_DOC
   implicit none
   integer                        :: i,j,k,l,kk
   complex*16               :: mo_alpha,mo_beta

   one_e_dm_ao_alpha_kpts = (0.d0,0.d0)
   one_e_dm_ao_beta_kpts = (0.d0,0.d0)
   do kk=1,kpt_num
     do k = 1, ao_num_per_kpt
       do l = 1, ao_num_per_kpt
         do i = 1, mo_num_per_kpt
           do j = 1, mo_num_per_kpt
             mo_alpha = one_e_dm_mo_alpha_average_kpts(j,i,kk)
             mo_beta = one_e_dm_mo_beta_average_kpts(j,i,kk)
             !    if(dabs(dm_mo).le.1.d-10)cycle
             one_e_dm_ao_alpha_kpts(l,k,kk) += dconjg(mo_coef_kpts(k,i,kk)) * mo_coef_kpts(l,j,kk) *  mo_alpha
             one_e_dm_ao_beta_kpts(l,k,kk) += dconjg(mo_coef_kpts(k,i,kk)) * mo_coef_kpts(l,j,kk)  *  mo_beta
           enddo
         enddo
       enddo
     enddo
   enddo

END_PROVIDER


