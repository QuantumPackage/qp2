 BEGIN_PROVIDER [ double precision, one_e_dm_mo_alpha_average, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, one_e_dm_mo_beta_average, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! $\alpha$ and $\beta$ one-body density matrix for each state
   END_DOC
   integer                        :: i
   one_e_dm_mo_alpha_average = 0.d0
   one_e_dm_mo_beta_average = 0.d0
   do i = 1,N_states
     one_e_dm_mo_alpha_average(:,:) += one_e_dm_mo_alpha(:,:,i) * state_average_weight(i)
     one_e_dm_mo_beta_average(:,:) += one_e_dm_mo_beta(:,:,i) * state_average_weight(i)
   enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, one_e_dm_mo_diff, (mo_num,mo_num,2:N_states) ]
  implicit none
  BEGIN_DOC
  ! Difference of the one-body density matrix with respect to the ground state
  END_DOC
  integer                        :: i,j, istate

  do istate=2,N_states
    do j=1,mo_num
      do i=1,mo_num
        one_e_dm_mo_diff(i,j,istate) =                            &
            one_e_dm_mo_alpha(i,j,istate) - one_e_dm_mo_alpha(i,j,1) +&
            one_e_dm_mo_beta (i,j,istate) - one_e_dm_mo_beta (i,j,1)
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, one_e_dm_mo_spin_index, (mo_num,mo_num,N_states,2) ]
  implicit none
  integer                        :: i,j,ispin,istate
  ispin = 1
  do istate = 1, N_states
    do j = 1, mo_num
      do i = 1, mo_num
        one_e_dm_mo_spin_index(i,j,istate,ispin) = one_e_dm_mo_alpha(i,j,istate)
      enddo
    enddo
  enddo

  ispin = 2
  do istate = 1, N_states
    do j = 1, mo_num
      do i = 1, mo_num
        one_e_dm_mo_spin_index(i,j,istate,ispin) = one_e_dm_mo_beta(i,j,istate)
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, one_e_dm_dagger_mo_spin_index, (mo_num,mo_num,N_states,2) ]
   implicit none
   integer                        :: i,j,ispin,istate
   ispin = 1
   do istate = 1, N_states
     do j = 1, mo_num
       one_e_dm_dagger_mo_spin_index(j,j,istate,ispin) = 1 - one_e_dm_mo_alpha(j,j,istate)
       do i = j+1, mo_num
         one_e_dm_dagger_mo_spin_index(i,j,istate,ispin) = -one_e_dm_mo_alpha(i,j,istate)
         one_e_dm_dagger_mo_spin_index(j,i,istate,ispin) = -one_e_dm_mo_alpha(i,j,istate)
       enddo
     enddo
   enddo

   ispin = 2
   do istate = 1, N_states
     do j = 1, mo_num
       one_e_dm_dagger_mo_spin_index(j,j,istate,ispin) = 1 - one_e_dm_mo_beta(j,j,istate)
       do i = j+1, mo_num
         one_e_dm_dagger_mo_spin_index(i,j,istate,ispin) = -one_e_dm_mo_beta(i,j,istate)
         one_e_dm_dagger_mo_spin_index(j,i,istate,ispin) = -one_e_dm_mo_beta(i,j,istate)
       enddo
     enddo
   enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, one_e_dm_mo_alpha, (mo_num,mo_num,N_states) ]
&BEGIN_PROVIDER [ double precision, one_e_dm_mo_beta, (mo_num,mo_num,N_states) ]
  implicit none
  BEGIN_DOC
  ! $\alpha$ and $\beta$ one-body density matrix for each state
  END_DOC

  integer                        :: j,k,l,m,k_a,k_b
  integer                        :: occ(N_int*bit_kind_size,2)
  double precision               :: ck, cl, ckl
  double precision               :: phase
  integer                        :: h1,h2,p1,p2,s1,s2, degree
  integer(bit_kind)              :: tmp_det(N_int,2), tmp_det2(N_int)
  integer                        :: exc(0:2,2),n_occ(2)
  double precision, allocatable  :: tmp_a(:,:,:), tmp_b(:,:,:)
  integer                        :: krow, kcol, lrow, lcol

  PROVIDE psi_det

  one_e_dm_mo_alpha = 0.d0
  one_e_dm_mo_beta  = 0.d0
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(j,k,k_a,k_b,l,m,occ,ck, cl, ckl,phase,h1,h2,p1,p2,s1,s2, degree,exc,&
      !$OMP  tmp_a, tmp_b, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2)&
      !$OMP SHARED(psi_det,psi_coef,N_int,N_states,elec_alpha_num,  &
      !$OMP  elec_beta_num,one_e_dm_mo_alpha,one_e_dm_mo_beta,N_det,&
      !$OMP  mo_num,psi_bilinear_matrix_rows,psi_bilinear_matrix_columns,&
      !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns,&
      !$OMP  psi_bilinear_matrix_order_reverse, psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP  psi_bilinear_matrix_values, psi_bilinear_matrix_transp_values,&
      !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here)
  allocate(tmp_a(mo_num,mo_num,N_states), tmp_b(mo_num,mo_num,N_states) )
  tmp_a = 0.d0
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
      ck = psi_bilinear_matrix_values(k_a,m)*psi_bilinear_matrix_values(k_a,m)
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
        do m=1,N_states
          ckl = psi_bilinear_matrix_values(k_a,m)*psi_bilinear_matrix_values(l,m) * phase
          tmp_a(h1,p1,m) += ckl
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
  one_e_dm_mo_alpha(:,:,:) = one_e_dm_mo_alpha(:,:,:) + tmp_a(:,:,:)
  !$OMP END CRITICAL
  deallocate(tmp_a)

  tmp_b = 0.d0
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
      ck = psi_bilinear_matrix_transp_values(k_b,m)*psi_bilinear_matrix_transp_values(k_b,m)
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
          ckl = psi_bilinear_matrix_transp_values(k_b,m)*psi_bilinear_matrix_transp_values(l,m) * phase
          tmp_b(h1,p1,m) += ckl
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
  one_e_dm_mo_beta(:,:,:)  = one_e_dm_mo_beta(:,:,:)  + tmp_b(:,:,:)
  !$OMP END CRITICAL

  deallocate(tmp_b)
  !$OMP END PARALLEL

END_PROVIDER

BEGIN_PROVIDER [ double precision, one_e_dm_mo, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! One-body density matrix
   END_DOC
   one_e_dm_mo = one_e_dm_mo_alpha_average + one_e_dm_mo_beta_average
END_PROVIDER

BEGIN_PROVIDER [ double precision, one_e_spin_density_mo, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! $\rho(\alpha) - \rho(\beta)$
   END_DOC
   one_e_spin_density_mo = one_e_dm_mo_alpha_average - one_e_dm_mo_beta_average
END_PROVIDER

subroutine set_natural_mos
   implicit none
   BEGIN_DOC
   ! Set natural orbitals, obtained by diagonalization of the one-body density matrix
   ! in the |MO| basis
   END_DOC
   character*(64)                 :: label
   double precision, allocatable  :: tmp(:,:)

   label = "Natural"
    integer :: i,j,iorb,jorb
    do i = 1, n_virt_orb
     iorb = list_virt(i)
     do j = 1, n_core_inact_act_orb
      jorb = list_core_inact_act(j)
      if(one_e_dm_mo(iorb,jorb).ne. 0.d0)then
        print*,'AHAHAH'
        print*,iorb,jorb,one_e_dm_mo(iorb,jorb)
        stop
      endif
     enddo
    enddo
   call mo_as_svd_vectors_of_mo_matrix_eig(one_e_dm_mo,size(one_e_dm_mo,1),mo_num,mo_num,mo_occ,label)
   soft_touch mo_occ

end
subroutine save_natural_mos
   implicit none
   BEGIN_DOC
   ! Save natural orbitals, obtained by diagonalization of the one-body density matrix in
   ! the |MO| basis
   END_DOC
   call set_natural_mos
   call nullify_small_elements(ao_num,mo_num,mo_coef,size(mo_coef,1),1.d-10)
   call orthonormalize_mos
   call save_mos
end


BEGIN_PROVIDER [ double precision, c0_weight, (N_states) ]
   implicit none
   BEGIN_DOC
   ! Weight of the states in the selection : $\frac{1}{c_0^2}$.
   END_DOC
   if (N_states > 1) then
     integer                        :: i
     double precision               :: c
     do i=1,N_states
       c0_weight(i) = 1.d-31
       c = maxval(psi_coef(:,i) * psi_coef(:,i))
       c0_weight(i) = 1.d0/(c+1.d-20)
     enddo
     c = 1.d0/minval(c0_weight(:))
     do i=1,N_states
       c0_weight(i) = c0_weight(i) * c
     enddo
   else
     c0_weight(:) = 1.d0
   endif

END_PROVIDER


BEGIN_PROVIDER [ double precision, state_average_weight, (N_states) ]
   implicit none
   BEGIN_DOC
   ! Weights in the state-average calculation of the density matrix
   END_DOC
   logical                        :: exists

   state_average_weight(:) = 1.d0
   if (weight_one_e_dm == 0) then
     state_average_weight(:) = c0_weight(:)
   else if (weight_one_e_dm == 1) then
     state_average_weight(:) = 1.d0/N_states
   else
     call ezfio_has_determinants_state_average_weight(exists)
     if (exists) then
       call ezfio_get_determinants_state_average_weight(state_average_weight)
     endif
   endif
   state_average_weight(:) = state_average_weight(:)+1.d-31
   state_average_weight(:) = state_average_weight(:)/(sum(state_average_weight(:)))
END_PROVIDER


BEGIN_PROVIDER [ double precision, one_e_spin_density_ao, (ao_num,ao_num) ]
   BEGIN_DOC
   ! One body spin density matrix on the |AO| basis : $\rho_{AO}(\alpha) - \rho_{AO}(\beta)$
   END_DOC
   implicit none
   integer                        :: i,j,k,l
   double precision               :: dm_mo

   one_e_spin_density_ao = 0.d0
   do k = 1, ao_num
     do l = 1, ao_num
       do i = 1, mo_num
         do j = 1, mo_num
           dm_mo = one_e_spin_density_mo(j,i)
           !    if(dabs(dm_mo).le.1.d-10)cycle
           one_e_spin_density_ao(l,k) += mo_coef(k,i) * mo_coef(l,j) * dm_mo

         enddo
       enddo
     enddo
   enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, one_e_dm_ao_alpha, (ao_num,ao_num) ]
&BEGIN_PROVIDER [ double precision, one_e_dm_ao_beta, (ao_num,ao_num) ]
   BEGIN_DOC
   ! One body density matrix on the |AO| basis : $\rho_{AO}(\alpha), \rho_{AO}(\beta)$.
   END_DOC
   implicit none
   integer                        :: i,j,k,l
   double precision               :: mo_alpha,mo_beta

   one_e_dm_ao_alpha = 0.d0
   one_e_dm_ao_beta = 0.d0
   do k = 1, ao_num
     do l = 1, ao_num
       do i = 1, mo_num
         do j = 1, mo_num
           mo_alpha = one_e_dm_mo_alpha_average(j,i)
           mo_beta = one_e_dm_mo_beta_average(j,i)
           !    if(dabs(dm_mo).le.1.d-10)cycle
           one_e_dm_ao_alpha(l,k) += mo_coef(k,i) * mo_coef(l,j) *  mo_alpha
           one_e_dm_ao_beta(l,k) += mo_coef(k,i) * mo_coef(l,j)  *  mo_beta
         enddo
       enddo
     enddo
   enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, one_e_dm_ao, (ao_num, ao_num)]
 implicit none
   BEGIN_DOC
   !  one_e_dm_ao = one_e_dm_ao_alpha + one_e_dm_ao_beta 
   END_DOC
    one_e_dm_ao = one_e_dm_ao_alpha + one_e_dm_ao_beta 
END_PROVIDER 


subroutine get_occupation_from_dets(istate,occupation)
  implicit none
  double precision, intent(out)  :: occupation(mo_num)
  integer, intent(in)            :: istate
  BEGIN_DOC
  ! Returns the average occupation of the MOs
  END_DOC
  integer                        :: i,j, ispin
  integer                        :: list(N_int*bit_kind_size,2)
  integer                        :: n_elements(2)
  double precision               :: c, norm_2
  ASSERT (istate > 0)
  ASSERT (istate <= N_states)

  occupation = 0.d0
  double precision, external :: u_dot_u

  norm_2 = 1.d0/u_dot_u(psi_coef(1,istate),N_det)

  do i=1,N_det
    c = psi_coef(i,istate)*psi_coef(i,istate)*norm_2
    call bitstring_to_list_ab(psi_det(1,1,i), list, n_elements, N_int)
    do ispin=1,2
      do j=1,n_elements(ispin)
        ASSERT ( list(j,ispin) < mo_num )
        occupation( list(j,ispin) ) += c
      enddo
    enddo
  enddo
end

