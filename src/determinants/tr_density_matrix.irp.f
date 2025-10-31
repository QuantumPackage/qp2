BEGIN_PROVIDER [double precision, one_e_tr_dm_mo, (mo_num, mo_num, N_states, N_states)]
 use OMP_LIB
  implicit none

  BEGIN_DOC
  ! One body transition density matrix for all pairs of states n and m, < Psi^n | a_i^\dagger a_a | Psi^m >
  END_DOC

  integer                        :: j,l,m,k_a,k_b,n
  integer                        :: occ(N_int*bit_kind_size,2)
  double precision               :: ck, ckl
  double precision               :: phase
  integer                        :: h1,h2,p1,p2,degree
  integer(bit_kind)              :: tmp_det(N_int,2), tmp_det2(N_int)
  integer                        :: exc(0:2,2),n_occ(2)
  double precision, allocatable  :: tmp_a(:,:,:,:), tmp_b(:,:,:,:)
  integer                        :: krow, kcol, lrow, lcol

  double precision :: mem_tot_tr_dm
  integer :: nthreads
  nthreads = OMP_GET_MAX_THREADS()
  mem_tot_tr_dm = 8.d0*mo_num*mo_num*n_states*n_states*1d-9*(nthreads+1)
  if(mem_tot_tr_dm .gt. 0.9d0 * qp_max_mem) then
   print*,'Warning ! you are providing one_e_tr_dm_mo in parallel ! '
   print*,'Using that amount of threads ',nthreads
   print*,'Each thread needs that amount of Gb ',8.d0*mo_num*mo_num*n_states*n_states*1d-9
   print*,'So it takes in total that amount of Gb ',mem_tot_tr_dm
   print*,'The max memory for the calculation is ',qp_max_mem
   print*,'It may crash because of lack of memory ...'
  endif

  PROVIDE psi_det_alpha_unique psi_det_beta_unique

  one_e_tr_dm_mo = 0d0
 
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(j,k_a,k_b,l,m,n,occ,ck, ckl,phase,h1,h2,p1,p2,degree,exc,&
      !$OMP  tmp_a, tmp_b, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2)&
      !$OMP SHARED(N_int,N_states,elec_alpha_num,  &
      !$OMP  elec_beta_num,one_e_tr_dm_mo,N_det,&
      !$OMP  mo_num,psi_bilinear_matrix_rows,psi_bilinear_matrix_columns,&
      !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns,&
      !$OMP  psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP  psi_bilinear_matrix_values, psi_bilinear_matrix_transp_values,&
      !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here)
  allocate(tmp_a(mo_num,mo_num,N_states,N_states), tmp_b(mo_num,mo_num,N_states,N_states) )
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
      do n = 1, N_states
        ck = psi_bilinear_matrix_values(k_a,m)*psi_bilinear_matrix_values(k_a,n)
        do l=1,elec_alpha_num
          j = occ(l,1)
          tmp_a(j,j,m,n) += ck
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
        do m=1,N_states
          do n = 1, N_states
            ckl = psi_bilinear_matrix_values(k_a,m)*psi_bilinear_matrix_values(l,n) * phase
            tmp_a(h1,p1,m,n) += ckl
            ckl = psi_bilinear_matrix_values(k_a,n)*psi_bilinear_matrix_values(l,m) * phase
            tmp_a(p1,h1,m,n) += ckl
          enddo
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
  one_e_tr_dm_mo(:,:,:,:) = one_e_tr_dm_mo(:,:,:,:) + tmp_a(:,:,:,:)
  !$OMP END CRITICAL
  deallocate(tmp_a)
  !$OMP BARRIER

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
      do n = 1, N_states
        ck = psi_bilinear_matrix_transp_values(k_b,m)*psi_bilinear_matrix_transp_values(k_b,n)
        do l=1,elec_beta_num
          j = occ(l,2)
          tmp_b(j,j,m,n) += ck
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
        do m=1,N_states
          do n = 1, N_states
            ckl = psi_bilinear_matrix_transp_values(k_b,m)*psi_bilinear_matrix_transp_values(l,n) * phase
            tmp_b(h1,p1,m,n) += ckl
            ckl = psi_bilinear_matrix_transp_values(k_b,n)*psi_bilinear_matrix_transp_values(l,m) * phase
            tmp_b(p1,h1,m,n) += ckl
          enddo
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
  one_e_tr_dm_mo(:,:,:,:)  = one_e_tr_dm_mo(:,:,:,:)  + tmp_b(:,:,:,:)
  !$OMP END CRITICAL

  deallocate(tmp_b)
  !$OMP END PARALLEL

END_PROVIDER
 BEGIN_PROVIDER [ double precision, one_e_tr_dm_mo_alpha, (mo_num,mo_num,N_states,N_states) ]
&BEGIN_PROVIDER [ double precision, one_e_tr_dm_mo_beta, (mo_num,mo_num,N_states,N_states) ]
  implicit none
  BEGIN_DOC
  ! $\alpha$ and $\beta$ one-body transition density matrices for all pairs of states
  END_DOC

  integer                        :: j,l,m,n,k_a,k_b
  integer                        :: occ(N_int*bit_kind_size,2)
  double precision               :: ck, ckl
  double precision               :: phase
  integer                        :: h1,h2,p1,p2,degree
  integer(bit_kind)              :: tmp_det(N_int,2), tmp_det2(N_int)
  integer                        :: exc(0:2,2),n_occ(2)
  double precision, allocatable  :: tmp_a(:,:,:,:), tmp_b(:,:,:,:)
  integer                        :: krow, kcol, lrow, lcol

  PROVIDE psi_det_alpha_unique psi_det_beta_unique

  one_e_tr_dm_mo_alpha = 0.d0
  one_e_tr_dm_mo_beta  = 0.d0
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(j,k_a,k_b,l,m,n,occ,ck, ckl,phase,h1,h2,p1,p2,degree,exc,&
      !$OMP  tmp_a, tmp_b, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2)&
      !$OMP SHARED(N_int,N_states,elec_alpha_num,  &
      !$OMP  elec_beta_num,one_e_tr_dm_mo_alpha,one_e_tr_dm_mo_beta,N_det,&
      !$OMP  mo_num,psi_bilinear_matrix_rows,psi_bilinear_matrix_columns,&
      !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns,&
      !$OMP  psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP  psi_bilinear_matrix_values, psi_bilinear_matrix_transp_values,&
      !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here)
  allocate(tmp_a(mo_num,mo_num,N_states,N_states), tmp_b(mo_num,mo_num,N_states,N_states) )
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
      do n = 1, N_states
        ck = psi_bilinear_matrix_values(k_a,m)*psi_bilinear_matrix_values(k_a,n)
        do l=1,elec_alpha_num
          j = occ(l,1)
          tmp_a(j,j,m,n) += ck
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
        do m=1,N_states
          do n = 1, N_states
            ckl = psi_bilinear_matrix_values(k_a,m)*psi_bilinear_matrix_values(l,n) * phase
            tmp_a(h1,p1,m,n) += ckl
            tmp_a(p1,h1,m,n) += ckl
          enddo
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
  one_e_tr_dm_mo_alpha(:,:,:,:) = one_e_tr_dm_mo_alpha(:,:,:,:) + tmp_a(:,:,:,:)
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
      do n = 1, N_states
        ck = psi_bilinear_matrix_transp_values(k_b,m)*psi_bilinear_matrix_transp_values(k_b,n)
        do l=1,elec_beta_num
          j = occ(l,2)
          tmp_b(j,j,m,n) += ck
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
        do m=1,N_states
          do n = 1, N_states
            ckl = psi_bilinear_matrix_transp_values(k_b,m)*psi_bilinear_matrix_transp_values(l,n) * phase
            tmp_b(h1,p1,m,n) += ckl
            tmp_b(p1,h1,m,n) += ckl
          enddo
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
  one_e_tr_dm_mo_beta(:,:,:,:)  = one_e_tr_dm_mo_beta(:,:,:,:)  + tmp_b(:,:,:,:)
  !$OMP END CRITICAL

  deallocate(tmp_b)
  !$OMP END PARALLEL

END_PROVIDER

