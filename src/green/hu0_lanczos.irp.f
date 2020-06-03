! modified from H_S2_u_0_nstates_openmp in Davidson/u0Hu0.irp.f

subroutine H_u_0_openmp(v_0,u_0,sze)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0 = H|u_0>
  !
  ! Assumes that the determinants are in psi_det
  !
  ! istart, iend, ishift, istep are used in ZMQ parallelization.
  END_DOC
  integer :: N_st=1
  integer, intent(in)            :: sze
  complex*16, intent(inout)  :: v_0(sze), u_0(sze)
  integer :: k
  call cdset_order(u_0(1),psi_bilinear_matrix_order,N_det)
  v_0 = (0.d0,0.d0)
  
  call h_u_0_openmp_work(v_0,u_0,sze,1,N_det,0,1)

  call cdset_order(v_0(1),psi_bilinear_matrix_order_reverse,N_det)
  call cdset_order(u_0(1),psi_bilinear_matrix_order_reverse,N_det)

end


subroutine H_u_0_openmp_work(v_t,u_t,sze,istart,iend,ishift,istep)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_t = H|u_t>
  !
  ! Default should be 1,N_det,0,1
  END_DOC
  integer :: N_st=1
  integer, intent(in)            :: sze,istart,iend,ishift,istep
  complex*16, intent(in)   :: u_t(N_det)
  complex*16, intent(out)  :: v_t(sze)

  
  PROVIDE ref_bitmask_energy N_int

  select case (N_int)
    case (1)
      call H_u_0_openmp_work_1(v_t,u_t,sze,istart,iend,ishift,istep)
    case (2)
      call H_u_0_openmp_work_2(v_t,u_t,sze,istart,iend,ishift,istep)
    case (3)
      call H_u_0_openmp_work_3(v_t,u_t,sze,istart,iend,ishift,istep)
    case (4)
      call H_u_0_openmp_work_4(v_t,u_t,sze,istart,iend,ishift,istep)
    case default
      call H_u_0_openmp_work_N_int(v_t,u_t,sze,istart,iend,ishift,istep)
  end select
end
BEGIN_TEMPLATE

subroutine H_u_0_openmp_work_$N_int(v_t,u_t,sze,istart,iend,ishift,istep)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_t = H|u_t>
  !
  ! Default should be 1,N_det,0,1
  END_DOC
  integer :: N_st=1
  integer, intent(in)            :: sze,istart,iend,ishift,istep
  complex*16, intent(in)   :: u_t(N_det)
  complex*16, intent(out)  :: v_t(sze) 

  complex*16               :: hij
  double precision :: hii 
  integer                        :: i,j,k,l
  integer                        :: k_a, k_b, l_a, l_b, m_a, m_b
  integer                        :: istate
  integer                        :: krow, kcol, krow_b, kcol_b
  integer                        :: lrow, lcol
  integer                        :: mrow, mcol
  integer(bit_kind)              :: spindet($N_int)
  integer(bit_kind)              :: tmp_det($N_int,2)
  integer(bit_kind)              :: tmp_det2($N_int,2)
  integer(bit_kind)              :: tmp_det3($N_int,2)
  integer(bit_kind), allocatable :: buffer(:,:)
  integer                        :: n_doubles
  integer, allocatable           :: doubles(:)
  integer, allocatable           :: singles_a(:)
  integer, allocatable           :: singles_b(:)
  integer, allocatable           :: idx(:), idx0(:)
  integer                        :: maxab, n_singles_a, n_singles_b, kcol_prev
  integer*8                      :: k8

  maxab = max(N_det_alpha_unique, N_det_beta_unique)+1
  allocate(idx0(maxab))
  
  do i=1,maxab
    idx0(i) = i
  enddo

  ! Prepare the array of all alpha single excitations
  ! -------------------------------------------------

  PROVIDE N_int nthreads_davidson
  !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads_davidson)        &
      !$OMP   SHARED(psi_bilinear_matrix_rows, N_det,                &
      !$OMP          psi_bilinear_matrix_columns,                    &
      !$OMP          psi_det_alpha_unique, psi_det_beta_unique,      &
      !$OMP          n_det_alpha_unique, n_det_beta_unique, N_int,   &
      !$OMP          psi_bilinear_matrix_transp_rows,                &
      !$OMP          psi_bilinear_matrix_transp_columns,             &
      !$OMP          psi_bilinear_matrix_transp_order, N_st,         &
      !$OMP          psi_bilinear_matrix_order_transp_reverse,       &
      !$OMP          psi_bilinear_matrix_columns_loc,                &
      !$OMP          psi_bilinear_matrix_transp_rows_loc,            &
      !$OMP          istart, iend, istep, irp_here, v_t,             &
      !$OMP          ishift, idx0, u_t, maxab)                       &
      !$OMP   PRIVATE(krow, kcol, tmp_det, spindet, k_a, k_b, i,     &
      !$OMP          lcol, lrow, l_a, l_b,                           &
      !$OMP          buffer, doubles, n_doubles,                     &
      !$OMP          tmp_det2, hii, hij, idx, l, kcol_prev,          &
      !$OMP          singles_a, n_singles_a, singles_b,              &
      !$OMP          n_singles_b, k8)
  
  ! Alpha/Beta double excitations
  ! =============================
    
  allocate( buffer($N_int,maxab),                                     &
      singles_a(maxab),                                              &
      singles_b(maxab),                                              &
      doubles(maxab),                                                &
      idx(maxab))

  kcol_prev=-1

  ASSERT (iend <= N_det)
  ASSERT (istart > 0)
  ASSERT (istep  > 0)

  !$OMP DO SCHEDULE(dynamic,64)
  do k_a=istart+ishift,iend,istep

    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)

    if (kcol /= kcol_prev) then
      call get_all_spin_singles_$N_int(                              &
          psi_det_beta_unique, idx0,                                 &
          tmp_det(1,2), N_det_beta_unique,                           &
          singles_b, n_singles_b)
    endif
    kcol_prev = kcol

    ! Loop over singly excited beta columns
    ! -------------------------------------

    do i=1,n_singles_b
      lcol = singles_b(i)

      tmp_det2(1:$N_int,2) = psi_det_beta_unique(1:$N_int, lcol)

      l_a = psi_bilinear_matrix_columns_loc(lcol)
      ASSERT (l_a <= N_det)

      do j=1,psi_bilinear_matrix_columns_loc(lcol+1) - l_a
        lrow = psi_bilinear_matrix_rows(l_a)
        ASSERT (lrow <= N_det_alpha_unique)

        buffer(1:$N_int,j) = psi_det_alpha_unique(1:$N_int, lrow)

        ASSERT (l_a <= N_det)
        idx(j) = l_a
        l_a = l_a+1
      enddo
      j = j-1

      call get_all_spin_singles_$N_int(                              &
          buffer, idx, tmp_det(1,1), j,                              &
          singles_a, n_singles_a )

      ! Loop over alpha singles
      ! -----------------------

      do k = 1,n_singles_a
        l_a = singles_a(k)
        ASSERT (l_a <= N_det)

        lrow = psi_bilinear_matrix_rows(l_a)
        ASSERT (lrow <= N_det_alpha_unique)

        tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
        call i_h_j_double_alpha_beta_complex(tmp_det,tmp_det2,$N_int,hij)
        v_t(k_a) = v_t(k_a) + hij * u_t(l_a)
      enddo
    enddo
  enddo
  !$OMP END DO 

  !$OMP DO SCHEDULE(dynamic,64)
  do k_a=istart+ishift,iend,istep


    ! Single and double alpha excitations
    ! ===================================
    
    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    
    ! Initial determinant is at k_b in beta-major representation
    ! ----------------------------------------------------------------------
    
    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
    ASSERT (k_b <= N_det)

    spindet(1:$N_int) = tmp_det(1:$N_int,1)
    
    ! Loop inside the beta column to gather all the connected alphas
    lcol = psi_bilinear_matrix_columns(k_a)
    l_a = psi_bilinear_matrix_columns_loc(lcol)
    do i=1,N_det_alpha_unique
      if (l_a > N_det) exit
      lcol = psi_bilinear_matrix_columns(l_a)
      if (lcol /= kcol) exit
      lrow = psi_bilinear_matrix_rows(l_a)
      ASSERT (lrow <= N_det_alpha_unique)

      buffer(1:$N_int,i) = psi_det_alpha_unique(1:$N_int, lrow)
      idx(i) = l_a
      l_a = l_a+1
    enddo
    i = i-1
    
    !call get_all_spin_singles_and_doubles_$N_int(                    &
    !    buffer, idx, spindet, i,                                     &
    !    singles_a, doubles, n_singles_a, n_doubles )
    call get_all_spin_singles_and_doubles(                           &
        buffer, idx, spindet, $N_int, i,                             &
        singles_a, doubles, n_singles_a, n_doubles )

    ! Compute Hij for all alpha singles
    ! ----------------------------------

    tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    do i=1,n_singles_a
      l_a = singles_a(i)
      ASSERT (l_a <= N_det)

      lrow = psi_bilinear_matrix_rows(l_a)
      ASSERT (lrow <= N_det_alpha_unique)

      tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
      call i_h_j_single_spin_complex( tmp_det, tmp_det2, $N_int, 1, hij)

      v_t(k_a) = v_t(k_a) + hij * u_t(l_a)
    enddo

    
    ! Compute Hij for all alpha doubles
    ! ----------------------------------
    
    do i=1,n_doubles
      l_a = doubles(i)
      ASSERT (l_a <= N_det)

      lrow = psi_bilinear_matrix_rows(l_a)
      ASSERT (lrow <= N_det_alpha_unique)

      call i_h_j_double_spin_complex( tmp_det(1,1), psi_det_alpha_unique(1, lrow), $N_int, hij)
      v_t(k_a) = v_t(k_a) + hij * u_t(l_a)
    enddo
    

    ! Single and double beta excitations
    ! ==================================

    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    
    spindet(1:$N_int) = tmp_det(1:$N_int,2)
    
    ! Initial determinant is at k_b in beta-major representation
    ! -----------------------------------------------------------------------

    k_b = psi_bilinear_matrix_order_transp_reverse(k_a) 
    ASSERT (k_b <= N_det)
    
    ! Loop inside the alpha row to gather all the connected betas
    lrow = psi_bilinear_matrix_transp_rows(k_b)
    l_b = psi_bilinear_matrix_transp_rows_loc(lrow)
    do i=1,N_det_beta_unique
      if (l_b > N_det) exit
      lrow = psi_bilinear_matrix_transp_rows(l_b)
      if (lrow /= krow) exit
      lcol = psi_bilinear_matrix_transp_columns(l_b)
      ASSERT (lcol <= N_det_beta_unique)

      buffer(1:$N_int,i) = psi_det_beta_unique(1:$N_int, lcol)
      idx(i) = l_b
      l_b = l_b+1
    enddo
    i = i-1
  
    !call get_all_spin_singles_and_doubles_$N_int(                    &
    !    buffer, idx, spindet, i,                                     &
    !    singles_b, doubles, n_singles_b, n_doubles )
    call get_all_spin_singles_and_doubles(                           &
        buffer, idx, spindet, $N_int, i,                             &
        singles_b, doubles, n_singles_b, n_doubles )
    
    ! Compute Hij for all beta singles
    ! ----------------------------------
    
    tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    do i=1,n_singles_b
      l_b = singles_b(i)
      ASSERT (l_b <= N_det)

      lcol = psi_bilinear_matrix_transp_columns(l_b)
      ASSERT (lcol <= N_det_beta_unique)

      tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, lcol)
      call i_h_j_single_spin_complex( tmp_det, tmp_det2, $N_int, 2, hij)
      l_a = psi_bilinear_matrix_transp_order(l_b)
      ASSERT (l_a <= N_det)
      v_t(k_a) = v_t(k_a) + hij * u_t(l_a)
    enddo
    
    ! Compute Hij for all beta doubles
    ! ----------------------------------
    
    do i=1,n_doubles
      l_b = doubles(i)
      ASSERT (l_b <= N_det)

      lcol = psi_bilinear_matrix_transp_columns(l_b)
      ASSERT (lcol <= N_det_beta_unique)

      call i_h_j_double_spin_complex( tmp_det(1,2), psi_det_beta_unique(1, lcol), $N_int, hij)
      l_a = psi_bilinear_matrix_transp_order(l_b)
      ASSERT (l_a <= N_det)

      v_t(k_a) = v_t(k_a) + hij * u_t(l_a)
    enddo


    ! Diagonal contribution
    ! =====================

    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    
    !double precision, external :: diag_H_mat_elem, diag_S_mat_elem
    double precision, external :: diag_H_mat_elem
    hii = diag_H_mat_elem(tmp_det,$N_int) 
    v_t(k_a) = v_t(k_a) + hii * u_t(k_a)

  end do
  !$OMP END DO 
  deallocate(buffer, singles_a, singles_b, doubles, idx)
  !$OMP END PARALLEL

end

SUBST [ N_int ]

1;;
2;;
3;;
4;;
N_int;;

END_TEMPLATE

