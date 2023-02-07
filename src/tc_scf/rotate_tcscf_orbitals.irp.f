
! ---

program rotate_tcscf_orbitals

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  bi_ortho = .True.
  touch bi_ortho

  call maximize_overlap()

end

! ---

subroutine maximize_overlap()

  implicit none
  integer                       :: i, m, n
  double precision              :: accu_d, accu_nd
  double precision, allocatable :: C(:,:), R(:,:), L(:,:), W(:,:), e(:)
  double precision, allocatable :: S(:,:)

  n = ao_num
  m = mo_num

  allocate(L(n,m), R(n,m), C(n,m), W(n,n), e(m))
  L = mo_l_coef
  R = mo_r_coef
  C = mo_coef
  W = ao_overlap

  print*, ' fock matrix diag elements'
  do i = 1, m
    e(i) = Fock_matrix_tc_mo_tot(i,i)
    print*, e(i)
  enddo

  ! ---
   
  print *, ' overlap before :'
  print *, ' '

  allocate(S(m,m)) 

  call LTxSxR(n, m, L, W, R, S)
  !print*, " L.T x R"
  !do i = 1, m
  !  write(*, '(100(F16.10,X))') S(i,i)
  !enddo
  call LTxSxR(n, m, L, W, C, S)
  print*, " L.T x C"
  do i = 1, m
    write(*, '(100(F16.10,X))') S(i,:)
  enddo
  call LTxSxR(n, m, C, W, R, S)
  print*, " C.T x R"
  do i = 1, m
    write(*, '(100(F16.10,X))') S(i,:)
  enddo

  deallocate(S)

  ! ---

  call rotate_degen_eigvec_to_maximize_overlap(n, m, e, C, W, L, R)

  ! ---
   
  print *, ' overlap after :'
  print *, ' '

  allocate(S(m,m)) 

  call LTxSxR(n, m, L, W, R, S)
  !print*, " L.T x R"
  !do i = 1, m
  !  write(*, '(100(F16.10,X))') S(i,i)
  !enddo
  call LTxSxR(n, m, L, W, C, S)
  print*, " L.T x C"
  do i = 1, m
    write(*, '(100(F16.10,X))') S(i,:)
  enddo
  call LTxSxR(n, m, C, W, R, S)
  print*, " C.T x R"
  do i = 1, m
    write(*, '(100(F16.10,X))') S(i,:)
  enddo

  deallocate(S)

  ! ---

  mo_l_coef = L
  mo_r_coef = R
  call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
  call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)

  ! ---

  deallocate(L, R, C, W, e)

end subroutine maximize_overlap

! ---

subroutine rotate_degen_eigvec_to_maximize_overlap(n, m, e0, C0, W0, L0, R0)

  implicit none

  integer,          intent(in)    :: n, m
  double precision, intent(in)    :: e0(m), W0(n,n), C0(n,m)
  double precision, intent(inout) :: L0(n,m), R0(n,m)


  integer                         :: i, j, k, kk, mm, id1, tot_deg
  double precision                :: ei, ej, de, de_thr
  integer,          allocatable   :: deg_num(:)
  double precision, allocatable   :: L(:,:), R(:,:), C(:,:), Lnew(:,:), Rnew(:,:), tmp(:,:)
  !double precision, allocatable   :: S(:,:), Snew(:,:), T(:,:), Ttmp(:,:), Stmp(:,:)
  double precision, allocatable   :: S(:,:), Snew(:,:), T(:,:), Ttmp(:,:), Stmp(:,:)
  !real*8                          :: S(m,m), Snew(m,m), T(m,m)

  id1 = 700
  allocate(S(id1,id1), Snew(id1,id1), T(id1,id1))

  ! ---

  allocate( deg_num(m) )
  do i = 1, m
    deg_num(i) = 1
  enddo

  de_thr = thr_degen_tc

  do i = 1, m-1
    ei = e0(i)

    ! already considered in degen vectors
    if(deg_num(i).eq.0) cycle

    do j = i+1, m
      ej = e0(j)
      de = dabs(ei - ej)

      if(de .lt. de_thr) then
        deg_num(i) = deg_num(i) + 1
        deg_num(j) = 0
      endif

    enddo
  enddo

  tot_deg = 0
  do i = 1, m
    if(deg_num(i).gt.1) then
      print *, ' degen on', i, deg_num(i)
      tot_deg = tot_deg + 1
    endif
  enddo

  if(tot_deg .eq. 0) then
    print *, ' no degen'
    return
  endif

  ! ---

  do i = 1, m
    mm = deg_num(i)

    if(mm .gt. 1) then

      allocate(L(n,mm), R(n,mm), C(n,mm))
      do j = 1, mm
        L(1:n,j) = L0(1:n,i+j-1)
        R(1:n,j) = R0(1:n,i+j-1)
        C(1:n,j) = C0(1:n,i+j-1)
      enddo

      ! ---

      ! C.T x W0 x R
      allocate(tmp(mm,n), Stmp(mm,mm))
      call dgemm( 'T', 'N', mm, n, n, 1.d0       &
                , C, size(C, 1), W0, size(W0, 1) &
                , 0.d0, tmp, size(tmp, 1) )
      call dgemm( 'N', 'N', mm, mm, n, 1.d0        &
                , tmp, size(tmp, 1), R, size(R, 1) &
                , 0.d0, Stmp, size(Stmp, 1) )
      deallocate(C, tmp)

      S = 0.d0
      do k = 1, mm
        do kk = 1, mm
          S(kk,k) = Stmp(kk,k)
        enddo
      enddo
      deallocate(Stmp)

      !print*, " overlap bef"
      !do k = 1, mm
      !  write(*, '(100(F16.10,X))') (S(k,kk), kk=1, mm)
      !enddo
    
      T    = 0.d0
      Snew = 0.d0
      call maxovl(mm, mm, S, T, Snew)

      !print*, " overlap aft"
      !do k = 1, mm
      !  write(*, '(100(F16.10,X))') (Snew(k,kk), kk=1, mm)
      !enddo

      allocate(Ttmp(mm,mm))
      Ttmp(1:mm,1:mm) = T(1:mm,1:mm)

      allocate(Lnew(n,mm), Rnew(n,mm))
      call dgemm( 'N', 'N', n, mm, mm, 1.d0               &
                , R, size(R, 1), Ttmp(1,1), size(Ttmp, 1) &
                , 0.d0, Rnew, size(Rnew, 1) )
      call dgemm( 'N', 'N', n, mm, mm, 1.d0               &
                , L, size(L, 1), Ttmp(1,1), size(Ttmp, 1) &
                , 0.d0, Lnew, size(Lnew, 1) )

      deallocate(L, R)
      deallocate(Ttmp)

      ! ---

      do j = 1, mm
        L0(1:n,i+j-1) = Lnew(1:n,j)
        R0(1:n,i+j-1) = Rnew(1:n,j)
      enddo
      deallocate(Lnew, Rnew)

    endif
  enddo

  deallocate(S, Snew, T)

end subroutine rotate_degen_eigvec_to_maximize_overlap

! ---

subroutine fix_right_to_one()

  implicit none
  integer                       :: i, j, m, n, mm, tot_deg
  double precision              :: accu_d, accu_nd
  double precision              :: de_thr, ei, ej, de
  integer,          allocatable :: deg_num(:)
  double precision, allocatable :: R0(:,:), L0(:,:), W(:,:), e0(:)
  double precision, allocatable :: R(:,:), L(:,:), S(:,:), Stmp(:,:), tmp(:,:)

  n = ao_num
  m = mo_num

  allocate(L0(n,m), R0(n,m), W(n,n), e0(m))
  L0 = mo_l_coef
  R0 = mo_r_coef
  W  = ao_overlap

  print*, ' fock matrix diag elements'
  do i = 1, m
    e0(i) = Fock_matrix_tc_mo_tot(i,i)
    print*, e0(i)
  enddo

  ! ---

  allocate( deg_num(m) )
  do i = 1, m
    deg_num(i) = 1
  enddo

  de_thr = 1d-6

  do i = 1, m-1
    ei = e0(i)

    ! already considered in degen vectors
    if(deg_num(i).eq.0) cycle

    do j = i+1, m
      ej = e0(j)
      de = dabs(ei - ej)

      if(de .lt. de_thr) then
        deg_num(i) = deg_num(i) + 1
        deg_num(j) = 0
      endif

    enddo
  enddo

  deallocate(e0)

  tot_deg = 0
  do i = 1, m
    if(deg_num(i).gt.1) then
      print *, ' degen on', i, deg_num(i)
      tot_deg = tot_deg + 1
    endif
  enddo

  if(tot_deg .eq. 0) then
    print *, ' no degen'
    return
  endif

  ! ---

  do i = 1, m
    mm = deg_num(i)

    if(mm .gt. 1) then

      allocate(L(n,mm), R(n,mm))
      do j = 1, mm
        L(1:n,j) = L0(1:n,i+j-1)
        R(1:n,j) = R0(1:n,i+j-1)
      enddo

      ! ---

      call impose_weighted_orthog_svd(n, mm, W, R)
      call impose_weighted_biorthog_qr(n, mm, thresh_biorthog_diag, thresh_biorthog_nondiag, R, W, L)

      ! ---

      do j = 1, mm
        L0(1:n,i+j-1) = L(1:n,j)
        R0(1:n,i+j-1) = R(1:n,j)
      enddo
      deallocate(L, R)

    endif
  enddo

  call check_weighted_biorthog_binormalize(n, m, L0, W, R0, thresh_biorthog_diag, thresh_biorthog_nondiag, .true.)

  deallocate(W, deg_num)

  mo_l_coef = L0
  mo_r_coef = R0
  deallocate(L0, R0)

  call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
  call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
  print *, ' orbitals are rotated '

  return
end subroutine fix_right_to_one

! ---
