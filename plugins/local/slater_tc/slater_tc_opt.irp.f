
! ---

subroutine provide_all_three_ints_bi_ortho()

  BEGIN_DOC
  ! routine that provides all necessary three-electron integrals
  END_DOC

  implicit none
  double precision :: t1, t2

  print *, ' start provide_all_three_ints_bi_ortho'
  call wall_time(t1)

  if(three_body_h_tc) then

    if(three_e_3_idx_term) then
      PROVIDE three_e_3_idx_direct_bi_ort three_e_3_idx_cycle_1_bi_ort three_e_3_idx_cycle_2_bi_ort
      PROVIDE three_e_3_idx_exch23_bi_ort three_e_3_idx_exch13_bi_ort three_e_3_idx_exch12_bi_ort
    endif

    if(three_e_4_idx_term) then
      PROVIDE three_e_4_idx_direct_bi_ort three_e_4_idx_cycle_1_bi_ort three_e_4_idx_exch23_bi_ort three_e_4_idx_exch13_bi_ort
    endif
    if(pure_three_body_h_tc)then
     provide three_body_ints_bi_ort
    endif

    if(.not. double_normal_ord .and. three_e_5_idx_term) then
      PROVIDE three_e_5_idx_direct_bi_ort
    elseif(double_normal_ord .and. (.not. three_e_5_idx_term)) then
      PROVIDE normal_two_body_bi_orth
    endif

  endif

  call wall_time(t2)
  print *, ' end provide_all_three_ints_bi_ortho after (min) = ', (t2-t1)/60.d0

  return
end

! ---

subroutine htilde_mu_mat_opt_bi_ortho_tot(key_j, key_i, Nint, htot)

  implicit none

  BEGIN_DOC
  !
  ! <key_j |H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
  !!
  ! Returns the total matrix element
  !! WARNING !!
  !
  ! Non hermitian !!
  !
  END_DOC

  use bitmasks
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out) :: htot
  double precision              :: hmono, htwoe, hthree

  call htilde_mu_mat_opt_bi_ortho(key_j, key_i, Nint, hmono, htwoe, hthree, htot)

end

! ---

subroutine htilde_mu_mat_opt_bi_ortho(key_j, key_i, Nint, hmono, htwoe, hthree, htot)

  BEGIN_DOC
  !
  ! <key_j |H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
  !!
  ! Returns the detail of the matrix element in terms of single, two and three electron contribution.
  !! WARNING !!
  !
  ! Non hermitian !!
  !
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out) :: hmono, htwoe, hthree, htot

  integer                       :: degree

  PROVIDE pure_three_body_h_tc

  hmono  = 0.d0
  htwoe  = 0.d0
  htot   = 0.d0
  hthree = 0.d0

  call get_excitation_degree(key_i, key_j, degree, Nint)

  if(.not.pure_three_body_h_tc) then

    if(degree .gt. 2) return

    if(degree == 0) then
      call diag_htilde_mu_mat_fock_bi_ortho  (Nint, key_i, hmono, htwoe, hthree, htot)
    else if (degree == 1) then
      call single_htilde_mu_mat_fock_bi_ortho(Nint, key_j, key_i, hmono, htwoe, hthree, htot)
    else if(degree == 2) then
      call double_htilde_mu_mat_fock_bi_ortho(Nint, key_j, key_i, hmono, htwoe, hthree, htot)
    endif

  else

    if(degree .gt. 3) return

    if(degree == 0) then
      call diag_htilde_mu_mat_fock_bi_ortho  (Nint, key_i, hmono, htwoe, hthree, htot)
    else if (degree == 1) then
      call single_htilde_mu_mat_fock_bi_ortho(Nint, key_j, key_i, hmono, htwoe, hthree, htot)
    else if(degree == 2) then
      call double_htilde_mu_mat_fock_bi_ortho(Nint, key_j, key_i, hmono, htwoe, hthree, htot)
    else
      call triple_htilde_mu_mat_fock_bi_ortho(Nint, key_j, key_i, hmono, htwoe, hthree, htot)
    endif

  endif

  if(degree==0) then
    htot += nuclear_repulsion
  endif

end

! ---

subroutine htilde_mu_mat_opt_bi_ortho_no_3e(key_j, key_i, Nint, htot)

  BEGIN_DOC
  !
  ! <key_j |H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
  !!
  ! Returns the detail of the matrix element WITHOUT ANY CONTRIBUTION FROM THE THREE ELECTRON TERMS
  !! WARNING !!
  !
  ! Non hermitian !!
  !
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out) :: htot
  integer                       :: degree

  htot = 0.d0

  call get_excitation_degree(key_i, key_j, degree, Nint)
  if(degree.gt.2) return

  if(degree == 0) then
    call diag_htilde_mu_mat_fock_bi_ortho_no_3e(Nint, key_i,htot)
  else if (degree == 1) then
    call single_htilde_mu_mat_fock_bi_ortho_no_3e(Nint,key_j, key_i , htot)
  else if(degree == 2) then
    call double_htilde_mu_mat_fock_bi_ortho_no_3e(Nint, key_j, key_i, htot)
  endif

  if(degree==0) then
    htot += nuclear_repulsion
  endif

end

! ---

subroutine htilde_mu_mat_opt_bi_ortho_no_3e_both(key_j, key_i, Nint, hji,hij)

  BEGIN_DOC
  !
  ! <key_j |H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
  !!
  ! Returns the detail of the matrix element WITHOUT ANY CONTRIBUTION FROM THE THREE ELECTRON TERMS
  !! WARNING !!
  !
  ! Non hermitian !!
  !
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out) :: hji,hij
  integer                       :: degree

  hji = 0.d0
  hij = 0.d0

  call get_excitation_degree(key_i, key_j, degree, Nint)
  if(degree.gt.2) return

  if(degree == 0) then
    call diag_htilde_mu_mat_fock_bi_ortho_no_3e(Nint, key_i,hji)
    hij = hji
  else if (degree == 1) then
    call single_htilde_mu_mat_fock_bi_ortho_no_3e_both(Nint,key_j, key_i , hji,hij)
  else if(degree == 2) then
    call double_htilde_mu_mat_fock_bi_ortho_no_3e_both(Nint, key_j, key_i, hji,hij)
  endif

  if(degree==0) then
    hji += nuclear_repulsion
    hij += nuclear_repulsion
  endif

end

! ---

