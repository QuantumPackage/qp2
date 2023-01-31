
! ---

subroutine htilde_mu_mat_bi_ortho_tot(key_j, key_i, Nint, htot)

  BEGIN_DOC
  ! <key_j | H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer, intent(in)           :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2),key_i(Nint,2)
  double precision, intent(out) :: htot
  double precision              :: hmono, htwoe, hthree
  integer :: degree

  call get_excitation_degree(key_j, key_i, degree, Nint)
  if(degree.gt.2)then
    htot = 0.d0
  else
    call htilde_mu_mat_bi_ortho(key_j, key_i, Nint, hmono, htwoe, hthree, htot)
  endif

end subroutine htilde_mu_mat_bi_ortho_tot

! --

subroutine htilde_mu_mat_bi_ortho(key_j, key_i, Nint, hmono, htwoe, hthree, htot)

  BEGIN_DOC
  !
  ! <key_j | H_tilde | key_i> where |key_j> is developed on the LEFT basis and |key_i> is developed on the RIGHT basis
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

  hmono  = 0.d0
  htwoe  = 0.d0
  htot   = 0.d0
  hthree = 0.D0

  call get_excitation_degree(key_i, key_j, degree, Nint)
  if(degree.gt.2) return

  if(degree == 0)then
    call diag_htilde_mu_mat_bi_ortho(Nint, key_i, hmono, htwoe, htot)
  else if (degree == 1)then
    call single_htilde_mu_mat_bi_ortho(Nint, key_j, key_i, hmono, htwoe, htot)
  else if(degree == 2)then
    call double_htilde_mu_mat_bi_ortho(Nint, key_j, key_i, hmono, htwoe, htot)
  endif

  if(three_body_h_tc) then
    if(degree == 2) then
      if(.not.double_normal_ord) then
        call double_htilde_three_body_ints_bi_ort(Nint, key_j, key_i, hthree)
      endif
    else if(degree == 1) then
      call single_htilde_three_body_ints_bi_ort(Nint, key_j, key_i, hthree)
    else if(degree == 0) then
      call diag_htilde_three_body_ints_bi_ort(Nint, key_i, hthree)
    endif
  endif

  htot = hmono + htwoe + hthree
  if(degree==0) then
    htot += nuclear_repulsion
  endif
 
end

! ---

subroutine diag_htilde_mu_mat_bi_ortho(Nint, key_i, hmono, htwoe, htot)

  BEGIN_DOC
  !  diagonal element of htilde ONLY FOR ONE- AND TWO-BODY TERMS 
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in)  :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2)
  double precision, intent(out)  :: hmono,htwoe,htot
  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: Ne(2), i, j, ii, jj, ispin, jspin, k, kk
  double precision               :: get_mo_two_e_integral_tc_int
  integer(bit_kind)              :: key_i_core(Nint,2)

!  PROVIDE mo_two_e_integrals_tc_int_in_map mo_bi_ortho_tc_two_e
!
!  PROVIDE mo_integrals_erf_map core_energy nuclear_repulsion core_bitmask
!  PROVIDE core_fock_operator
!
!  PROVIDE j1b_gauss

!  if(core_tc_op)then
!   print*,'core_tc_op not already taken into account for bi ortho'
!   print*,'stopping ...'
!   stop
!   do i = 1, Nint
!    key_i_core(i,1) = xor(key_i(i,1),core_bitmask(i,1))
!    key_i_core(i,2) = xor(key_i(i,2),core_bitmask(i,2))
!   enddo
!   call bitstring_to_list_ab(key_i_core, occ, Ne, Nint)
!   hmono = core_energy - nuclear_repulsion
!  else
   call bitstring_to_list_ab(key_i, occ, Ne, Nint)
   hmono = 0.d0
!  endif
  htwoe= 0.d0
  htot = 0.d0

  do ispin = 1, 2 
   do i = 1, Ne(ispin) ! 
    ii = occ(i,ispin) 
    hmono += mo_bi_ortho_tc_one_e(ii,ii)

!    if(j1b_gauss .eq. 1) then
!      print*,'j1b not implemented for bi ortho TC'
!      print*,'stopping  ....'
!      stop
!      !hmono += mo_j1b_gauss_hermI  (ii,ii) &
!      !       + mo_j1b_gauss_hermII (ii,ii) &
!      !       + mo_j1b_gauss_nonherm(ii,ii)
!    endif

!    if(core_tc_op)then
!   print*,'core_tc_op not already taken into account for bi ortho'
!   print*,'stopping ...'
!   stop
!     hmono += core_fock_operator(ii,ii) ! add the usual Coulomb - Exchange from the core 
!    endif
   enddo
  enddo


   ! alpha/beta two-body
   ispin = 1
   jspin = 2 
   do i = 1, Ne(ispin) ! electron 1 (so it can be associated to mu(r1))
    ii = occ(i,ispin) 
    do j = 1, Ne(jspin) ! electron 2 
     jj = occ(j,jspin) 
     htwoe += mo_bi_ortho_tc_two_e(jj,ii,jj,ii) 
    enddo
   enddo
 
   ! alpha/alpha two-body
   do i = 1, Ne(ispin)
    ii = occ(i,ispin) 
    do j = i+1, Ne(ispin)
     jj = occ(j,ispin) 
     htwoe += mo_bi_ortho_tc_two_e(ii,jj,ii,jj) - mo_bi_ortho_tc_two_e(ii,jj,jj,ii)
    enddo
   enddo
 
   ! beta/beta two-body
   do i = 1, Ne(jspin)
    ii = occ(i,jspin) 
    do j = i+1, Ne(jspin)
     jj = occ(j,jspin) 
     htwoe += mo_bi_ortho_tc_two_e(ii,jj,ii,jj) - mo_bi_ortho_tc_two_e(ii,jj,jj,ii)
    enddo
   enddo
  htot = hmono + htwoe 

end



subroutine double_htilde_mu_mat_bi_ortho(Nint, key_j, key_i, hmono, htwoe, htot)

  BEGIN_DOC
  ! <key_j | H_tilde | key_i> for double excitation  ONLY FOR ONE- AND TWO-BODY TERMS 
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint 
  integer(bit_kind), intent(in) :: key_j(Nint,2), key_i(Nint,2)
  double precision, intent(out) :: hmono, htwoe, htot
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2), i, j, ii, jj, ispin, jspin, k, kk
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  integer                       :: other_spin(2)
  integer(bit_kind)             :: key_i_core(Nint,2)
  double precision              :: get_mo_two_e_integral_tc_int,phase

!  PROVIDE mo_two_e_integrals_tc_int_in_map mo_bi_ortho_tc_two_e

  other_spin(1) = 2
  other_spin(2) = 1

  call get_excitation_degree(key_i, key_j, degree, Nint)

  hmono = 0.d0
  htwoe= 0.d0
  htot = 0.d0

  if(degree.ne.2)then
   return
  endif

!  if(core_tc_op)then
!   print*,'core_tc_op not already taken into account for bi ortho'
!   print*,'stopping ...'
!   stop
!   do i = 1, Nint
!    key_i_core(i,1) = xor(key_i(i,1),core_bitmask(i,1))
!    key_i_core(i,2) = xor(key_i(i,2),core_bitmask(i,2))
!   enddo
!   call bitstring_to_list_ab(key_i_core, occ, Ne, Nint)
!  else
   call bitstring_to_list_ab(key_i, occ, Ne, Nint)
!  endif
  call get_double_excitation(key_i, key_j, exc, phase, Nint)
  call decode_exc(exc, 2, h1, p1, h2, p2, s1, s2)

  if(s1.ne.s2)then
   ! opposite spin two-body 
!   key_j, key_i
    htwoe  = mo_bi_ortho_tc_two_e(p2,p1,h2,h1) 
    if(double_normal_ord.and.+Ne(1).gt.2)then
     htwoe += normal_two_body_bi_orth(p2,h2,p1,h1)!!! WTF ???
    endif
  else
   ! same spin two-body 
   ! direct terms 
   htwoe  = mo_bi_ortho_tc_two_e(p2,p1,h2,h1)  
   ! exchange terms 
   htwoe -= mo_bi_ortho_tc_two_e(p1,p2,h2,h1) 
   if(double_normal_ord.and.+Ne(1).gt.2)then
    htwoe -= normal_two_body_bi_orth(h2,p1,h1,p2)!!! WTF ???
    htwoe += normal_two_body_bi_orth(h1,p1,h2,p2)!!! WTF ???
   endif
  endif
  htwoe *= phase
  htot =  htwoe 

end


subroutine single_htilde_mu_mat_bi_ortho(Nint, key_j, key_i, hmono, htwoe, htot)

  BEGIN_DOC
  ! <key_j | H_tilde | key_i> for single excitation ONLY FOR ONE- AND TWO-BODY TERMS 
  !!
  !! WARNING !!
  ! 
  ! Non hermitian !!
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2), key_i(Nint,2)
  double precision, intent(out) :: hmono, htwoe, htot
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2), i, j, ii, jj, ispin, jspin, k, kk
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  double precision              :: get_mo_two_e_integral_tc_int, phase
  double precision              :: direct_int, exchange_int_12, exchange_int_23, exchange_int_13
  integer                       :: other_spin(2)
  integer(bit_kind)             :: key_j_core(Nint,2), key_i_core(Nint,2)

!  PROVIDE mo_two_e_integrals_tc_int_in_map mo_bi_ortho_tc_two_e
!
!  PROVIDE core_bitmask core_fock_operator mo_integrals_erf_map

!  PROVIDE j1b_gauss

  other_spin(1) = 2
  other_spin(2) = 1

  hmono = 0.d0
  htwoe= 0.d0
  htot = 0.d0
  call get_excitation_degree(key_i, key_j, degree, Nint)
  if(degree.ne.1)then
   return
  endif
!  if(core_tc_op)then
!   print*,'core_tc_op not already taken into account for bi ortho'
!   print*,'stopping ...'
!   stop
!   do i = 1, Nint
!    key_i_core(i,1) = xor(key_i(i,1),core_bitmask(i,1))
!    key_i_core(i,2) = xor(key_i(i,2),core_bitmask(i,2))
!    key_j_core(i,1) = xor(key_j(i,1),core_bitmask(i,1))
!    key_j_core(i,2) = xor(key_j(i,2),core_bitmask(i,2))
!   enddo
!   call bitstring_to_list_ab(key_i_core, occ, Ne, Nint)
!  else
   call bitstring_to_list_ab(key_i, occ, Ne, Nint)
!  endif

  call get_single_excitation(key_i, key_j, exc, phase, Nint)
  call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)
!  if(h1==14.and.p1==2)then
!   print*,'h1,p1 old = ',h1,p1
!  endif

  hmono = mo_bi_ortho_tc_one_e(p1,h1) * phase

!  if(j1b_gauss .eq. 1) then
!     print*,'j1b not implemented for bi ortho TC'
!     print*,'stopping  ....'
!     stop
!    !hmono += ( mo_j1b_gauss_hermI  (h1,p1) &
!    !         + mo_j1b_gauss_hermII (h1,p1) &
!    !         + mo_j1b_gauss_nonherm(h1,p1) ) * phase
!  endif

!  if(core_tc_op)then
!   print*,'core_tc_op not already taken into account for bi ortho'
!   print*,'stopping ...'
!   stop
!   hmono += phase * core_fock_operator(h1,p1)
!  endif
  
   ! alpha/beta two-body 
   ispin = other_spin(s1)
   if(s1==1)then
    ! single alpha 
    do i = 1, Ne(ispin) ! electron 2 
     ii = occ(i,ispin) 
     htwoe += mo_bi_ortho_tc_two_e(ii,p1,ii,h1) 
    enddo
   else
    ! single beta 
    do i = 1, Ne(ispin) ! electron 1 
     ii = occ(i,ispin) 
     htwoe += mo_bi_ortho_tc_two_e(p1,ii,h1,ii) 
    enddo
   endif
!   ! same spin two-body 
   do i = 1, Ne(s1)
    ii = occ(i,s1) 
    ! (h1p1|ii ii) - (h1 ii| p1 ii)
    htwoe += mo_bi_ortho_tc_two_e(ii,p1,ii,h1) - mo_bi_ortho_tc_two_e(p1,ii,ii,h1) 
   enddo
   
  htwoe  *= phase
  htot = hmono + htwoe 

end


