
! --

subroutine hmat_bi_ortho(key_j, key_i, Nint, hmono, htwoe, htot)

  BEGIN_DOC
  !
  ! < key_j | H | key_i > where | key_j > is developed on the LEFT basis and | key_i > is developed on the RIGHT basis
  ! 
  END_DOC

  use bitmasks

  implicit none

  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out) :: hmono, htwoe, htot

  integer                       :: degree 

  hmono = 0.d0
  htwoe = 0.d0
  htot  = 0.d0

  call get_excitation_degree(key_i, key_j, degree, Nint)
  if(degree .gt. 2) return

  if(degree == 0) then

    call diag_hmat_bi_ortho(Nint, key_i, hmono, htwoe)
    htot = htot + nuclear_repulsion

  else if (degree == 1) then

    call single_hmat_bi_ortho(Nint, key_j, key_i, hmono, htwoe)

  else if(degree == 2) then

    call double_hmat_bi_ortho(Nint, key_j, key_i, hmono, htwoe)

  endif

  htot += hmono + htwoe

  return
end subroutine hmat_bi_ortho

! ---

subroutine diag_hmat_bi_ortho(Nint, key_i, hmono, htwoe)

  use bitmasks

  implicit none

  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_i(Nint,2)
  double precision, intent(out) :: hmono, htwoe

  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: Ne(2), i, j, ii, jj, ispin, jspin

  hmono = 0.d0
  htwoe = 0.d0

  call bitstring_to_list_ab(key_i, occ, Ne, Nint)

  do ispin = 1, 2 
    do i = 1, Ne(ispin)
      ii = occ(i,ispin) 
      hmono += mo_bi_ortho_one_e(ii,ii)
    enddo
  enddo

  ! alpha/beta two-body
  ispin = 1
  jspin = 2 
  do i = 1, Ne(ispin) ! electron 1
    ii = occ(i,ispin) 
    do j = 1, Ne(jspin) ! electron 2 
      jj = occ(j,jspin) 
      htwoe += mo_bi_ortho_coul_e(jj,ii,jj,ii) 
    enddo
  enddo
 
  ! alpha/alpha two-body
  do i = 1, Ne(ispin)
    ii = occ(i,ispin) 
    do j = i+1, Ne(ispin)
      jj = occ(j,ispin) 
      htwoe += mo_bi_ortho_coul_e(ii,jj,ii,jj) - mo_bi_ortho_coul_e(ii,jj,jj,ii)
    enddo
  enddo
 
  ! beta/beta two-body
  do i = 1, Ne(jspin)
    ii = occ(i,jspin) 
    do j = i+1, Ne(jspin)
      jj = occ(j,jspin) 
      htwoe += mo_bi_ortho_coul_e(ii,jj,ii,jj) - mo_bi_ortho_coul_e(ii,jj,jj,ii)
    enddo
  enddo

  return
end subroutine diag_hmat_bi_ortho

! ---

subroutine single_hmat_bi_ortho(Nint, key_j, key_i, hmono, htwoe)

  BEGIN_DOC
  !
  ! < key_j | H | key_i > for single excitation 
  ! 
  END_DOC

  use bitmasks

  implicit none

  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_j(Nint,2), key_i(Nint,2)
  double precision, intent(out) :: hmono, htwoe

  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2), i, j, ii, ispin, jspin
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  integer                       :: other_spin(2)
  double precision              :: phase

  other_spin(1) = 2
  other_spin(2) = 1

  hmono = 0.d0
  htwoe = 0.d0

  call get_excitation_degree(key_i, key_j, degree, Nint)
  if(degree .ne. 1) then
    return
  endif

  call bitstring_to_list_ab(key_i, occ, Ne, Nint)

  call get_single_excitation(key_i, key_j, exc, phase, Nint)
  call decode_exc(exc, 1, h1, p1, h2, p2, s1, s2)

  hmono = mo_bi_ortho_one_e(p1,h1) * phase

  ! alpha/beta two-body 
  ispin = other_spin(s1)
  if(s1 == 1) then

    ! single alpha 
    do i = 1, Ne(ispin) ! electron 2 
      ii = occ(i,ispin) 
      htwoe += mo_bi_ortho_coul_e(ii,p1,ii,h1) 
    enddo

  else

    ! single beta 
    do i = 1, Ne(ispin) ! electron 1 
      ii = occ(i,ispin) 
      htwoe += mo_bi_ortho_coul_e(p1,ii,h1,ii) 
    enddo

  endif

  ! same spin two-body 
  do i = 1, Ne(s1)
    ii = occ(i,s1) 
    ! ( h1 p1 |ii ii ) - ( h1 ii | p1 ii )
    htwoe += mo_bi_ortho_coul_e(ii,p1,ii,h1) - mo_bi_ortho_coul_e(p1,ii,ii,h1) 
  enddo
   
  htwoe *= phase

end subroutine single_hmat_bi_ortho

! ---

subroutine double_hmat_bi_ortho(Nint, key_j, key_i, hmono, htwoe)

  BEGIN_DOC
  !
  ! < key_j | H | key_i> for double excitation
  ! 
  END_DOC

  use bitmasks

  implicit none

  integer,           intent(in) :: Nint 
  integer(bit_kind), intent(in) :: key_j(Nint,2), key_i(Nint,2)
  double precision, intent(out) :: hmono, htwoe

  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2), i, j, ii, ispin, jspin
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  integer                       :: other_spin(2)
  double precision              :: phase

  other_spin(1) = 2
  other_spin(2) = 1

  call get_excitation_degree(key_i, key_j, degree, Nint)

  hmono = 0.d0
  htwoe = 0.d0

  if(degree .ne. 2) then
    return
  endif

  call bitstring_to_list_ab(key_i, occ, Ne, Nint)

  call get_double_excitation(key_i, key_j, exc, phase, Nint)
  call decode_exc(exc, 2, h1, p1, h2, p2, s1, s2)

  if(s1 .ne. s2) then

    htwoe = mo_bi_ortho_coul_e(p2,p1,h2,h1) 

  else

    ! same spin two-body 

    !                    direct terms                 exchange terms 
    htwoe = mo_bi_ortho_coul_e(p2,p1,h2,h1) - mo_bi_ortho_coul_e(p1,p2,h2,h1) 

  endif

  htwoe *= phase

end subroutine double_hmat_bi_ortho

! ---


