
! ---

 BEGIN_PROVIDER [ double precision, ref_tc_energy_tot]
&BEGIN_PROVIDER [ double precision, ref_tc_energy_1e]
&BEGIN_PROVIDER [ double precision, ref_tc_energy_2e]
&BEGIN_PROVIDER [ double precision, ref_tc_energy_3e]

  BEGIN_DOC
  !
  ! Various component of the TC energy for the reference "HF" Slater determinant
  !
  END_DOC 

  implicit none
  double precision :: hmono, htwoe, htot, hthree

  PROVIDE N_int
  PROVIDE HF_bitmask
  PROVIDE mo_l_coef mo_r_coef

  call diag_htc_bi_orth_2e_brute(N_int, HF_bitmask, hmono, htwoe, htot)

  ref_tc_energy_1e = hmono
  ref_tc_energy_2e = htwoe 

  if(three_body_h_tc) then
    call diag_htc_bi_orth_3e_brute(N_int, HF_bitmask, hthree)
    ref_tc_energy_3e = hthree
  else
    ref_tc_energy_3e = 0.d0
  endif

  ref_tc_energy_tot = ref_tc_energy_1e + ref_tc_energy_2e + ref_tc_energy_3e + nuclear_repulsion

  if(noL_standard) then
    PROVIDE noL_0e
    ref_tc_energy_tot += noL_0e
  endif

END_PROVIDER 

! ---

subroutine diag_htilde_mu_mat_fock_bi_ortho(Nint, det_in, hmono, htwoe, hthree, htot)

  BEGIN_DOC
  !
  ! Computes $\langle i|H|i \rangle$.
  !
  END_DOC

  implicit none
  integer,           intent(in)  :: Nint
  integer(bit_kind), intent(in)  :: det_in(Nint,2)
  double precision,  intent(out) :: hmono, htwoe, htot, hthree

  integer(bit_kind)              :: hole(Nint,2)
  integer(bit_kind)              :: particle(Nint,2)
  integer                        :: i, nexc(2), ispin
  integer                        :: occ_particle(Nint*bit_kind_size,2)
  integer                        :: occ_hole(Nint*bit_kind_size,2)
  integer(bit_kind)              :: det_tmp(Nint,2)
  integer                        :: na, nb

  ASSERT (Nint > 0)
  ASSERT (sum(popcnt(det_in(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(det_in(:,2))) == elec_beta_num)

  nexc(1) = 0
  nexc(2) = 0
  do i = 1, Nint
    hole(i,1)     = xor(det_in(i,1),ref_bitmask(i,1))
    hole(i,2)     = xor(det_in(i,2),ref_bitmask(i,2))
    particle(i,1) = iand(hole(i,1),det_in(i,1))
    particle(i,2) = iand(hole(i,2),det_in(i,2))
    hole(i,1)     = iand(hole(i,1),ref_bitmask(i,1))
    hole(i,2)     = iand(hole(i,2),ref_bitmask(i,2))
    nexc(1)       = nexc(1) + popcnt(hole(i,1))
    nexc(2)       = nexc(2) + popcnt(hole(i,2))
  enddo

  if (nexc(1)+nexc(2) == 0) then
    hmono  = ref_tc_energy_1e
    htwoe  = ref_tc_energy_2e
    hthree = ref_tc_energy_3e
    htot   = ref_tc_energy_tot
    return
  endif

  !call debug_det(det_in,Nint)
  integer  :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(particle, occ_particle, tmp, Nint)
  ASSERT (tmp(1) == nexc(1)) ! Number of particles alpha
  ASSERT (tmp(2) == nexc(2)) ! Number of particle beta 
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(hole, occ_hole, tmp, Nint)
  ASSERT (tmp(1) == nexc(1)) ! Number of holes alpha
  ASSERT (tmp(2) == nexc(2)) ! Number of holes beta 

  hmono  = ref_tc_energy_1e
  htwoe  = ref_tc_energy_2e 
  hthree = ref_tc_energy_3e

  det_tmp = ref_bitmask

  do ispin = 1, 2
    na = elec_num_tab(ispin)
    nb = elec_num_tab(iand(ispin,1)+1)
    do i = 1, nexc(ispin)
      !DIR$ FORCEINLINE
      call ac_tc_operator(occ_particle(i,ispin), ispin, det_tmp, hmono, htwoe, hthree, Nint, na, nb)
      !DIR$ FORCEINLINE
      call a_tc_operator (occ_hole    (i,ispin), ispin, det_tmp, hmono, htwoe, hthree, Nint, na, nb)
    enddo
  enddo

  htot = hmono + htwoe + hthree + nuclear_repulsion

  if(noL_standard) then
    PROVIDE noL_0e
    htot += noL_0e
  endif

end

! ---

subroutine ac_tc_operator(iorb, ispin, key, hmono, htwoe, hthree, Nint, na, nb)

  BEGIN_DOC
  !
  ! Routine that computes one- and two-body energy corresponding 
  ! 
  ! to the ADDITION of an electron in an orbital 'iorb' of spin 'ispin' 
  ! 
  ! onto a determinant 'key'.
  !
  ! in output, the determinant key is changed by the ADDITION of that electron 
  !
  ! and the quantities hmono,htwoe,hthree are INCREMENTED 
  !
  END_DOC

  use bitmasks
  implicit none
  integer,           intent(in)    :: iorb, ispin, Nint
  integer,           intent(inout) :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision,  intent(inout) :: hmono, htwoe, hthree

  integer                          :: occ(Nint*bit_kind_size,2)
  integer                          :: other_spin
  integer                          :: k, l, i, jj, mm, j, m
  integer                          :: tmp(2)
  double precision                 :: direct_int, exchange_int
  
  if (iorb < 1) then
    print *,  irp_here, ': iorb < 1'
    print *,  iorb, mo_num
    stop -1
  endif
  if (iorb > mo_num) then
    print *,  irp_here, ': iorb > mo_num'
    print *,  iorb, mo_num
    stop -1
  endif

  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)

  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(key, occ, tmp, Nint)
  ASSERT (tmp(1) == elec_alpha_num)
  ASSERT (tmp(2) == elec_beta_num)

  k = shiftr(iorb-1,bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - shiftl(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  key(k,ispin) = ibset(key(k,ispin),l)
  other_spin = iand(ispin,1)+1

  hmono = hmono + mo_bi_ortho_tc_one_e(iorb,iorb)

  ! Same spin
  do i = 1, na
    htwoe = htwoe + mo_bi_ortho_tc_two_e_jj_anti(occ(i,ispin),iorb)
  enddo

  ! Opposite spin
  do i = 1, nb
    htwoe = htwoe + mo_bi_ortho_tc_two_e_jj(occ(i,other_spin),iorb)
  enddo

  if(three_body_h_tc .and. (elec_num.gt.2) .and. three_e_3_idx_term) then
    !!!!! 3-e part 

    !! same-spin/same-spin
    do j = 1, na
      jj = occ(j,ispin)
      do m = j+1, na
        mm = occ(m,ispin)
        hthree += three_e_diag_parrallel_spin_prov(mm,jj,iorb)
      enddo
    enddo
    !! same-spin/oposite-spin
    do j = 1, na
      jj = occ(j,ispin)
      do m = 1, nb
        mm = occ(m,other_spin)
        direct_int   = three_e_3_idx_direct_bi_ort(mm,jj,iorb) ! USES 3-IDX TENSOR 
        exchange_int = three_e_3_idx_exch12_bi_ort(mm,jj,iorb) ! USES 3-IDX TENSOR 
        hthree += direct_int - exchange_int
      enddo
    enddo
    !! oposite-spin/opposite-spin
    do j = 1, nb
      jj = occ(j,other_spin) 
      do m = j+1, nb 
        mm = occ(m,other_spin) 
        direct_int   = three_e_3_idx_direct_bi_ort(mm,jj,iorb) ! USES 3-IDX TENSOR 
        exchange_int = three_e_3_idx_exch23_bi_ort(mm,jj,iorb) ! USES 3-IDX TENSOR 
        hthree += direct_int - exchange_int
      enddo
    enddo
  endif

  na = na + 1

end

! ---

subroutine a_tc_operator(iorb, ispin, key, hmono, htwoe, hthree, Nint, na, nb)

  use bitmasks
  implicit none

  BEGIN_DOC
  !
  ! Routine that computes one- and two-body energy corresponding 
  ! 
  ! to the REMOVAL of an electron in an orbital 'iorb' of spin 'ispin' 
  ! 
  ! onto a determinant 'key'.
  !
  ! in output, the determinant key is changed by the REMOVAL of that electron 
  !
  ! and the quantities hmono,htwoe,hthree are INCREMENTED 
  !
  END_DOC

  integer,           intent(in)    :: iorb, ispin, Nint
  integer,           intent(inout) :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision,  intent(inout) :: hmono,htwoe,hthree
  
  double precision                 :: direct_int, exchange_int
  integer                          :: occ(Nint*bit_kind_size,2)
  integer                          :: other_spin
  integer                          :: k, l, i, jj, mm, j, m
  integer                          :: tmp(2)

  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)

  k = shiftr(iorb-1,bit_kind_shift)+1
  ASSERT (k>0)
  l = iorb - shiftl(k-1,bit_kind_shift)-1
  key(k,ispin) = ibclr(key(k,ispin),l)
  other_spin = iand(ispin,1)+1

  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(key, occ, tmp, Nint)
  na = na-1

  hmono = hmono - mo_bi_ortho_tc_one_e(iorb,iorb)

  ! Same spin
  do i = 1, na
    htwoe = htwoe - mo_bi_ortho_tc_two_e_jj_anti(occ(i,ispin),iorb)
  enddo

  ! Opposite spin
  do i = 1, nb
    htwoe = htwoe - mo_bi_ortho_tc_two_e_jj(occ(i,other_spin),iorb)
  enddo

  if(three_body_h_tc .and. elec_num.gt.2 .and. three_e_3_idx_term) then
    !!!!! 3-e part 

    !! same-spin/same-spin
    do j = 1, na
      jj = occ(j,ispin)
      do m = j+1, na
        mm = occ(m,ispin)
        hthree -= three_e_diag_parrallel_spin_prov(mm,jj,iorb)
      enddo
    enddo
    !! same-spin/oposite-spin
    do j = 1, na
      jj = occ(j,ispin)
      do m = 1, nb
        mm = occ(m,other_spin)
        direct_int   = three_e_3_idx_direct_bi_ort(mm,jj,iorb) ! USES 3-IDX TENSOR 
        exchange_int = three_e_3_idx_exch12_bi_ort(mm,jj,iorb) ! USES 3-IDX TENSOR 
        hthree -= (direct_int - exchange_int)
      enddo
    enddo 
    !! oposite-spin/opposite-spin
    do j = 1, nb
      jj = occ(j,other_spin) 
      do m = j+1, nb 
        mm = occ(m,other_spin) 
        direct_int   = three_e_3_idx_direct_bi_ort(mm,jj,iorb) ! USES 3-IDX TENSOR 
        exchange_int = three_e_3_idx_exch23_bi_ort(mm,jj,iorb) ! USES 3-IDX TENSOR 
        hthree -= (direct_int - exchange_int)
      enddo
    enddo
  endif

end

! ---

subroutine diag_htilde_mu_mat_fock_bi_ortho_no_3e(Nint, det_in,htot)

  BEGIN_DOC
  ! Computes $\langle i|H|i \rangle$. WITHOUT ANY CONTRIBUTIONS FROM 3E TERMS
  END_DOC

  implicit none
  integer,           intent(in)  :: Nint
  integer(bit_kind), intent(in)  :: det_in(Nint,2)
  double precision,  intent(out) :: htot
  double precision               :: hmono, htwoe
  integer(bit_kind)              :: hole(Nint,2)
  integer(bit_kind)              :: particle(Nint,2)
  integer                        :: i, nexc(2), ispin
  integer                        :: occ_particle(Nint*bit_kind_size,2)
  integer                        :: occ_hole(Nint*bit_kind_size,2)
  integer(bit_kind)              :: det_tmp(Nint,2)
  integer                        :: na, nb

  ASSERT (Nint > 0)
  ASSERT (sum(popcnt(det_in(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(det_in(:,2))) == elec_beta_num)


  nexc(1) = 0
  nexc(2) = 0
  do i=1,Nint
    hole(i,1)     = xor(det_in(i,1),ref_bitmask(i,1))
    hole(i,2)     = xor(det_in(i,2),ref_bitmask(i,2))
    particle(i,1) = iand(hole(i,1),det_in(i,1))
    particle(i,2) = iand(hole(i,2),det_in(i,2))
    hole(i,1)     = iand(hole(i,1),ref_bitmask(i,1))
    hole(i,2)     = iand(hole(i,2),ref_bitmask(i,2))
    nexc(1)       = nexc(1) + popcnt(hole(i,1))
    nexc(2)       = nexc(2) + popcnt(hole(i,2))
  enddo

  if(nexc(1)+nexc(2) == 0) then
    hmono = ref_tc_energy_1e
    htwoe = ref_tc_energy_2e
    htot  = ref_tc_energy_tot
    return
  endif

  !call debug_det(det_in,Nint)
  integer  :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(particle, occ_particle, tmp, Nint)
  ASSERT (tmp(1) == nexc(1)) ! Number of particles alpha
  ASSERT (tmp(2) == nexc(2)) ! Number of particle beta 
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(hole, occ_hole, tmp, Nint)
  ASSERT (tmp(1) == nexc(1)) ! Number of holes alpha
  ASSERT (tmp(2) == nexc(2)) ! Number of holes beta 

  det_tmp = ref_bitmask
  
  hmono = ref_tc_energy_1e
  htwoe = ref_tc_energy_2e 
  do ispin=1,2
    na = elec_num_tab(ispin)
    nb = elec_num_tab(iand(ispin,1)+1)
    do i=1,nexc(ispin)
      !DIR$ FORCEINLINE
      call ac_tc_operator_no_3e( occ_particle(i,ispin), ispin, det_tmp, hmono,htwoe, Nint,na,nb)
      !DIR$ FORCEINLINE
      call a_tc_operator_no_3e ( occ_hole    (i,ispin), ispin, det_tmp, hmono,htwoe, Nint,na,nb)
    enddo
  enddo
  htot = hmono+htwoe
end

subroutine ac_tc_operator_no_3e(iorb,ispin,key,hmono,htwoe,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Routine that computes one- and two-body energy corresponding 
  ! 
  ! to the ADDITION of an electron in an orbital 'iorb' of spin 'ispin' 
  ! 
  ! onto a determinant 'key'.
  !
  ! in output, the determinant key is changed by the ADDITION of that electron 
  !
  ! and the quantities hmono,htwoe are INCREMENTED 
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hmono,htwoe

  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i,jj,mm,j,m
  double precision ::  direct_int, exchange_int
  

  if (iorb < 1) then
    print *,  irp_here, ': iorb < 1'
    print *,  iorb, mo_num
    stop -1
  endif
  if (iorb > mo_num) then
    print *,  irp_here, ': iorb > mo_num'
    print *,  iorb, mo_num
    stop -1
  endif

  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)

  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(key, occ, tmp, Nint)
  ASSERT (tmp(1) == elec_alpha_num)
  ASSERT (tmp(2) == elec_beta_num)

  k = shiftr(iorb-1,bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - shiftl(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  key(k,ispin) = ibset(key(k,ispin),l)
  other_spin = iand(ispin,1)+1

  hmono = hmono + mo_bi_ortho_tc_one_e(iorb,iorb)

  ! Same spin
  do i=1,na
    htwoe = htwoe + mo_bi_ortho_tc_two_e_jj_anti(occ(i,ispin),iorb)
  enddo

  ! Opposite spin
  do i=1,nb
    htwoe = htwoe + mo_bi_ortho_tc_two_e_jj(occ(i,other_spin),iorb)
  enddo

  na = na+1
end

subroutine a_tc_operator_no_3e(iorb,ispin,key,hmono,htwoe,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Routine that computes one- and two-body energy corresponding 
  ! 
  ! to the REMOVAL of an electron in an orbital 'iorb' of spin 'ispin' 
  ! 
  ! onto a determinant 'key'.
  !
  ! in output, the determinant key is changed by the REMOVAL of that electron 
  !
  ! and the quantities hmono,htwoe are INCREMENTED 
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hmono,htwoe
  
  double precision  :: direct_int, exchange_int
  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i,jj,mm,j,m
  integer                        :: tmp(2)

  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)

  k = shiftr(iorb-1,bit_kind_shift)+1
  ASSERT (k>0)
  l = iorb - shiftl(k-1,bit_kind_shift)-1
  key(k,ispin) = ibclr(key(k,ispin),l)
  other_spin = iand(ispin,1)+1

  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(key, occ, tmp, Nint)
  na = na-1

  hmono = hmono - mo_bi_ortho_tc_one_e(iorb,iorb)

  ! Same spin
  do i = 1, na
    htwoe = htwoe- mo_bi_ortho_tc_two_e_jj_anti(occ(i,ispin),iorb)
  enddo

  ! Opposite spin
  do i = 1, nb
    htwoe = htwoe- mo_bi_ortho_tc_two_e_jj(occ(i,other_spin),iorb)
  enddo

end

! ---

subroutine diag_htc_bi_orth_2e_brute(Nint, key_i, hmono, htwoe, htot)

  BEGIN_DOC
  !
  ! diagonal element of htilde ONLY FOR ONE- AND TWO-BODY TERMS 
  !
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

  PROVIDE mo_bi_ortho_tc_two_e

  hmono = 0.d0
  htwoe = 0.d0
  htot  = 0.d0

  call bitstring_to_list_ab(key_i, occ, Ne, Nint)

  do ispin = 1, 2
    do i = 1, Ne(ispin)
      ii = occ(i,ispin)
      hmono += mo_bi_ortho_tc_one_e(ii,ii)
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

! ---                                                                                           

subroutine diag_htc_bi_orth_3e_brute(Nint, key_i, hthree)

  BEGIN_DOC
  !  diagonal element of htilde ONLY FOR THREE-BODY TERMS WITH BI ORTHONORMAL ORBITALS
  END_DOC

  use bitmasks

  implicit none
  integer,           intent(in) :: Nint
  integer(bit_kind), intent(in) :: key_i(Nint,2)
  double precision, intent(out) :: hthree
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2),i,j,ii,jj,ispin,jspin,m,mm
  integer(bit_kind)             :: key_i_core(Nint,2)
  double precision              :: direct_int, exchange_int, ref
  double precision, external    :: sym_3_e_int_from_6_idx_tensor
  double precision, external    :: three_e_diag_parrallel_spin

  PROVIDE mo_l_coef mo_r_coef

  if(core_tc_op) then
    do i = 1, Nint
      key_i_core(i,1) = xor(key_i(i,1), core_bitmask(i,1))
      key_i_core(i,2) = xor(key_i(i,2), core_bitmask(i,2))
    enddo
    call bitstring_to_list_ab(key_i_core, occ, Ne, Nint)
  else
    call bitstring_to_list_ab(key_i, occ, Ne, Nint)
  endif

  hthree = 0.d0

  if((Ne(1)+Ne(2)) .ge. 3) then

    ! alpha/alpha/beta three-body
    do i = 1, Ne(1)
      ii = occ(i,1)
      do j = i+1, Ne(1)
        jj = occ(j,1)
        do m = 1, Ne(2)
          mm = occ(m,2)
          !direct_int   = three_body_ints_bi_ort(mm,jj,ii,mm,jj,ii) !uses the 6-idx tensor
          !exchange_int = three_body_ints_bi_ort(mm,jj,ii,mm,ii,jj) !uses the 6-idx tensor
          direct_int   = three_e_3_idx_direct_bi_ort(mm,jj,ii)      !uses 3-idx tensor
          exchange_int = three_e_3_idx_exch12_bi_ort(mm,jj,ii)      !uses 3-idx tensor
          hthree      += direct_int - exchange_int
        enddo
      enddo
    enddo

    ! beta/beta/alpha three-body
    do i = 1, Ne(2)
      ii = occ(i,2)
      do j = i+1, Ne(2)
        jj = occ(j,2)
        do m = 1, Ne(1)
          mm = occ(m,1)
          !direct_int   = three_body_ints_bi_ort(mm,jj,ii,mm,jj,ii) !uses the 6-idx tensor
          !exchange_int = three_body_ints_bi_ort(mm,jj,ii,mm,ii,jj) !uses the 6-idx tensor
          direct_int   = three_e_3_idx_direct_bi_ort(mm,jj,ii)
          exchange_int = three_e_3_idx_exch12_bi_ort(mm,jj,ii)
          hthree      += direct_int - exchange_int
        enddo
      enddo
    enddo

    ! alpha/alpha/alpha three-body
    do i = 1, Ne(1)
      ii = occ(i,1) ! 1
      do j = i+1, Ne(1)
        jj = occ(j,1) ! 2
        do m = j+1, Ne(1)
          mm = occ(m,1) ! 3
          !hthree += sym_3_e_int_from_6_idx_tensor(mm,jj,ii,mm,jj,ii) !uses the 6 idx tensor
          hthree += three_e_diag_parrallel_spin(mm,jj,ii)             !uses only 3-idx tensors
        enddo
      enddo
    enddo

    ! beta/beta/beta three-body
    do i = 1, Ne(2)
      ii = occ(i,2) ! 1
      do j = i+1, Ne(2)
        jj = occ(j,2) ! 2
        do m = j+1, Ne(2)
          mm = occ(m,2) ! 3
          !hthree += sym_3_e_int_from_6_idx_tensor(mm,jj,ii,mm,jj,ii) !uses the 6 idx tensor
          hthree += three_e_diag_parrallel_spin(mm,jj,ii)             !uses only 3-idx tensors
        enddo
      enddo
    enddo

  endif

end



BEGIN_PROVIDER [ double precision, three_e_diag_parrallel_spin_prov, (mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator ON A BI ORTHONORMAL BASIS 
  !
  ! three_e_diag_parrallel_spin_prov(m,j,i) = All combinations of the form <mji|-L|mji> for same spin matrix elements  
  ! 
  ! notice the -1 sign: in this way three_e_diag_parrallel_spin_prov can be directly used to compute Slater rules with a + sign
  !
  END_DOC

  implicit none
  integer          :: i, j, m
  double precision :: integral, wall1, wall0, three_e_diag_parrallel_spin

  three_e_diag_parrallel_spin_prov = 0.d0
  print *, ' Providing the three_e_diag_parrallel_spin_prov ...'

 integral = three_e_diag_parrallel_spin(1,1,1) ! to provide all stuffs
  call wall_time(wall0)
 !$OMP PARALLEL                 &
 !$OMP DEFAULT (NONE)           &
 !$OMP PRIVATE (i,j,m,integral) & 
 !$OMP SHARED (mo_num,three_e_diag_parrallel_spin_prov)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do j = 1, mo_num
      do m = j, mo_num
        three_e_diag_parrallel_spin_prov(m,j,i) =  three_e_diag_parrallel_spin(m,j,i)
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  do i = 1, mo_num
    do j = 1, mo_num
      do m = 1, j
        three_e_diag_parrallel_spin_prov(m,j,i) = three_e_diag_parrallel_spin_prov(j,m,i)
      enddo
    enddo
  enddo

  call wall_time(wall1)
  print *, ' wall time for three_e_diag_parrallel_spin_prov', wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, three_e_single_parrallel_spin_prov, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF SINGLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_single_parrallel_spin_prov(m,j,k,i) = All combination of <mjk|-L|mji> for same spin matrix elements 
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  !
  END_DOC

 implicit none
 integer          :: i, j, k, m
 double precision :: integral, wall1, wall0, three_e_single_parrallel_spin

  three_e_single_parrallel_spin_prov = 0.d0
  print *, ' Providing the three_e_single_parrallel_spin_prov ...'

  integral = three_e_single_parrallel_spin(1,1,1,1)
  call wall_time(wall0)
 !$OMP PARALLEL                   &
 !$OMP DEFAULT (NONE)             &
 !$OMP PRIVATE (i,j,k,m,integral) & 
 !$OMP SHARED (mo_num,three_e_single_parrallel_spin_prov)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do m = 1, mo_num
          three_e_single_parrallel_spin_prov(m,j,k,i) = three_e_single_parrallel_spin(m,j,k,i)
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_single_parrallel_spin_prov', wall1 - wall0

END_PROVIDER 


! ---

BEGIN_PROVIDER [ double precision, three_e_double_parrallel_spin_prov, (mo_num, mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! matrix element of the -L  three-body operator FOR THE DIRECT TERMS OF DOUBLE EXCITATIONS AND BI ORTHO MOs
  !
  ! three_e_double_parrallel_spin_prov(m,l,j,k,i) = <mlk|-L|mji> ::: notice that i is the RIGHT MO and k is the LEFT MO
  !
  ! notice the -1 sign: in this way three_e_3_idx_direct_bi_ort can be directly used to compute Slater rules with a + sign
  END_DOC

  implicit none
  integer          :: i, j, k, m, l
  double precision :: integral, wall1, wall0, three_e_double_parrallel_spin

  three_e_double_parrallel_spin_prov = 0.d0
  print *, ' Providing the three_e_double_parrallel_spin_prov ...'
  call wall_time(wall0)

 integral = three_e_double_parrallel_spin(1,1,1,1,1)
 !$OMP PARALLEL                     &
 !$OMP DEFAULT (NONE)               &
 !$OMP PRIVATE (i,j,k,m,l,integral) & 
 !$OMP SHARED (mo_num,three_e_double_parrallel_spin_prov)
 !$OMP DO SCHEDULE (dynamic)
  do i = 1, mo_num
    do k = 1, mo_num
      do j = 1, mo_num
        do l = 1, mo_num
          do m = 1, mo_num
            three_e_double_parrallel_spin_prov(m,l,j,k,i) = three_e_double_parrallel_spin(m,l,j,k,i)
          enddo
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print *, ' wall time for three_e_double_parrallel_spin_prov', wall1 - wall0

END_PROVIDER 

