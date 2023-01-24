 BEGIN_PROVIDER [ double precision, ref_tc_energy_tot]
&BEGIN_PROVIDER [ double precision, ref_tc_energy_1e]
&BEGIN_PROVIDER [ double precision, ref_tc_energy_2e]
&BEGIN_PROVIDER [ double precision, ref_tc_energy_3e]
 implicit none
 BEGIN_DOC
! Various component of the TC energy for the reference "HF" Slater determinant
 END_DOC 
 double precision :: hmono, htwoe, htot, hthree
 call diag_htilde_mu_mat_bi_ortho(N_int,HF_bitmask , hmono, htwoe, htot)
 ref_tc_energy_1e = hmono
 ref_tc_energy_2e = htwoe 
 if(three_body_h_tc)then
  call diag_htilde_three_body_ints_bi_ort(N_int, HF_bitmask, hthree)
  ref_tc_energy_3e = hthree
 else
  ref_tc_energy_3e = 0.d0
 endif
 ref_tc_energy_tot = ref_tc_energy_1e + ref_tc_energy_2e + ref_tc_energy_3e
 END_PROVIDER 

subroutine diag_htilde_mu_mat_fock_bi_ortho(Nint, det_in, hmono, htwoe, hthree, htot)
  implicit none
  BEGIN_DOC
  ! Computes $\langle i|H|i \rangle$.
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)
  double precision, intent(out)  :: hmono,htwoe,htot,hthree

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

  if (nexc(1)+nexc(2) == 0) then
    hmono = ref_tc_energy_1e
    htwoe = ref_tc_energy_2e
    hthree= ref_tc_energy_3e
    htot = ref_tc_energy_tot
    return
  endif

  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
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
  hthree= ref_tc_energy_3e
  do ispin=1,2
    na = elec_num_tab(ispin)
    nb = elec_num_tab(iand(ispin,1)+1)
    do i=1,nexc(ispin)
      !DIR$ FORCEINLINE
      call ac_tc_operator( occ_particle(i,ispin), ispin, det_tmp, hmono,htwoe,hthree, Nint,na,nb)
      !DIR$ FORCEINLINE
      call a_tc_operator ( occ_hole    (i,ispin), ispin, det_tmp, hmono,htwoe,hthree, Nint,na,nb)
    enddo
  enddo
  htot = hmono+htwoe+hthree
end

subroutine ac_tc_operator(iorb,ispin,key,hmono,htwoe,hthree,Nint,na,nb)
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
  ! and the quantities hmono,htwoe,hthree are INCREMENTED 
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hmono,htwoe,hthree

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

  if(three_body_h_tc)then
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
      direct_int = three_e_3_idx_direct_bi_ort(mm,jj,iorb) ! USES 3-IDX TENSOR 
      exchange_int = three_e_3_idx_exch23_bi_ort(mm,jj,iorb) ! USES 3-IDX TENSOR 
      hthree += direct_int - exchange_int
     enddo
    enddo
  endif

  na = na+1
end

subroutine a_tc_operator(iorb,ispin,key,hmono,htwoe,hthree,Nint,na,nb)
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
  ! and the quantities hmono,htwoe,hthree are INCREMENTED 
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hmono,htwoe,hthree
  
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
  do i=1,na
    htwoe= htwoe- mo_bi_ortho_tc_two_e_jj_anti(occ(i,ispin),iorb)
  enddo

  ! Opposite spin
  do i=1,nb
    htwoe= htwoe- mo_bi_ortho_tc_two_e_jj(occ(i,other_spin),iorb)
  enddo

  if(three_body_h_tc)then
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


subroutine diag_htilde_mu_mat_fock_bi_ortho_no_3e(Nint, det_in,htot)
  implicit none
  BEGIN_DOC
  ! Computes $\langle i|H|i \rangle$. WITHOUT ANY CONTRIBUTIONS FROM 3E TERMS
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)
  double precision, intent(out)  :: htot
  double precision :: hmono,htwoe

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

  if (nexc(1)+nexc(2) == 0) then
    hmono = ref_tc_energy_1e
    htwoe = ref_tc_energy_2e
    htot = ref_tc_energy_tot
    return
  endif

  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
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
  do i=1,na
    htwoe= htwoe- mo_bi_ortho_tc_two_e_jj_anti(occ(i,ispin),iorb)
  enddo

  ! Opposite spin
  do i=1,nb
    htwoe= htwoe- mo_bi_ortho_tc_two_e_jj(occ(i,other_spin),iorb)
  enddo

end

