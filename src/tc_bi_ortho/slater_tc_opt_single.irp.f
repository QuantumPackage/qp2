

subroutine  single_htilde_mu_mat_fock_bi_ortho (Nint, key_j, key_i, hmono, htwoe, hthree, htot)
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
  double precision, intent(out) :: hmono, htwoe, hthree, htot
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2), i, j, ii, jj, ispin, jspin, k, kk
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  double precision              :: get_mo_two_e_integral_tc_int, phase
  double precision              :: direct_int, exchange_int_12, exchange_int_23, exchange_int_13
  integer                       :: other_spin(2)
  integer(bit_kind)             :: key_j_core(Nint,2), key_i_core(Nint,2)

  other_spin(1) = 2
  other_spin(2) = 1

  hmono  = 0.d0
  htwoe  = 0.d0
  hthree = 0.d0
  htot   = 0.d0
  call get_excitation_degree(key_i, key_j, degree, Nint)
  if(degree.ne.1)then
   return
  endif
  call bitstring_to_list_ab(key_i, occ, Ne, Nint)

  call get_single_excitation(key_i, key_j, exc, phase, Nint)
  call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)
  call get_single_excitation_from_fock_tc(key_i,key_j,h1,p1,s1,phase,hmono,htwoe,hthree,htot)
end


subroutine get_single_excitation_from_fock_tc(key_i,key_j,h,p,spin,phase,hmono,htwoe,hthree,htot)
 use bitmasks
 implicit none
 integer,intent(in) :: h,p,spin
 double precision, intent(in)  :: phase
 integer(bit_kind), intent(in) :: key_i(N_int,2), key_j(N_int,2)
 double precision, intent(out) :: hmono,htwoe,hthree,htot
 integer(bit_kind) :: differences(N_int,2)
 integer(bit_kind) :: hole(N_int,2)
 integer(bit_kind) :: partcl(N_int,2)
 integer :: occ_hole(N_int*bit_kind_size,2)
 integer :: occ_partcl(N_int*bit_kind_size,2)
 integer :: n_occ_ab_hole(2),n_occ_ab_partcl(2)
 integer :: i0,i
 double precision :: buffer_c(mo_num),buffer_x(mo_num)
 do i=1, mo_num
   buffer_c(i) = tc_2e_3idx_coulomb_integrals(i,p,h)
   buffer_x(i) = tc_2e_3idx_exchange_integrals(i,p,h)
 enddo
 do i = 1, N_int
  differences(i,1) = xor(key_i(i,1),ref_closed_shell_bitmask(i,1))
  differences(i,2) = xor(key_i(i,2),ref_closed_shell_bitmask(i,2))
  hole(i,1) = iand(differences(i,1),ref_closed_shell_bitmask(i,1))
  hole(i,2) = iand(differences(i,2),ref_closed_shell_bitmask(i,2))
  partcl(i,1) = iand(differences(i,1),key_i(i,1))
  partcl(i,2) = iand(differences(i,2),key_i(i,2))
 enddo
 call bitstring_to_list_ab(hole, occ_hole, n_occ_ab_hole, N_int)
 call bitstring_to_list_ab(partcl, occ_partcl, n_occ_ab_partcl, N_int)
 hmono = mo_bi_ortho_tc_one_e(p,h)
 htwoe = fock_op_2_e_tc_closed_shell(p,h)
 ! holes :: direct terms
 do i0 = 1, n_occ_ab_hole(1)
  i = occ_hole(i0,1)
  htwoe -= buffer_c(i)
 enddo
 do i0 = 1, n_occ_ab_hole(2)
  i = occ_hole(i0,2)
  htwoe -= buffer_c(i)
 enddo

 ! holes :: exchange terms
 do i0 = 1, n_occ_ab_hole(spin)
  i = occ_hole(i0,spin)
  htwoe += buffer_x(i)
 enddo

 ! particles :: direct terms
 do i0 = 1, n_occ_ab_partcl(1)
  i = occ_partcl(i0,1)
  htwoe += buffer_c(i)
 enddo
 do i0 = 1, n_occ_ab_partcl(2)
  i = occ_partcl(i0,2)
  htwoe += buffer_c(i)
 enddo

 ! particles :: exchange terms
 do i0 = 1, n_occ_ab_partcl(spin)
  i = occ_partcl(i0,spin)
  htwoe -= buffer_x(i)
 enddo
 hthree = 0.d0
 if (three_body_h_tc)then
  call three_comp_fock_elem(key_i,h,p,spin,hthree)
 endif


 htwoe = htwoe * phase
 hmono = hmono * phase
 hthree = hthree * phase
 htot  = htwoe + hmono + hthree

end

subroutine three_comp_fock_elem(key_i,h_fock,p_fock,ispin_fock,hthree)
 implicit none
 integer,intent(in) :: h_fock,p_fock,ispin_fock
 integer(bit_kind), intent(in) :: key_i(N_int,2)
 double precision, intent(out) :: hthree
 integer :: nexc(2),i,ispin,na,nb
 integer(bit_kind) :: hole(N_int,2)
 integer(bit_kind) :: particle(N_int,2)
 integer :: occ_hole(N_int*bit_kind_size,2)
 integer :: occ_particle(N_int*bit_kind_size,2)
 integer :: n_occ_ab_hole(2),n_occ_ab_particle(2)
 integer(bit_kind)              :: det_tmp(N_int,2)


  nexc(1) = 0
  nexc(2) = 0
  !! Get all the holes and particles of key_i with respect to the ROHF determinant
  do i=1,N_int
    hole(i,1)     = xor(key_i(i,1),ref_bitmask(i,1))
    hole(i,2)     = xor(key_i(i,2),ref_bitmask(i,2))
    particle(i,1) = iand(hole(i,1),key_i(i,1))
    particle(i,2) = iand(hole(i,2),key_i(i,2))
    hole(i,1)     = iand(hole(i,1),ref_bitmask(i,1))
    hole(i,2)     = iand(hole(i,2),ref_bitmask(i,2))
    nexc(1)       = nexc(1) + popcnt(hole(i,1))
    nexc(2)       = nexc(2) + popcnt(hole(i,2))
  enddo
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(particle, occ_particle, tmp, N_int)
  ASSERT (tmp(1) == nexc(1)) ! Number of particles alpha
  ASSERT (tmp(2) == nexc(2)) ! Number of particle beta 
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(hole, occ_hole, tmp, N_int)
  ASSERT (tmp(1) == nexc(1)) ! Number of holes alpha
  ASSERT (tmp(2) == nexc(2)) ! Number of holes beta 

  !! Initialize the matrix element with the reference ROHF Slater determinant Fock element
  if(ispin_fock==1)then
   hthree = fock_a_tot_3e_bi_orth(p_fock,h_fock) 
  else 
   hthree = fock_b_tot_3e_bi_orth(p_fock,h_fock) 
  endif
  det_tmp = ref_bitmask
  do ispin=1,2
    na = elec_num_tab(ispin)
    nb = elec_num_tab(iand(ispin,1)+1)
    do i=1,nexc(ispin)
      !DIR$ FORCEINLINE
      call fock_ac_tc_operator( occ_particle(i,ispin), ispin, det_tmp, h_fock,p_fock, ispin_fock, hthree, N_int,na,nb)
      !DIR$ FORCEINLINE
      call fock_a_tc_operator ( occ_hole    (i,ispin), ispin, det_tmp, h_fock,p_fock, ispin_fock, hthree, N_int,na,nb)
    enddo
  enddo
end

subroutine fock_ac_tc_operator(iorb,ispin,key, h_fock,p_fock, ispin_fock,hthree,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Routine that computes the contribution to the three-electron part of the Fock operator 
  !
  ! a^dagger_{p_fock} a_{h_fock} of spin ispin_fock
  ! 
  ! on top of a determinant 'key' on which you ADD an electron of spin ispin in orbital iorb
  ! 
  ! in output, the determinant key is changed by the ADDITION of that electron 
  !
  ! the output hthree is INCREMENTED
  END_DOC
  integer, intent(in)              :: iorb, ispin, Nint, h_fock,p_fock, ispin_fock
  integer, intent(inout)           :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout)  :: hthree

  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i,jj,j
  double precision :: direct_int, exchange_int
  

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


  !! spin of other electrons == ispin 
  if(ispin == ispin_fock)then
   !! in what follows :: jj == other electrons in the determinant 
   !!                 :: iorb == electron that has been added of spin ispin
   !!                 :: p_fock, h_fock == hole particle of spin ispin_fock
   !! jj = ispin = ispin_fock >> pure parallel spin
   do j = 1, na
    jj = occ(j,ispin)
    hthree += three_e_single_parrallel_spin_prov(jj,iorb,p_fock,h_fock)
   enddo
   !! spin of jj == other spin than ispin AND ispin_fock
   !! exchange between the iorb and (h_fock, p_fock)
   do j = 1, nb
    jj = occ(j,other_spin) 
    direct_int = three_e_4_idx_direct_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    exchange_int = three_e_4_idx_exch12_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    hthree += direct_int - exchange_int
   enddo
  else !! ispin NE to ispin_fock
   !! jj = ispin BUT NON EQUAL TO ispin_fock 
   !! exchange between the jj and iorb
   do j = 1, na
    jj = occ(j,ispin)
    direct_int   = three_e_4_idx_direct_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    exchange_int = three_e_4_idx_exch23_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    hthree += direct_int - exchange_int
   enddo
   !! jj = other_spin than ispin BUT jj == ispin_fock
   !! exchange between jj and (h_fock,p_fock)
   do j = 1, nb
    jj = occ(j,other_spin) 
    direct_int = three_e_4_idx_direct_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    exchange_int = three_e_4_idx_exch13_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    hthree += direct_int - exchange_int
   enddo
  endif

  na = na+1
end

subroutine fock_a_tc_operator(iorb,ispin,key, h_fock,p_fock, ispin_fock,hthree,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Routine that computes the contribution to the three-electron part of the Fock operator 
  !
  ! a^dagger_{p_fock} a_{h_fock} of spin ispin_fock
  ! 
  ! on top of a determinant 'key' on which you REMOVE an electron of spin ispin in orbital iorb
  ! 
  ! in output, the determinant key is changed by the REMOVAL of that electron 
  !
  ! the output hthree is INCREMENTED
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint, h_fock,p_fock, ispin_fock
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hthree
  
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
  !! spin of other electrons == ispin 
  if(ispin == ispin_fock)then
   !! in what follows :: jj == other electrons in the determinant 
   !!                 :: iorb == electron that has been added of spin ispin
   !!                 :: p_fock, h_fock == hole particle of spin ispin_fock
   !! jj = ispin = ispin_fock >> pure parallel spin
   do j = 1, na
    jj = occ(j,ispin)
    hthree -= three_e_single_parrallel_spin_prov(jj,iorb,p_fock,h_fock)
   enddo
   !! spin of jj == other spin than ispin AND ispin_fock
   !! exchange between the iorb and (h_fock, p_fock)
   do j = 1, nb
    jj = occ(j,other_spin) 
    direct_int = three_e_4_idx_direct_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    exchange_int = three_e_4_idx_exch12_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    hthree -= direct_int - exchange_int
   enddo
  else !! ispin NE to ispin_fock
   !! jj = ispin BUT NON EQUAL TO ispin_fock 
   !! exchange between the jj and iorb
   do j = 1, na
    jj = occ(j,ispin)
    direct_int   = three_e_4_idx_direct_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    exchange_int = three_e_4_idx_exch23_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    hthree -= direct_int - exchange_int
   enddo
   !! jj = other_spin than ispin BUT jj == ispin_fock
   !! exchange between jj and (h_fock,p_fock)
   do j = 1, nb
    jj = occ(j,other_spin) 
    direct_int = three_e_4_idx_direct_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    exchange_int = three_e_4_idx_exch13_bi_ort(jj,iorb,p_fock,h_fock) ! USES 4-IDX TENSOR 
    hthree -= direct_int - exchange_int
   enddo
  endif

end


BEGIN_PROVIDER [double precision, fock_op_2_e_tc_closed_shell, (mo_num, mo_num) ]
 implicit none
 BEGIN_DOC
! Closed-shell part of the Fock operator for the TC operator
 END_DOC
 integer :: h0,p0,h,p,k0,k,i
 integer :: n_occ_ab(2)
 integer :: occ(N_int*bit_kind_size,2)
 integer :: n_occ_ab_virt(2)
 integer :: occ_virt(N_int*bit_kind_size,2)
 integer(bit_kind) :: key_test(N_int)
 integer(bit_kind) :: key_virt(N_int,2)
 double precision :: accu

 fock_op_2_e_tc_closed_shell = -1000.d0
 call bitstring_to_list_ab(ref_closed_shell_bitmask, occ, n_occ_ab, N_int)
 do i = 1, N_int
  key_virt(i,1) = full_ijkl_bitmask(i)
  key_virt(i,2) = full_ijkl_bitmask(i)
  key_virt(i,1) = xor(key_virt(i,1),ref_closed_shell_bitmask(i,1))
  key_virt(i,2) = xor(key_virt(i,2),ref_closed_shell_bitmask(i,2))
 enddo
 call bitstring_to_list_ab(key_virt, occ_virt, n_occ_ab_virt, N_int)
 ! docc ---> virt single excitations
 do h0 = 1,  n_occ_ab(1)
  h=occ(h0,1)
  do p0 = 1, n_occ_ab_virt(1)
   p = occ_virt(p0,1)
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * tc_2e_3idx_coulomb_integrals(k,p,h) - tc_2e_3idx_exchange_integrals(k,p,h)
   enddo
   fock_op_2_e_tc_closed_shell(p,h) = accu 
  enddo
 enddo

 do h0 = 1, n_occ_ab_virt(1)
  h = occ_virt(h0,1)
  do p0 = 1,  n_occ_ab(1)
   p=occ(p0,1)
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * tc_2e_3idx_coulomb_integrals(k,p,h) - tc_2e_3idx_exchange_integrals(k,p,h)
   enddo
   fock_op_2_e_tc_closed_shell(p,h) = accu 
  enddo
 enddo

 ! virt ---> virt single excitations
 do h0 = 1,  n_occ_ab_virt(1)
  h=occ_virt(h0,1)
  do p0 = 1, n_occ_ab_virt(1)
   p = occ_virt(p0,1)
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * tc_2e_3idx_coulomb_integrals(k,p,h) - tc_2e_3idx_exchange_integrals(k,p,h)
   enddo
   fock_op_2_e_tc_closed_shell(p,h) = accu 
  enddo
 enddo

 do h0 = 1, n_occ_ab_virt(1)
  h = occ_virt(h0,1)
  do p0 = 1,  n_occ_ab_virt(1)
   p=occ_virt(p0,1)
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * tc_2e_3idx_coulomb_integrals(k,p,h) - tc_2e_3idx_exchange_integrals(k,p,h)
   enddo
   fock_op_2_e_tc_closed_shell(p,h) = accu 
  enddo
 enddo


 ! docc ---> docc single excitations
 do h0 = 1,  n_occ_ab(1)
  h=occ(h0,1)
  do p0 = 1, n_occ_ab(1)
   p = occ(p0,1)
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * tc_2e_3idx_coulomb_integrals(k,p,h) - tc_2e_3idx_exchange_integrals(k,p,h)
   enddo
   fock_op_2_e_tc_closed_shell(p,h) = accu 
  enddo
 enddo

 do h0 = 1, n_occ_ab(1)
  h = occ(h0,1)
  do p0 = 1,  n_occ_ab(1)
   p=occ(p0,1)
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * tc_2e_3idx_coulomb_integrals(k,p,h) - tc_2e_3idx_exchange_integrals(k,p,h)
   enddo
   fock_op_2_e_tc_closed_shell(p,h) = accu 
  enddo
 enddo

! do i = 1, mo_num
!  write(*,'(100(F10.5,X))')fock_op_2_e_tc_closed_shell(:,i)
! enddo

END_PROVIDER


subroutine  single_htilde_mu_mat_fock_bi_ortho_no_3e(Nint, key_j, key_i, htot)
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
  double precision, intent(out) ::  htot
  double precision :: hmono, htwoe
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2), i, j, ii, jj, ispin, jspin, k, kk
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  double precision              :: get_mo_two_e_integral_tc_int, phase
  double precision              :: direct_int, exchange_int_12, exchange_int_23, exchange_int_13
  integer                       :: other_spin(2)
  integer(bit_kind)             :: key_j_core(Nint,2), key_i_core(Nint,2)

  other_spin(1) = 2
  other_spin(2) = 1

  hmono  = 0.d0
  htwoe  = 0.d0
  htot   = 0.d0
  call get_excitation_degree(key_i, key_j, degree, Nint)
  if(degree.ne.1)then
   return
  endif
  call bitstring_to_list_ab(key_i, occ, Ne, Nint)

  call get_single_excitation(key_i, key_j, exc, phase, Nint)
  call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)
  call get_single_excitation_from_fock_tc_no_3e(key_i,key_j,h1,p1,s1,phase,hmono,htwoe,htot)
end


subroutine get_single_excitation_from_fock_tc_no_3e(key_i,key_j,h,p,spin,phase,hmono,htwoe,htot)
 use bitmasks
 implicit none
 integer,intent(in) :: h,p,spin
 double precision, intent(in)  :: phase
 integer(bit_kind), intent(in) :: key_i(N_int,2), key_j(N_int,2)
 double precision, intent(out) :: hmono,htwoe,htot
 integer(bit_kind) :: differences(N_int,2)
 integer(bit_kind) :: hole(N_int,2)
 integer(bit_kind) :: partcl(N_int,2)
 integer :: occ_hole(N_int*bit_kind_size,2)
 integer :: occ_partcl(N_int*bit_kind_size,2)
 integer :: n_occ_ab_hole(2),n_occ_ab_partcl(2)
 integer :: i0,i
 double precision :: buffer_c(mo_num),buffer_x(mo_num)
 do i=1, mo_num
   buffer_c(i) = tc_2e_3idx_coulomb_integrals(i,p,h)
   buffer_x(i) = tc_2e_3idx_exchange_integrals(i,p,h)
 enddo
 do i = 1, N_int
  differences(i,1) = xor(key_i(i,1),ref_closed_shell_bitmask(i,1))
  differences(i,2) = xor(key_i(i,2),ref_closed_shell_bitmask(i,2))
  hole(i,1) = iand(differences(i,1),ref_closed_shell_bitmask(i,1))
  hole(i,2) = iand(differences(i,2),ref_closed_shell_bitmask(i,2))
  partcl(i,1) = iand(differences(i,1),key_i(i,1))
  partcl(i,2) = iand(differences(i,2),key_i(i,2))
 enddo
 call bitstring_to_list_ab(hole, occ_hole, n_occ_ab_hole, N_int)
 call bitstring_to_list_ab(partcl, occ_partcl, n_occ_ab_partcl, N_int)
 hmono = mo_bi_ortho_tc_one_e(p,h)
 htwoe = fock_op_2_e_tc_closed_shell(p,h)
 ! holes :: direct terms
 do i0 = 1, n_occ_ab_hole(1)
  i = occ_hole(i0,1)
  htwoe -= buffer_c(i)
 enddo
 do i0 = 1, n_occ_ab_hole(2)
  i = occ_hole(i0,2)
  htwoe -= buffer_c(i)
 enddo

 ! holes :: exchange terms
 do i0 = 1, n_occ_ab_hole(spin)
  i = occ_hole(i0,spin)
  htwoe += buffer_x(i)
 enddo

 ! particles :: direct terms
 do i0 = 1, n_occ_ab_partcl(1)
  i = occ_partcl(i0,1)
  htwoe += buffer_c(i)
 enddo
 do i0 = 1, n_occ_ab_partcl(2)
  i = occ_partcl(i0,2)
  htwoe += buffer_c(i)
 enddo

 ! particles :: exchange terms
 do i0 = 1, n_occ_ab_partcl(spin)
  i = occ_partcl(i0,spin)
  htwoe -= buffer_x(i)
 enddo
 htwoe = htwoe * phase
 hmono = hmono * phase
 htot  = htwoe + hmono 

end

