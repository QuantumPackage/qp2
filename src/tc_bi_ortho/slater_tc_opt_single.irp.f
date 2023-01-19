

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
 htwoe = htwoe * phase
 hmono = hmono * phase
 htot  = htwoe + hmono + hthree

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

