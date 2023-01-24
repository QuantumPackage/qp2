
subroutine double_htilde_mu_mat_fock_bi_ortho(Nint, key_j, key_i, hmono, htwoe, hthree, htot)

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
  double precision, intent(out) :: hmono, htwoe, hthree, htot
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2), i, j, ii, jj, ispin, jspin, k, kk
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  double precision              :: get_mo_two_e_integral_tc_int,phase


  call get_excitation_degree(key_i, key_j, degree, Nint)

  hmono  = 0.d0
  htwoe  = 0.d0
  hthree = 0.d0
  htot   = 0.d0

  if(degree.ne.2)then
   return
  endif
  integer :: degree_i,degree_j
  call get_excitation_degree(ref_bitmask,key_i,degree_i,N_int)
  call get_excitation_degree(ref_bitmask,key_j,degree_j,N_int)
  call get_double_excitation(key_i, key_j, exc, phase, Nint)
  call decode_exc(exc, 2, h1, p1, h2, p2, s1, s2)

  if(s1.ne.s2)then
   ! opposite spin two-body 
    htwoe  = mo_bi_ortho_tc_two_e(p2,p1,h2,h1) 
    if(three_body_h_tc)then
     if(.not.double_normal_ord)then
      if(degree_i>degree_j)then
       call three_comp_two_e_elem(key_j,h1,h2,p1,p2,s1,s2,hthree)
      else
       call three_comp_two_e_elem(key_i,h1,h2,p1,p2,s1,s2,hthree)
      endif
     elseif(double_normal_ord.and.elec_num+elec_num.gt.2)then
      htwoe += normal_two_body_bi_orth(p2,h2,p1,h1)!!! WTF ???
     endif
    endif
  else
   ! same spin two-body 
   ! direct terms 
   htwoe  = mo_bi_ortho_tc_two_e(p2,p1,h2,h1)  
   ! exchange terms 
   htwoe -= mo_bi_ortho_tc_two_e(p1,p2,h2,h1) 
   if(three_body_h_tc)then
    if(.not.double_normal_ord)then
     if(degree_i>degree_j)then
      call three_comp_two_e_elem(key_j,h1,h2,p1,p2,s1,s2,hthree)
     else
      call three_comp_two_e_elem(key_i,h1,h2,p1,p2,s1,s2,hthree)
     endif
    elseif(double_normal_ord.and.elec_num+elec_num.gt.2)then
      htwoe -= normal_two_body_bi_orth(h2,p1,h1,p2)!!! WTF ???
      htwoe += normal_two_body_bi_orth(h1,p1,h2,p2)!!! WTF ???
    endif
   endif
  endif
  hthree *= phase
  htwoe  *= phase
  htot    =  htwoe + hthree

end



subroutine three_comp_two_e_elem(key_i,h1,h2,p1,p2,s1,s2,hthree)
 implicit none
 integer(bit_kind), intent(in) :: key_i(N_int,2)
 integer, intent(in) :: h1,h2,p1,p2,s1,s2
 double precision, intent(out) :: hthree
 integer :: nexc(2),i,ispin,na,nb
 integer(bit_kind) :: hole(N_int,2)
 integer(bit_kind) :: particle(N_int,2)
 integer :: occ_hole(N_int*bit_kind_size,2)
 integer :: occ_particle(N_int*bit_kind_size,2)
 integer :: n_occ_ab_hole(2),n_occ_ab_particle(2)
 integer(bit_kind)              :: det_tmp(N_int,2)
 integer :: ipart, ihole
 double precision :: direct_int, exchange_int

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
  if(s1==s2.and.s1==1)then
   !!!!!!!!!!!!!!!!!!!!!!!!!! alpha/alpha double exc
   hthree = eff_2_e_from_3_e_aa(p2,p1,h2,h1) 
   if(nexc(1)+nexc(2) ==0)return !! if you're on the reference determinant 
    !!!!!!!! the matrix element is already exact 
    !!!!!!!! else you need to take care of holes and particles 
    !!!!!!!!!!!!! Holes and particles !!!!!!!!!!!!!!!!!!!!!!!
    ispin = 1 ! i==alpha ==> pure same spin terms
    do i = 1, nexc(ispin) ! number of couple of holes/particles 
     ipart=occ_particle(i,ispin)
     hthree += three_e_double_parrallel_spin_prov(ipart,p2,h2,p1,h1)
     ihole=occ_hole(i,ispin)
     hthree -= three_e_double_parrallel_spin_prov(ihole,p2,h2,p1,h1)
    enddo
    ispin = 2 ! i==beta ==> alpha/alpha/beta terms
    do i = 1, nexc(ispin) ! number of couple of holes/particles 
     ! exchange between (h1,p1) and (h2,p2)
     ipart=occ_particle(i,ispin)
     direct_int  = three_e_5_idx_direct_bi_ort(ipart,p2,h2,p1,h1)
     exchange_int = three_e_5_idx_exch12_bi_ort(ipart,p2,h2,p1,h1)
     hthree += direct_int - exchange_int
     ihole=occ_hole(i,ispin)
     direct_int  = three_e_5_idx_direct_bi_ort(ihole,p2,h2,p1,h1)
     exchange_int = three_e_5_idx_exch12_bi_ort(ihole,p2,h2,p1,h1)
     hthree -= direct_int - exchange_int
    enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif(s1==s2.and.s1==2)then 
   !!!!!!!!!!!!!!!!!!!!!!!!!! beta/beta double exc
   hthree = eff_2_e_from_3_e_bb(p2,p1,h2,h1)
   if(nexc(1)+nexc(2) ==0)return !! if you're on the reference determinant 
   !!!!!!!! the matrix element is already exact 
   !!!!!!!! else you need to take care of holes and particles 
   !!!!!!!!!!!!! Holes and particles !!!!!!!!!!!!!!!!!!!!!!!
   ispin = 2 ! i==beta  ==> pure same spin terms
   do i = 1, nexc(ispin) ! number of couple of holes/particles 
    ipart=occ_particle(i,ispin)
    hthree += three_e_double_parrallel_spin_prov(ipart,p2,h2,p1,h1)
    ihole=occ_hole(i,ispin)
    hthree -= three_e_double_parrallel_spin_prov(ihole,p2,h2,p1,h1)
   enddo
   ispin = 1 ! i==alpha==> beta/beta/alpha terms
   do i = 1, nexc(ispin) ! number of couple of holes/particles 
    ! exchange between (h1,p1) and (h2,p2)
    ipart=occ_particle(i,ispin)
    direct_int  = three_e_5_idx_direct_bi_ort(ipart,p2,h2,p1,h1)
    exchange_int = three_e_5_idx_exch12_bi_ort(ipart,p2,h2,p1,h1)
    hthree += direct_int - exchange_int
    ihole=occ_hole(i,ispin)
    direct_int  = three_e_5_idx_direct_bi_ort(ihole,p2,h2,p1,h1)
    exchange_int = three_e_5_idx_exch12_bi_ort(ihole,p2,h2,p1,h1)
    hthree -= direct_int - exchange_int
   enddo
  else                         ! (h1,p1) == alpha/(h2,p2) == beta 
   hthree = eff_2_e_from_3_e_ab(p2,p1,h2,h1)
   if(nexc(1)+nexc(2) ==0)return !! if you're on the reference determinant 
   !!!!!!!! the matrix element is already exact 
   !!!!!!!! else you need to take care of holes and particles 
   !!!!!!!!!!!!! Holes and particles !!!!!!!!!!!!!!!!!!!!!!!
   ispin = 1 ! i==alpha ==> alpha/beta/alpha terms 
   do i = 1, nexc(ispin) ! number of couple of holes/particles 
    ! exchange between (h1,p1) and i
    ipart=occ_particle(i,ispin)
    direct_int  = three_e_5_idx_direct_bi_ort(ipart,p2,h2,p1,h1)
    exchange_int = three_e_5_idx_exch13_bi_ort(ipart,p2,h2,p1,h1)
    hthree += direct_int - exchange_int
    ihole=occ_hole(i,ispin)
    direct_int  = three_e_5_idx_direct_bi_ort(ihole,p2,h2,p1,h1)
    exchange_int = three_e_5_idx_exch13_bi_ort(ihole,p2,h2,p1,h1)
    hthree -= direct_int - exchange_int
   enddo
   ispin = 2 ! i==beta  ==> alpha/beta/beta  terms 
   do i = 1, nexc(ispin) ! number of couple of holes/particles 
    ! exchange between (h2,p2) and i
    ipart=occ_particle(i,ispin)
    direct_int  = three_e_5_idx_direct_bi_ort(ipart,p2,h2,p1,h1)
    exchange_int = three_e_5_idx_exch23_bi_ort(ipart,p2,h2,p1,h1)
    hthree += direct_int - exchange_int
    ihole=occ_hole(i,ispin)
    direct_int  = three_e_5_idx_direct_bi_ort(ihole,p2,h2,p1,h1)
    exchange_int = three_e_5_idx_exch23_bi_ort(ihole,p2,h2,p1,h1)
    hthree -= direct_int - exchange_int
   enddo
  endif
end


BEGIN_PROVIDER [ double precision, eff_2_e_from_3_e_ab, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! eff_2_e_from_3_e_ab(p2,p1,h2,h1) = Effective Two-electron operator for alpha/beta double excitations 
!
! from contraction with HF density = a^{dagger}_p1_alpha a^{dagger}_p2_beta a_h2_beta a_h1_alpha
 END_DOC
 integer :: i,h1,p1,h2,p2
 integer :: hh1,hh2,pp1,pp2,m,mm
 integer                        :: Ne(2)
 integer,           allocatable :: occ(:,:)
 double precision :: contrib
 allocate( occ(N_int*bit_kind_size,2) )
 call bitstring_to_list_ab(ref_bitmask,occ,Ne,N_int)
 call give_contrib_for_abab(1,1,1,1,occ,Ne,contrib)
 eff_2_e_from_3_e_ab = 0.d0
 !$OMP PARALLEL                                                                         &
 !$OMP DEFAULT (NONE)                                                                   &
 !$OMP PRIVATE (hh1, h1, hh2, h2, pp1, p1, pp2, p2, contrib) & 
 !$OMP SHARED (n_act_orb, list_act, Ne,occ, eff_2_e_from_3_e_ab)
 !$OMP DO SCHEDULE (static) 
  do hh1 = 1, n_act_orb !! alpha 
    h1 = list_act(hh1) 
    do hh2 = 1, n_act_orb !! beta 
      h2 = list_act(hh2) 
      do pp1 = 1, n_act_orb !! alpha
        p1 = list_act(pp1)
        do pp2 = 1, n_act_orb !! beta 
          p2 = list_act(pp2)
          call give_contrib_for_abab(h1,h2,p1,p2,occ,Ne,contrib)
          eff_2_e_from_3_e_ab(p2,p1,h2,h1) = contrib
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

END_PROVIDER 

subroutine give_contrib_for_abab(h1,h2,p1,p2,occ,Ne,contrib)
 implicit none
 BEGIN_DOC 
! gives the contribution for a double excitation (h1,p1)_alpha (h2,p2)_beta
!
! on top of a determinant whose occupied orbitals is in (occ, Ne)
 END_DOC
 integer, intent(in) :: h1,h2,p1,p2,occ(N_int*bit_kind_size,2),Ne(2)
 double precision, intent(out) :: contrib
 integer :: mm,m 
 double precision :: direct_int, exchange_int
 !! h1,p1 == alpha 
 !! h2,p2 == beta
 contrib = 0.d0
 do mm = 1, Ne(1) !! alpha 
   m = occ(mm,1)
   direct_int   = three_e_5_idx_direct_bi_ort(mm,p2,h2,p1,h1) 
   ! exchange between (h1,p1) and m
   exchange_int = three_e_5_idx_exch13_bi_ort(mm,p2,h2,p1,h1)
   contrib += direct_int - exchange_int
 enddo

 do mm = 1, Ne(2) !! beta
   m = occ(mm,2)
   direct_int   = three_e_5_idx_direct_bi_ort(mm,p2,h2,p1,h1) 
   ! exchange between (h2,p2) and m
   exchange_int = three_e_5_idx_exch23_bi_ort(mm,p2,h2,p1,h1)
   contrib += direct_int - exchange_int
 enddo
end

BEGIN_PROVIDER [ double precision, eff_2_e_from_3_e_aa, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! eff_2_e_from_3_e_ab(p2,p1,h2,h1) = Effective Two-electron operator for alpha/alpha double excitations 
!
! from contractionelec_alpha_num with HF density = a^{dagger}_p1_alpha a^{dagger}_p2_alpha a_h2_alpha a_h1_alpha
!
! WARNING :: to be coherent with the phase convention used in the Hamiltonian matrix elements, you must fulfill 
!
! ||||    h2>h1, p2>p1   ||||
 END_DOC
 integer :: i,h1,p1,h2,p2
 integer :: hh1,hh2,pp1,pp2,m,mm
 integer                        :: Ne(2)
 integer,           allocatable :: occ(:,:)
 double precision :: contrib
 allocate( occ(N_int*bit_kind_size,2) )
 call bitstring_to_list_ab(ref_bitmask,occ,Ne,N_int)
 call give_contrib_for_aaaa(1 ,1 ,1 ,1 ,occ,Ne,contrib)
 eff_2_e_from_3_e_aa = 100000000.d0
 !$OMP PARALLEL                                                                         &
 !$OMP DEFAULT (NONE)                                                                   &
 !$OMP PRIVATE (hh1, h1, hh2, h2, pp1, p1, pp2, p2, contrib) & 
 !$OMP SHARED (n_act_orb, list_act, Ne,occ, eff_2_e_from_3_e_aa)
 !$OMP DO SCHEDULE (static) 
  do hh1 = 1, n_act_orb !! alpha 
    h1 = list_act(hh1) 
    do hh2 = hh1+1, n_act_orb !! alpha
      h2 = list_act(hh2) 
      do pp1 = 1, n_act_orb !! alpha
        p1 = list_act(pp1)
        do pp2 = pp1+1, n_act_orb !! alpha
          p2 = list_act(pp2)
          call give_contrib_for_aaaa(h1,h2,p1,p2,occ,Ne,contrib)
          eff_2_e_from_3_e_aa(p2,p1,h2,h1) = contrib
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

END_PROVIDER 

subroutine give_contrib_for_aaaa(h1,h2,p1,p2,occ,Ne,contrib)
 implicit none
 BEGIN_DOC 
! gives the contribution for a double excitation (h1,p1)_alpha (h2,p2)_alpha
!
! on top of a determinant whose occupied orbitals is in (occ, Ne)
 END_DOC
 integer, intent(in) :: h1,h2,p1,p2,occ(N_int*bit_kind_size,2),Ne(2)
 double precision, intent(out) :: contrib
 integer :: mm,m 
 double precision :: direct_int, exchange_int
 !! h1,p1 == alpha 
 !! h2,p2 == alpha
 contrib = 0.d0
 do mm = 1, Ne(1) !! alpha ==> pure parallele spin contribution
   m = occ(mm,1)
   contrib += three_e_double_parrallel_spin_prov(m,p2,h2,p1,h1)
 enddo

 do mm = 1, Ne(2) !! beta
   m = occ(mm,2)
   direct_int   = three_e_5_idx_direct_bi_ort(mm,p2,h2,p1,h1) 
   ! exchange between (h1,p1) and (h2,p2)
   exchange_int = three_e_5_idx_exch12_bi_ort(mm,p2,h2,p1,h1)
   contrib += direct_int - exchange_int
 enddo
end


BEGIN_PROVIDER [ double precision, eff_2_e_from_3_e_bb, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! eff_2_e_from_3_e_ab(p2,p1,h2,h1) = Effective Two-electron operator for beta/beta double excitations 
!
! from contractionelec_beta_num with HF density = a^{dagger}_p1_beta a^{dagger}_p2_beta a_h2_beta a_h1_beta
!
! WARNING :: to be coherent with the phase convention used in the Hamiltonian matrix elements, you must fulfill 
!
! ||||    h2>h1, p2>p1   ||||
 END_DOC
 integer :: i,h1,p1,h2,p2
 integer :: hh1,hh2,pp1,pp2,m,mm
 integer                        :: Ne(2)
 integer,           allocatable :: occ(:,:)
 double precision :: contrib
 allocate( occ(N_int*bit_kind_size,2) )
 call bitstring_to_list_ab(ref_bitmask,occ,Ne,N_int)
 call give_contrib_for_bbbb(1,1 ,1 ,1 ,occ,Ne,contrib)
 eff_2_e_from_3_e_bb = 100000000.d0
 !$OMP PARALLEL                                                                         &
 !$OMP DEFAULT (NONE)                                                                   &
 !$OMP PRIVATE (hh1, h1, hh2, h2, pp1, p1, pp2, p2, contrib) & 
 !$OMP SHARED (n_act_orb, list_act, Ne,occ, eff_2_e_from_3_e_bb)
 !$OMP DO SCHEDULE (static) 
  do hh1 = 1, n_act_orb !! beta 
    h1 = list_act(hh1) 
    do hh2 = hh1+1, n_act_orb !! beta
      h2 = list_act(hh2) 
      do pp1 = 1, n_act_orb !! beta
        p1 = list_act(pp1)
        do pp2 = pp1+1, n_act_orb !! beta
          p2 = list_act(pp2)
          call give_contrib_for_bbbb(h1,h2,p1,p2,occ,Ne,contrib)
          eff_2_e_from_3_e_bb(p2,p1,h2,h1) = contrib
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

END_PROVIDER 

subroutine give_contrib_for_bbbb(h1,h2,p1,p2,occ,Ne,contrib)
 implicit none
 BEGIN_DOC 
! gives the contribution for a double excitation (h1,p1)_beta (h2,p2)_beta
!
! on top of a determinant whose occupied orbitals is in (occ, Ne)
 END_DOC
 integer, intent(in) :: h1,h2,p1,p2,occ(N_int*bit_kind_size,2),Ne(2)
 double precision, intent(out) :: contrib
 integer :: mm,m 
 double precision :: direct_int, exchange_int
 !! h1,p1 == beta
 !! h2,p2 == beta
 contrib = 0.d0
 do mm = 1, Ne(2) !! beta ==> pure parallele spin contribution
   m = occ(mm,1)
   contrib += three_e_double_parrallel_spin_prov(m,p2,h2,p1,h1)
 enddo

 do mm = 1, Ne(1) !! alpha
   m = occ(mm,1)
   direct_int   = three_e_5_idx_direct_bi_ort(mm,p2,h2,p1,h1) 
   ! exchange between (h1,p1) and (h2,p2)
   exchange_int = three_e_5_idx_exch12_bi_ort(mm,p2,h2,p1,h1)
   contrib += direct_int - exchange_int
 enddo
end


subroutine double_htilde_mu_mat_fock_bi_ortho_no_3e(Nint, key_j, key_i, htot)

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
  double precision, intent(out) :: htot
  double precision :: hmono, htwoe
  integer                       :: occ(Nint*bit_kind_size,2)
  integer                       :: Ne(2), i, j, ii, jj, ispin, jspin, k, kk
  integer                       :: degree,exc(0:2,2,2)
  integer                       :: h1, p1, h2, p2, s1, s2
  double precision              :: get_mo_two_e_integral_tc_int,phase


  call get_excitation_degree(key_i, key_j, degree, Nint)

  hmono  = 0.d0
  htwoe  = 0.d0
  htot   = 0.d0

  if(degree.ne.2)then
   return
  endif
  integer :: degree_i,degree_j
  call get_excitation_degree(ref_bitmask,key_i,degree_i,N_int)
  call get_excitation_degree(ref_bitmask,key_j,degree_j,N_int)
  call get_double_excitation(key_i, key_j, exc, phase, Nint)
  call decode_exc(exc, 2, h1, p1, h2, p2, s1, s2)

  if(s1.ne.s2)then
   ! opposite spin two-body 
    htwoe  = mo_bi_ortho_tc_two_e(p2,p1,h2,h1) 
  else
   ! same spin two-body 
   ! direct terms 
   htwoe  = mo_bi_ortho_tc_two_e(p2,p1,h2,h1)  
   ! exchange terms 
   htwoe -= mo_bi_ortho_tc_two_e(p1,p2,h2,h1) 
  endif
  htwoe  *= phase
  htot    =  htwoe 

end

