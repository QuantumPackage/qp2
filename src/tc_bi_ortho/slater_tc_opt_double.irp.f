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
 eff_2_e_from_3_e_ab = 0.d0
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
 eff_2_e_from_3_e_aa = 100000000.d0
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
 double precision :: three_e_double_parrallel_spin
 !! h1,p1 == alpha 
 !! h2,p2 == alpha
 contrib = 0.d0
 do mm = 1, Ne(1) !! alpha ==> pure parallele spin contribution
   m = occ(mm,1)
   contrib += three_e_double_parrallel_spin(m,p2,h2,p1,h1)
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
 eff_2_e_from_3_e_bb = 100000000.d0
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
 double precision :: three_e_double_parrallel_spin
 !! h1,p1 == beta
 !! h2,p2 == beta
 contrib = 0.d0
 do mm = 1, Ne(2) !! beta ==> pure parallele spin contribution
   m = occ(mm,1)
   contrib += three_e_double_parrallel_spin(m,p2,h2,p1,h1)
 enddo

 do mm = 1, Ne(1) !! alpha
   m = occ(mm,1)
   direct_int   = three_e_5_idx_direct_bi_ort(mm,p2,h2,p1,h1) 
   ! exchange between (h1,p1) and (h2,p2)
   exchange_int = three_e_5_idx_exch12_bi_ort(mm,p2,h2,p1,h1)
   contrib += direct_int - exchange_int
 enddo
end

