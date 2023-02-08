BEGIN_PROVIDER [ double precision, normal_two_body_bi_orth, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC 
  ! Normal ordering of the three body interaction on the HF density
  END_DOC 

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none

  integer :: i,h1,p1,h2,p2
  integer :: hh1,hh2,pp1,pp2
  integer                        :: Ne(2)
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)
  double precision :: hthree_aba,hthree_aaa,hthree_aab
  double precision :: wall0,wall1
 
  PROVIDE N_int

  allocate( occ(N_int*bit_kind_size,2) )
  allocate( key_i_core(N_int,2) )

  if(core_tc_op) then
    do i = 1, N_int
      key_i_core(i,1) = xor(ref_bitmask(i,1),core_bitmask(i,1))
      key_i_core(i,2) = xor(ref_bitmask(i,2),core_bitmask(i,2))
    enddo
    call bitstring_to_list_ab(key_i_core,occ,Ne,N_int)
  else
    call bitstring_to_list_ab(ref_bitmask,occ,Ne,N_int)
  endif

  normal_two_body_bi_orth = 0.d0
  print*,'Providing normal_two_body_bi_orth ...'
  call wall_time(wall0)

 !$OMP PARALLEL                                                                         &
 !$OMP DEFAULT (NONE)                                                                   &
 !$OMP PRIVATE (hh1, h1, hh2, h2, pp1, p1, pp2, p2, hthree_aba, hthree_aab, hthree_aaa) & 
 !$OMP SHARED (N_int, n_act_orb, list_act, Ne, occ, normal_two_body_bi_orth)
 !$OMP DO SCHEDULE (static) 
  do hh1 = 1, n_act_orb
    h1 = list_act(hh1) 
    do pp1 = 1, n_act_orb
      p1 = list_act(pp1)
      do hh2 = 1, n_act_orb
        h2 = list_act(hh2) 
        do pp2 = 1, n_act_orb
          p2 = list_act(pp2)
          ! opposite spin double excitations 
          call give_aba_contraction(N_int, h1, h2, p1, p2, Ne, occ, hthree_aba)
          ! same spin double excitations with opposite spin contributions 
          if(h1<h2.and.p1.gt.p2)then
           call give_aab_contraction(N_int, h2, h1, p1, p2, Ne, occ, hthree_aab) ! exchange h1<->h2
           ! same spin double excitations with same spin contributions 
           if(Ne(2).ge.3)then
             call give_aaa_contraction(N_int, h2, h1, p1, p2, Ne, occ, hthree_aaa) ! exchange h1<->h2
           else
             hthree_aaa = 0.d0
           endif
          else
           call give_aab_contraction(N_int, h1, h2, p1, p2, Ne, occ, hthree_aab)
           if(Ne(2).ge.3)then
             call give_aaa_contraction(N_int, h1, h2, p1, p2, Ne, occ, hthree_aaa)
           else
             hthree_aaa = 0.d0
           endif
          endif
          normal_two_body_bi_orth(p2,h2,p1,h1) = 0.5d0*(hthree_aba + hthree_aab + hthree_aaa)
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO
 !$OMP END PARALLEL

  call wall_time(wall1)
  print*,'Wall time for normal_two_body_bi_orth ',wall1-wall0

  deallocate( occ )
  deallocate( key_i_core )

END_PROVIDER 



subroutine give_aba_contraction(Nint, h1, h2, p1, p2, Ne, occ, hthree)

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer, intent(in)           :: Nint, h1, h2, p1, p2
  integer, intent(in)           :: Ne(2), occ(Nint*bit_kind_size,2)
  double precision, intent(out) :: hthree
  integer                       :: ii, i
  double precision              :: int_direct, int_exc_12, int_exc_13, integral

  !!!! double alpha/beta
  hthree = 0.d0
  do ii = 1, Ne(2) ! purely closed shell part 
    i = occ(ii,2)
    call give_integrals_3_body_bi_ort(i ,p2,p1,i,h2,h1,integral)
    int_direct = -1.d0 * integral
    call give_integrals_3_body_bi_ort(p1,p2, i,i,h2,h1,integral)
    int_exc_13 = -1.d0 * integral
    call give_integrals_3_body_bi_ort(p2, i,p1,i,h2,h1,integral)
    int_exc_12 = -1.d0 * integral
    hthree += 2.d0 * int_direct - 1.d0 * ( int_exc_13 + int_exc_12)
  enddo
  do ii = Ne(2) + 1, Ne(1) ! purely open-shell part 
   i = occ(ii,1)
    call give_integrals_3_body_bi_ort(i ,p2,p1,i,h2,h1,integral)
    int_direct = -1.d0 * integral
    call give_integrals_3_body_bi_ort(p1,p2, i,i,h2,h1,integral)
    int_exc_13 = -1.d0 * integral
    call give_integrals_3_body_bi_ort(p2, i,p1,i,h2,h1,integral)
    int_exc_12 = -1.d0 * integral
    hthree += 1.d0 * int_direct - 0.5d0* ( int_exc_13 + int_exc_12)
  enddo

end subroutine give_aba_contraction



BEGIN_PROVIDER [ double precision, normal_two_body_bi_orth_ab, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  ! Normal ordered two-body sector of the three-body terms for opposite spin double excitations 
  END_DOC

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: h1, p1, h2, p2, i
  integer                        :: hh1, hh2, pp1, pp2
  integer                        :: Ne(2)
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)
  double precision               :: hthree

  PROVIDE N_int

  allocate( key_i_core(N_int,2) )
  allocate( occ(N_int*bit_kind_size,2) )

  if(core_tc_op)then
   do i = 1, N_int
    key_i_core(i,1) = xor(ref_bitmask(i,1),core_bitmask(i,1))
    key_i_core(i,2) = xor(ref_bitmask(i,2),core_bitmask(i,2))
   enddo
   call bitstring_to_list_ab(key_i_core,occ,Ne,N_int)
  else
   call bitstring_to_list_ab(ref_bitmask,occ,Ne,N_int)
  endif
  normal_two_body_bi_orth_ab = 0.d0
  do hh1 = 1, n_act_orb
   h1 = list_act(hh1) 
   do pp1 = 1, n_act_orb
    p1 = list_act(pp1)
    do hh2 = 1, n_act_orb
     h2 = list_act(hh2) 
     do pp2 = 1, n_act_orb
      p2 = list_act(pp2)
      call give_aba_contraction(N_int, h1, h2, p1, p2, Ne, occ, hthree)
      normal_two_body_bi_orth_ab(p2,h2,p1,h1) = hthree    
     enddo
    enddo
   enddo
  enddo

  deallocate( key_i_core )
  deallocate( occ )

END_PROVIDER 



BEGIN_PROVIDER [ double precision, normal_two_body_bi_orth_aa_bb, (n_act_orb, n_act_orb, n_act_orb, n_act_orb)]

  BEGIN_DOC
  ! Normal ordered two-body sector of the three-body terms for same spin double excitations 
  END_DOC

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: i,ii,j,h1,p1,h2,p2
  integer                        :: hh1,hh2,pp1,pp2
  integer                        :: Ne(2)
  integer,           allocatable :: occ(:,:)
  integer(bit_kind), allocatable :: key_i_core(:,:)
  double precision               :: hthree_aab, hthree_aaa

  PROVIDE N_int

  allocate( key_i_core(N_int,2) )
  allocate( occ(N_int*bit_kind_size,2) )

  if(core_tc_op)then
    do i = 1, N_int
      key_i_core(i,1) = xor(ref_bitmask(i,1),core_bitmask(i,1))
      key_i_core(i,2) = xor(ref_bitmask(i,2),core_bitmask(i,2))
    enddo
    call bitstring_to_list_ab(key_i_core, occ, Ne, N_int)
  else
    call bitstring_to_list_ab(ref_bitmask, occ, Ne, N_int)
  endif

  normal_two_body_bi_orth_aa_bb = 0.d0
  do hh1 = 1, n_act_orb
    h1 = list_act(hh1) 
    do pp1 = 1 , n_act_orb
      p1 = list_act(pp1)
      do hh2 = 1, n_act_orb
        h2 = list_act(hh2) 
        do pp2 = 1 , n_act_orb
          p2 = list_act(pp2)
          if(h1<h2.and.p1.gt.p2)then
           call give_aab_contraction(N_int, h2, h1, p1, p2, Ne, occ, hthree_aab) ! exchange h1<->h2
           if(Ne(2).ge.3)then
             call give_aaa_contraction(N_int, h2, h1, p1, p2, Ne, occ, hthree_aaa) ! exchange h1<->h2
           else
             hthree_aaa = 0.d0
           endif
          else
           call give_aab_contraction(N_int, h1, h2, p1, p2, Ne, occ, hthree_aab)
           if(Ne(2).ge.3)then
             call give_aaa_contraction(N_int, h1, h2, p1, p2, Ne, occ, hthree_aaa)
           else
             hthree_aaa = 0.d0
           endif
          endif
          normal_two_body_bi_orth_aa_bb(p2,h2,p1,h1) = hthree_aab + hthree_aaa
        enddo
      enddo
    enddo
  enddo

  deallocate( key_i_core )
  deallocate( occ )

END_PROVIDER 



subroutine give_aaa_contraction(Nint, h1, h2, p1, p2, Ne, occ, hthree)

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer, intent(in)           :: Nint, h1, h2, p1, p2
  integer, intent(in)           :: Ne(2), occ(Nint*bit_kind_size,2)
  double precision, intent(out) :: hthree
  integer                       :: ii,i
  double precision              :: int_direct,int_exc_12,int_exc_13,int_exc_23
  double precision              :: integral,int_exc_l,int_exc_ll

  hthree = 0.d0
  do ii = 1, Ne(2) ! purely closed shell part 
    i = occ(ii,2)
    call give_integrals_3_body_bi_ort(i ,p2,p1,i,h2,h1,integral)
    int_direct = -1.d0 * integral
    call give_integrals_3_body_bi_ort(p2,p1,i ,i,h2,h1,integral)
    int_exc_l = -1.d0 * integral
    call give_integrals_3_body_bi_ort(p1,i ,p2,i,h2,h1,integral)
    int_exc_ll= -1.d0 * integral
    call give_integrals_3_body_bi_ort(p2,i ,p1,i,h2,h1,integral)
    int_exc_12= -1.d0 * integral
    call give_integrals_3_body_bi_ort(p1,p2, i,i,h2,h1,integral)
    int_exc_13= -1.d0 * integral
    call give_integrals_3_body_bi_ort(i ,p1,p2,i,h2,h1,integral)
    int_exc_23= -1.d0 * integral

    hthree +=  1.d0 * int_direct + int_exc_l + int_exc_ll -( int_exc_12+ int_exc_13+ int_exc_23  )
  enddo
  do ii = Ne(2)+1,Ne(1) ! purely open-shell part 
    i = occ(ii,1)
    call give_integrals_3_body_bi_ort(i ,p2,p1,i,h2,h1,integral)
    int_direct = -1.d0 * integral
    call give_integrals_3_body_bi_ort(p2,p1,i ,i,h2,h1,integral)
    int_exc_l = -1.d0 * integral
    call give_integrals_3_body_bi_ort(p1,i ,p2,i,h2,h1,integral)
    int_exc_ll= -1.d0 * integral
    call give_integrals_3_body_bi_ort(p2,i ,p1,i,h2,h1,integral)
    int_exc_12= -1.d0 * integral
    call give_integrals_3_body_bi_ort(p1,p2, i,i,h2,h1,integral)
    int_exc_13= -1.d0 * integral
    call give_integrals_3_body_bi_ort(i ,p1,p2,i,h2,h1,integral)
    int_exc_23= -1.d0 * integral

    hthree +=  1.d0 * int_direct + 0.5d0 * (int_exc_l + int_exc_ll -( int_exc_12+ int_exc_13+ int_exc_23  ))
  enddo

end subroutine give_aaa_contraction



subroutine give_aab_contraction(Nint, h1, h2, p1, p2, Ne, occ, hthree)
  implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
  integer, intent(in)           :: Nint, h1, h2, p1, p2
  integer, intent(in)           :: Ne(2), occ(Nint*bit_kind_size,2)
  double precision, intent(out) :: hthree
  integer                       :: ii, i
  double precision              :: int_direct, int_exc_12, int_exc_13, int_exc_23
  double precision              :: integral, int_exc_l, int_exc_ll

  hthree = 0.d0
  do ii = 1, Ne(2) ! purely closed shell part 
    i = occ(ii,2)
    call give_integrals_3_body_bi_ort(p2,p1,i,h2,h1,i,integral)
    int_direct = -1.d0 * integral
    call give_integrals_3_body_bi_ort(p1,p2,i,h2,h1,i,integral)
    int_exc_23= -1.d0 * integral
    hthree  +=  1.d0 * int_direct - int_exc_23
  enddo

end subroutine give_aab_contraction
