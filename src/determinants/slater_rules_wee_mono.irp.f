
subroutine i_Wee_j_single(key_i,key_j,Nint,spin,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$  where $i$ and $j$ are determinants differing by a
  ! single excitation.
  END_DOC
  integer, intent(in)            :: Nint, spin
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij

  integer                        :: exc(0:2,2)
  double precision               :: phase

  PROVIDE big_array_exchange_integrals mo_two_e_integrals_in_map

  call get_single_excitation_spin(key_i(1,spin),key_j(1,spin),exc,phase,Nint)
  call single_excitation_wee(key_i,key_j,exc(1,1),exc(1,2),spin,phase,hij)
end


double precision function diag_wee_mat_elem(det_in,Nint)
  implicit none
  BEGIN_DOC
  ! Computes $\langle i|H|i \rangle$.
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)

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

  diag_wee_mat_elem = ref_bitmask_two_e_energy
  if (nexc(1)+nexc(2) == 0) then
    return
  endif

  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(particle, occ_particle, tmp, Nint)
  ASSERT (tmp(1) == nexc(1))
  ASSERT (tmp(2) == nexc(2))
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(hole, occ_hole, tmp, Nint)
  ASSERT (tmp(1) == nexc(1))
  ASSERT (tmp(2) == nexc(2))

  det_tmp = ref_bitmask
  do ispin=1,2
    na = elec_num_tab(ispin)
    nb = elec_num_tab(iand(ispin,1)+1)
    do i=1,nexc(ispin)
      !DIR$ FORCEINLINE
      call ac_operator_two_e( occ_particle(i,ispin), ispin, det_tmp, diag_wee_mat_elem, Nint,na,nb)
      !DIR$ FORCEINLINE
      call a_operator_two_e ( occ_hole    (i,ispin), ispin, det_tmp, diag_wee_mat_elem, Nint,na,nb)
    enddo
  enddo
end


subroutine a_operator_two_e(iorb,ispin,key,hjj,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Needed for :c:func:`diag_Wee_mat_elem`.
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hjj

  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i
  integer                        :: tmp(2)

  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)

  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k > 0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  key(k,ispin) = ibclr(key(k,ispin),l)
  other_spin = iand(ispin,1)+1

  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(key, occ, tmp, Nint)
  na = na-1

  ! Same spin
  do i=1,na
    hjj = hjj - mo_two_e_integrals_jj_anti(occ(i,ispin),iorb)
  enddo

  ! Opposite spin
  do i=1,nb
    hjj = hjj - mo_two_e_integrals_jj(occ(i,other_spin),iorb)
  enddo

end


subroutine ac_operator_two_e(iorb,ispin,key,hjj,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Needed for :c:func:`diag_Wee_mat_elem`.
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hjj

  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i

  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)

  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(key, occ, tmp, Nint)
  ASSERT (tmp(1) == elec_alpha_num)
  ASSERT (tmp(2) == elec_beta_num)

  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k > 0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  key(k,ispin) = ibset(key(k,ispin),l)
  other_spin = iand(ispin,1)+1


  ! Same spin
  do i=1,na
    hjj = hjj + mo_two_e_integrals_jj_anti(occ(i,ispin),iorb)
  enddo

  ! Opposite spin
  do i=1,nb
    hjj = hjj + mo_two_e_integrals_jj(occ(i,other_spin),iorb)
  enddo
  na = na+1
end



subroutine i_H_j_mono_spin_one_e(key_i,key_j,Nint,spin,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$  where $i$ and $j$ are determinants differing by
  ! a single excitation.
  END_DOC
  integer, intent(in)            :: Nint, spin
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij

  integer                        :: exc(0:2,2)
  double precision               :: phase

  call get_single_excitation_spin(key_i(1,spin),key_j(1,spin),exc,phase,Nint)
  integer :: m,p
  m = exc(1,1)
  p = exc(1,2)
  hij = phase * mo_one_e_integrals(m,p)
end


double precision function diag_H_mat_elem_one_e(det_in,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $\langle i|H|i \rangle$.
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)

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

  diag_H_mat_elem_one_e = 0.d0

  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(det_in, occ_particle, tmp, Nint)
  do ispin = 1,2
   do i = 1, tmp(ispin)
    diag_H_mat_elem_one_e +=  mo_one_e_integrals(occ_particle(i,ispin),occ_particle(i,ispin))
   enddo
  enddo

end

subroutine i_H_j_one_e(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$  where $i$ and $j$ are determinants.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij

  integer :: degree,m,p
  double precision :: diag_H_mat_elem_one_e,phase
  integer                        :: exc(0:2,2,2)
  call get_excitation_degree(key_i,key_j,degree,Nint)
  hij = 0.d0
  if(degree>1)then
   return
  endif
  if(degree==0)then
   hij = diag_H_mat_elem_one_e(key_i,N_int)
  else
   call get_single_excitation(key_i,key_j,exc,phase,Nint)
   if (exc(0,1,1) == 1) then
     ! Mono alpha
     m = exc(1,1,1)
     p = exc(1,2,1)
   else
     ! Mono beta
     m = exc(1,1,2)
     p = exc(1,2,2)
   endif
   hij = phase * mo_one_e_integrals(m,p)
  endif

end

subroutine i_H_j_two_e(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$  where $i$ and $j$ are determinants.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij

  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: get_two_e_integral
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  double precision               :: diag_H_mat_elem, phase,phase_2
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals ref_bitmask_two_e_energy

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)

  hij = 0.d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    case (2)
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Mono alpha, mono beta
        if(exc(1,1,1) == exc(1,2,2) )then
         hij = phase * big_array_exchange_integrals(exc(1,1,1),exc(1,1,2),exc(1,2,1))
        else if (exc(1,2,1) ==exc(1,1,2))then
         hij = phase * big_array_exchange_integrals(exc(1,2,1),exc(1,1,1),exc(1,2,2))
        else
         hij = phase*get_two_e_integral(                          &
             exc(1,1,1),                                              &
             exc(1,1,2),                                              &
             exc(1,2,1),                                              &
             exc(1,2,2) ,mo_integrals_map)
        endif
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_two_e_integral(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_map) -                          &
            get_two_e_integral(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_two_e_integral(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_map) -                          &
            get_two_e_integral(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_map) )
      endif
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      if (exc(0,1,1) == 1) then
        ! Mono alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        spin = 1
      else
        ! Mono beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
      endif
      call single_excitation_wee(key_i,key_j,p,m,spin,phase,hij)
    case (0)
      double precision :: diag_wee_mat_elem
      hij = diag_wee_mat_elem(key_i,Nint)
  end select
end

