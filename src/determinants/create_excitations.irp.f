subroutine do_single_excitation(key_in,i_hole,i_particle,ispin,i_ok)
  implicit none
  BEGIN_DOC
  ! Apply the single excitation operator : a^{dager}_(i_particle) a_(i_hole) of spin = ispin
  ! on key_in
  ! ispin = 1  == alpha
  ! ispin = 2  == beta
  ! i_ok = 1  == the excitation is possible
  ! i_ok = -1 == the excitation is not possible
  END_DOC
  integer, intent(in)            :: i_hole,i_particle,ispin
  integer(bit_kind), intent(inout) :: key_in(N_int,2)
  integer, intent(out)           :: i_ok
  integer                        :: k,j,i
  integer(bit_kind)              :: mask
  use bitmasks
  ASSERT (i_hole > 0 )
  ASSERT (i_particle <= mo_num)
  i_ok = 1
  ! hole
  k = shiftr(i_hole-1,bit_kind_shift)+1
  j = i_hole-shiftl(k-1,bit_kind_shift)-1
  mask = ibset(0_bit_kind,j)
! check whether position j is occupied
  if (iand(key_in(k,ispin),mask) /= 0_bit_kind) then
   key_in(k,ispin) = ibclr(key_in(k,ispin),j)
  else 
   i_ok= -1
   return
  end if

  ! particle
  k = shiftr(i_particle-1,bit_kind_shift)+1
  j = i_particle-shiftl(k-1,bit_kind_shift)-1
  mask = ibset(0_bit_kind,j)
  if (iand(key_in(k,ispin),mask) == 0_bit_kind) then
   key_in(k,ispin) = ibset(key_in(k,ispin),j)
  else 
   i_ok= -1
   return
  end if

!  integer                        :: n_elec_tmp
!  n_elec_tmp = 0
!  do i = 1, N_int
!    n_elec_tmp += popcnt(key_in(i,1)) + popcnt(key_in(i,2))
!  enddo
!  if(n_elec_tmp .ne. elec_num)then
!    print*, n_elec_tmp,elec_num
!    call debug_det(key_in,N_int)
!    stop -1
!  endif
end


subroutine build_singly_excited_wavefunction(i_hole,i_particle,ispin,det_out,coef_out)
  implicit none
  BEGIN_DOC
  ! Applies the single excitation operator : a^{dager}_(i_particle) a_(i_hole) of
  ! spin = ispin to the current wave function (psi_det, psi_coef)
  END_DOC
  integer, intent(in)            :: i_hole,i_particle,ispin
  integer(bit_kind), intent(out) :: det_out(N_int,2,N_det)
  double precision, intent(out)  :: coef_out(N_det,N_states)

  integer :: k
  integer :: i_ok
  double precision :: phase
  do k=1,N_det
    coef_out(k,:) = psi_coef(k,:)
    det_out(:,:,k) = psi_det(:,:,k)
    call do_single_excitation(det_out(1,1,k),i_hole,i_particle,ispin,i_ok)
    if (i_ok == 1) then
      call get_phase(psi_det(1,1,k), det_out(1,1,k),phase,N_int)
      coef_out(k,:) = phase * coef_out(k,:)
    else
      coef_out(k,:) = 0.d0
      det_out(:,:,k) = psi_det(:,:,k)
    endif
  enddo
end

logical function is_spin_flip_possible(key_in,i_flip,ispin)
  implicit none
  BEGIN_DOC
  ! returns |true| if the spin-flip of spin ispin in the orbital i_flip is possible
  ! on key_in
  END_DOC
  integer, intent(in)            :: i_flip,ispin
  integer(bit_kind), intent(in)  :: key_in(N_int,2)
  integer                        :: k,j,i
  integer(bit_kind)              :: key_tmp(N_int,2)
  is_spin_flip_possible = .False.
  key_tmp = 0_bit_kind
  k = shiftr(i_flip-1,bit_kind_shift)+1
  j = i_flip-shiftl(k-1,bit_kind_shift)-1
  key_tmp(k,1) = ibset(key_tmp(k,1),j)
  integer                        :: other_spin(2)
  other_spin(1) = 2
  other_spin(2) = 1
  if(popcnt(iand(key_tmp(k,1),key_in(k,ispin))) == 1 .and. popcnt(iand(key_tmp(k,1),key_in(k,other_spin(ispin)))) == 0  )then
    ! There is a spin "ispin" in the orbital i_flip   AND  There is no electron of opposit spin in the same orbital "i_flip"
    is_spin_flip_possible = .True.
    return
  else
    return
  endif
end

