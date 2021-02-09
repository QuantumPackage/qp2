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
    ! There is a spin "ispin" in the orbital i_flip   AND
    ! There is no electron of opposit spin in the same orbital "i_flip"
    is_spin_flip_possible = .True.
    return
  else
    return
  endif
end

subroutine do_single_excitation_cfg(key_in,key_out,i_hole,i_particle,ok)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Applies the single excitation operator to a configuration
  ! If the excitation is possible, ok is True
  END_DOC
  integer, intent(in)            :: i_hole,i_particle
  integer(bit_kind), intent(in)  :: key_in(N_int,2)
  logical , intent(out)          :: ok
  integer                        :: k,j,i
  integer(bit_kind)              :: mask
  integer(bit_kind)              :: key_out(N_int,2)

  ASSERT (i_hole > 0)
  ASSERT (i_particle <= mo_num)

  ok = .True.
  key_out(:,:) = key_in(:,:)

  ! hole
  k = shiftr(i_hole-1,bit_kind_shift)+1
  j = i_hole-shiftl(k-1,bit_kind_shift)-1
  mask = ibset(0_bit_kind,j)

  ! Check if the position j is singly occupied
  ! 1  ->  0  (SOMO)
  ! 0      0  (DOMO)
  if (iand(key_out(k,1),mask) /= 0_bit_kind) then
     key_out(k,1) = ibclr(key_out(k,1),j)

  ! Check if the position j is doubly occupied
  ! 0  ->  1  (SOMO)
  ! 1      0  (DOMO)
  else if (iand(key_out(k,2),mask) /= 0_bit_kind) then
     key_out(k,1) = ibset(key_out(k,1),j)
     key_out(k,2) = ibclr(key_out(k,2),j)

  ! The position j is unoccupied: Not OK
  ! 0  ->  0  (SOMO)
  ! 0      0  (DOMO)
  else
     ok =.False.
     return
  endif


  ! particle
  k = shiftr(i_particle-1,bit_kind_shift)+1
  j = i_particle-shiftl(k-1,bit_kind_shift)-1
  mask = ibset(0_bit_kind,j)

  ! Check if the position j is singly occupied
  ! 1  ->  0  (SOMO)
  ! 0      1  (DOMO)
  if (iand(key_out(k,1),mask) /= 0_bit_kind) then
     key_out(k,1) = ibclr(key_out(k,1),j)
     key_out(k,2) = ibset(key_out(k,2),j)

  ! Check if the position j is doubly occupied : Not OK
  ! 0  ->  1  (SOMO)
  ! 1      0  (DOMO)
  else if (iand(key_out(k,2),mask) /= 0_bit_kind) then
     ok = .False.
     return

  ! Position at j is unoccupied
  ! 0  ->  0  (SOMO)
  ! 0      0  (DOMO)
  else
     key_out(k,1) = ibset(key_out(k,1),j)
  endif

end

subroutine do_single_excitation_cfg_with_type(key_in,key_out,i_hole,i_particle,ex_type,ok)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Applies the single excitation operator to a configuration
  ! Returns the type of excitation in ex_type
  ! where the following convention is used
  ! 1 = (SOMO -> SOMO) 1 change in Nsomo
  ! 2 = (DOMO ->  VMO) 1 change in Nsomo
  ! 3 = (SOMO ->  VMO) 0 change in Nsomo
  ! 4 = (DOMO -> SOMO) 0 change in Nsomo
  ! If the excitation is possible, ok is True
  END_DOC
  integer, intent(in)            :: i_hole,i_particle
  integer(bit_kind), intent(in)  :: key_in(N_int,2)
  integer , intent(out)          :: ex_type
  logical , intent(out)          :: ok
  integer                        :: k,j,i
  integer(bit_kind)              :: mask
  integer(bit_kind)              :: key_out(N_int,2)
  logical                        :: isholeSOMO
  logical                        :: isparticleSOMO
  logical                        :: isholeDOMO
  logical                        :: isparticleVMO
  isholeSOMO = .False.
  isholeDOMO = .False.
  isparticleSOMO = .False.
  isparticleVMO = .False.

  ASSERT (i_hole > 0)
  ASSERT (i_particle <= mo_num)

  ok = .True.
  key_out(:,:) = key_in(:,:)

  ! hole
  k = shiftr(i_hole-1,bit_kind_shift)+1
  j = i_hole-shiftl(k-1,bit_kind_shift)-1
  mask = ibset(0_bit_kind,j)

  ! Check if the position j is singly occupied
  ! 1  ->  0  (SOMO)
  ! 0      0  (DOMO)
  if (iand(key_out(k,1),mask) /= 0_bit_kind) then
     key_out(k,1) = ibclr(key_out(k,1),j)
     isholeSOMO = .True.

  ! Check if the position j is doubly occupied
  ! 0  ->  1  (SOMO)
  ! 1      0  (DOMO)
  else if (iand(key_out(k,2),mask) /= 0_bit_kind) then
     key_out(k,1) = ibset(key_out(k,1),j)
     key_out(k,2) = ibclr(key_out(k,2),j)
     isholeDOMO = .True.

  ! The position j is unoccupied: Not OK
  ! 0  ->  0  (SOMO)
  ! 0      0  (DOMO)
  else
     ok =.False.
     return
  endif


  ! particle
  k = shiftr(i_particle-1,bit_kind_shift)+1
  j = i_particle-shiftl(k-1,bit_kind_shift)-1
  mask = ibset(0_bit_kind,j)

  ! Check if the position j is singly occupied
  ! 1  ->  0  (SOMO)
  ! 0      1  (DOMO)
  if (iand(key_out(k,1),mask) /= 0_bit_kind) then
     key_out(k,1) = ibclr(key_out(k,1),j)
     key_out(k,2) = ibset(key_out(k,2),j)
     isparticleSOMO = .True.

  ! Check if the position j is doubly occupied : Not OK
  ! 0  ->  1  (SOMO)
  ! 1      0  (DOMO)
  else if (iand(key_out(k,2),mask) /= 0_bit_kind) then
     ok = .False.
     return

  ! Position at j is unoccupied
  ! 0  ->  0  (SOMO)
  ! 0      0  (DOMO)
  else
     key_out(k,1) = ibset(key_out(k,1),j)
     isparticleVMO = .True.
  endif

  if(isholeSOMO) then
      ! two possibilities
      ! particle is SOMO or VMO
      if(isparticleSOMO) then
          ! SOMO -> SOMO
          ex_type = 1
      else
          ! SOMO -> VMO
          ex_type = 3
      endif 
   else
       ! two possibilities
       ! particle is SOMO or VMO
       if(isparticleSOMO) then
           ! DOMO -> SOMO
           ex_type = 4
       else
           ! DOMO -> VMO
           ex_type = 2
       endif
    endif

end

subroutine generate_all_singles_cfg(cfg,singles,n_singles,Nint)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Generate all single excitation wrt a configuration
  !
  ! n_singles : on input, max number of singles :
  !   elec_alpha_num * (mo_num - elec_beta_num)
  !   on output, number of generated singles
  END_DOC
  integer, intent(in)            :: Nint
  integer, intent(inout)         :: n_singles
  integer(bit_kind), intent(in)  :: cfg(Nint,2)
  integer(bit_kind), intent(out) :: singles(Nint,2,*)

  integer           :: i,k, n_singles_ma, i_hole, i_particle
  integer(bit_kind) :: single(Nint,2)
  logical           :: i_ok

  n_singles = 0
  !TODO
  !Make list of Somo  and Domo for holes
  !Make list of Unocc and Somo for particles
  do i_hole = 1, mo_num
    do i_particle = 1, mo_num
      call do_single_excitation_cfg(cfg,single,i_hole,i_particle,i_ok)
      if (i_ok) then
        n_singles = n_singles + 1
        do k=1,Nint
          singles(k,1,n_singles) = single(k,1)
          singles(k,2,n_singles) = single(k,2)
        enddo
      endif
    enddo
  enddo
end

subroutine generate_all_singles_cfg_with_type(cfg,singles,idxs_singles,pq_singles,ex_type_singles,n_singles,Nint)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Generate all single excitation wrt a configuration
  !
  ! n_singles : on input, max number of singles :
  !   elec_alpha_num * (mo_num - elec_beta_num)
  !   on output, number of generated singles
  ! ex_type_singles : on output contains type of excitations  :
  !
  END_DOC
  integer, intent(in)            :: Nint
  integer, intent(inout)         :: n_singles
  integer, intent(out)           :: idxs_singles(*)
  integer, intent(out)           :: ex_type_singles(*)
  real*8 , intent(out)           :: pq_singles(2,*)
  integer(bit_kind), intent(in)  :: cfg(Nint,2)
  integer(bit_kind), intent(out) :: singles(Nint,2,*)

  integer           :: i,k, n_singles_ma, i_hole, i_particle, ex_type
  integer(bit_kind) :: single(Nint,2)
  logical           :: i_ok

  n_singles = 0
  !TODO
  !Make list of Somo  and Domo for holes
  !Make list of Unocc and Somo for particles
  do i_hole = 1, mo_num
    do i_particle = 1, mo_num
      call do_single_excitation_cfg_with_type(cfg,single,i_hole,i_particle,ex_type,i_ok)
      if (i_ok) then
        n_singles = n_singles + 1
        do k=1,Nint
          singles(k,1,n_singles) = single(k,1)
          singles(k,2,n_singles) = single(k,2)
          ex_type_singles(n_singles) = ex_type
          pq_singles(1,n_singles) = i_hole     ! p
          pq_singles(1,n_singles) = i_particle ! q
        enddo
      endif
    enddo
  enddo
end

