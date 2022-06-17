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

subroutine generate_all_singles_cfg_with_type(bit_tmp,cfgInp,singles,idxs_singles,pq_singles,ex_type_singles,n_singles,Nint)
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
  integer*8, intent(in)          :: bit_tmp(0:N_configuration+1)
  integer, intent(in)            :: Nint
  integer, intent(inout)         :: n_singles
  integer, intent(out)           :: idxs_singles(*)
  integer, intent(out)           :: ex_type_singles(*)
  integer, intent(out)           :: pq_singles(2,*)
  integer(bit_kind), intent(in)  :: cfgInp(Nint,2)
  integer(bit_kind), intent(out) :: singles(Nint,2,*)
  integer(bit_kind)              :: Jdet(Nint,2)

  integer           :: i,k, n_singles_ma, i_hole, i_particle, ex_type, addcfg
  integer           :: ii,kk
  integer(bit_kind) :: single(Nint,2)
  logical           :: i_ok


  n_singles = 0
  !TODO
  !Make list of Somo  and Domo for holes
  !Make list of Unocc and Somo for particles
  !do i_hole = 1+n_core_orb, n_core_orb + n_act_orb
  do ii = 1, n_act_orb
    i_hole = list_act(ii)
    !do i_particle = 1+n_core_orb, n_core_orb + n_act_orb
    do kk = 1, n_act_orb
      i_particle = list_act(kk)
      if(i_hole .EQ. i_particle) cycle
      addcfg = -1
      call do_single_excitation_cfg_with_type(cfgInp,single,i_hole,i_particle,ex_type,i_ok)
      if (i_ok) then
        call binary_search_cfg(single,addcfg,bit_tmp)
        if(addcfg .EQ. -1) cycle
        n_singles = n_singles + 1
        do k=1,Nint
          singles(k,1,n_singles) = single(k,1)
          singles(k,2,n_singles) = single(k,2)
          ex_type_singles(n_singles) = ex_type
          pq_singles(1,n_singles) = i_hole     ! p
          pq_singles(2,n_singles) = i_particle ! q
          idxs_singles(n_singles) = addcfg
        enddo
      endif
    enddo
  enddo
end

