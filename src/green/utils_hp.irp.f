subroutine print_mo_energies(key_ref,nint,nmo)
  use bitmasks
  BEGIN_DOC
  ! get mo energies for one det
  END_DOC
  implicit none
  integer, intent(in) :: nint, nmo
  integer(bit_kind), intent(in) :: key_ref(nint,2)
  double precision, allocatable :: e_mo(:,:)
  integer, allocatable :: occ(:,:),virt(:,:) !(nint*bit_kind_size,2)
  integer :: n_occ(2), n_virt(2)
  integer, parameter :: int_spin2(1:2) = (/2,1/)
  integer :: i,j,ispin,jspin,i0,j0,k
  integer(bit_kind), allocatable :: key_virt(:,:)
  integer, allocatable :: is_occ(:,:)


  allocate(occ(nint*bit_kind_size,2),virt(nint*bit_kind_size,2),key_virt(nint,2),e_mo(nmo,2),is_occ(nmo,2))
  is_occ=0

  call bitstring_to_list_ab(key_ref,occ,n_occ,nint)
  do i=1,nint
    do ispin=1,2
      key_virt(i,ispin)=xor(full_ijkl_bitmask(i),key_ref(i,ispin))
    enddo
  enddo
  call bitstring_to_list_ab(key_virt,virt,n_virt,nint)

  e_mo(1:nmo,1)=mo_one_e_integrals_diag(1:nmo)
  e_mo(1:nmo,2)=mo_one_e_integrals_diag(1:nmo)

  do ispin=1,2
    jspin=int_spin2(ispin)
    do i0=1,n_occ(ispin)
      i=occ(i0,ispin)
      is_occ(i,ispin)=1
      do j0=i0+1,n_occ(ispin)
        j=occ(j0,ispin)
        e_mo(i,ispin) = e_mo(i,ispin) + mo_two_e_integrals_jj_anti(i,j)
        e_mo(j,ispin) = e_mo(j,ispin) + mo_two_e_integrals_jj_anti(i,j)
      enddo
      do k=2,ispin
        do j0=1,n_occ(jspin)
          j=occ(j0,jspin)
          e_mo(i,ispin) = e_mo(i,ispin) + mo_two_e_integrals_jj(i,j)
          e_mo(j,jspin) = e_mo(j,jspin) + mo_two_e_integrals_jj(i,j) !can delete this and remove k level of loop
        enddo
      enddo
      do j0=1,n_virt(ispin)
        j=virt(j0,ispin)
        e_mo(j,ispin) = e_mo(j,ispin) + mo_two_e_integrals_jj_anti(i,j)
      enddo
      do j0=1,n_virt(jspin)
        j=virt(j0,jspin)
        e_mo(j,jspin) = e_mo(j,jspin) + mo_two_e_integrals_jj(i,j)
      enddo
    enddo
  enddo

  do i=1,nmo
    write(6,'(2(I5),2(E25.15))')is_occ(i,1),is_occ(i,2),e_mo(i,1),e_mo(i,2)
  enddo
  deallocate(occ,virt,key_virt,e_mo,is_occ)
end

subroutine get_mo_energies(key_ref,nint,nmo,e_mo)
  use bitmasks
  BEGIN_DOC
  ! get mo energies for one det
  END_DOC
  implicit none
  integer, intent(in) :: nint, nmo
  integer(bit_kind), intent(in) :: key_ref(nint,2)
  double precision, intent(out) :: e_mo(nmo,2)
  integer, allocatable :: occ(:,:),virt(:,:) !(nint*bit_kind_size,2)
  integer :: n_occ(2), n_virt(2)
  integer, parameter :: int_spin2(1:2) = (/2,1/)
  integer :: i,j,ispin,jspin,i0,j0,k
  integer(bit_kind), allocatable :: key_virt(:,:)


  allocate(occ(nint*bit_kind_size,2),virt(nint*bit_kind_size,2),key_virt(nint,2))

  call bitstring_to_list_ab(key_ref,occ,n_occ,nint)
  do i=1,nint
    do ispin=1,2
      key_virt(i,ispin)=xor(full_ijkl_bitmask(i),key_ref(i,ispin))
    enddo
  enddo
  call bitstring_to_list_ab(key_virt,virt,n_virt,nint)

  e_mo(1:nmo,1)=mo_one_e_integrals_diag(1:nmo)
  e_mo(1:nmo,2)=mo_one_e_integrals_diag(1:nmo)

  do ispin=1,2
    jspin=int_spin2(ispin)
    do i0=1,n_occ(ispin)
      i=occ(i0,ispin)
      do j0=i0+1,n_occ(ispin)
        j=occ(j0,ispin)
        e_mo(i,ispin) = e_mo(i,ispin) + mo_two_e_integrals_jj_anti(i,j)
        e_mo(j,ispin) = e_mo(j,ispin) + mo_two_e_integrals_jj_anti(i,j)
      enddo
      do k=2,ispin
        do j0=1,n_occ(jspin)
          j=occ(j0,jspin)
          e_mo(i,ispin) = e_mo(i,ispin) + mo_two_e_integrals_jj(i,j)
          e_mo(j,jspin) = e_mo(j,jspin) + mo_two_e_integrals_jj(i,j) !can delete this and remove k level of loop
        enddo
      enddo
      do j0=1,n_virt(ispin)
        j=virt(j0,ispin)
        e_mo(j,ispin) = e_mo(j,ispin) + mo_two_e_integrals_jj_anti(i,j)
      enddo
      do j0=1,n_virt(jspin)
        j=virt(j0,jspin)
        e_mo(j,jspin) = e_mo(j,jspin) + mo_two_e_integrals_jj(i,j)
      enddo
    enddo
  enddo

  deallocate(occ,virt,key_virt)
end

subroutine get_mask_phase_new(det1, pm, Nint)
  use bitmasks
  BEGIN_DOC
  ! phasemask copied from qp2
  ! return phasemask of det1 in pm
  END_DOC
  implicit none
  integer, intent(in) :: Nint
  integer(bit_kind), intent(in) :: det1(Nint,2)
  integer(bit_kind), intent(out) :: pm(Nint,2)
  integer(bit_kind) :: tmp1, tmp2
  integer :: i
  pm(1:Nint,1:2) = det1(1:Nint,1:2)
  tmp1 = 0_8
  tmp2 = 0_8
  do i=1,Nint
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 1))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 1))
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 2))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 2))
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 4))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 4))
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 8))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 8))
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 16))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 16))
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 32))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 32))
    pm(i,1) = ieor(pm(i,1), tmp1)
    pm(i,2) = ieor(pm(i,2), tmp2)
    if(iand(popcnt(det1(i,1)), 1) == 1) tmp1 = not(tmp1)
    if(iand(popcnt(det1(i,2)), 1) == 1) tmp2 = not(tmp2)
  end do
end subroutine

subroutine get_phase_hp(g_idx_int,g_idx_bit,g_spin,g_sign,det_in,g_det_phase,nint,n_g)
  use bitmasks
  implicit none
  integer, intent(in) :: nint,n_g
  integer, intent(in) :: g_idx_int(n_g), g_idx_bit(n_g),g_spin(n_g)
  double precision, intent(in) :: g_sign(n_g)
  integer(bit_kind), intent(in) :: det_in(nint,2)
  double precision, intent(out) :: g_det_phase(n_g)

  integer(bit_kind) :: tmp_spindet(nint), pm(nint,2)
  double precision, parameter :: phase_dble(0:1) = (/1.d0,-1.d0/)

  integer :: i
  logical :: is_allowed(n_g), all_banned, is_filled

  all_banned=.True.
  do i=1,n_g
    tmp_spindet(1:nint) = det_in(1:nint,g_spin(i))
    call spinorb_is_filled_int_bit(tmp_spindet,g_idx_int(i),g_idx_bit(i),nint,is_filled)
    is_allowed(i) = (.not.(((g_sign(i)<0).and.(.not.is_filled)).or.((g_sign(i)>0).and.(is_filled))))
    all_banned=(all_banned.and.(.not.is_allowed(i)))
  enddo

  if (all_banned) then
    g_det_phase(:)=0.d0
  else
    call get_mask_phase_new(det_in,pm,nint)
    do i=1,n_g
      if (is_allowed(i)) then
        g_det_phase(i) = phase_dble(popcnt(iand(ibset(0_bit_kind,g_idx_bit(i)),pm(g_idx_int(i),g_spin(i)))))
      else
        g_det_phase(i)=0.d0
      endif
    enddo
  endif
end

subroutine get_homo_lumo(key_ref,nint,nmo,idx_homo_lumo,spin_homo_lumo)
  use bitmasks
  implicit none
  integer, intent(in) :: nint,nmo
  integer(bit_kind), intent(in) :: key_ref(nint,2)
  integer, intent(out) :: idx_homo_lumo(2), spin_homo_lumo(2)

  double precision, allocatable :: e_mo(:,:)
  integer, allocatable :: occ(:,:),virt(:,:) !(nint*bit_kind_size,2)
  integer :: n_occ(2), n_virt(2)
  integer :: i,i0,ispin
  integer(bit_kind), allocatable :: key_virt(:,:)
  double precision :: maxocc(2), minvirt(2)
  integer :: imaxocc(2), iminvirt(2)

  allocate(e_mo(nmo,2),key_virt(nint,2),occ(nint*bit_kind_size,2),virt(nint*bit_kind_size,2))

  call get_mo_energies(key_ref,nint,nmo,e_mo)
  
  !allocate(occ(nint*bit_kind_size,2),virt(nint*bit_kind_size,2))

  call bitstring_to_list_ab(key_ref,occ,n_occ,nint)
  do i=1,nint
    do ispin=1,2
      key_virt(i,ispin)=xor(full_ijkl_bitmask(i),key_ref(i,ispin))
    enddo
  enddo
  call bitstring_to_list_ab(key_virt,virt,n_virt,nint)
  
  maxocc=-1.d20 !maybe use -1.d0*huge(0.d0)?
  minvirt=1.d20
  imaxocc=-1
  iminvirt=-1

  do ispin=1,2
    do i0=1,n_occ(ispin)
      i=occ(i0,ispin)
      if (e_mo(i,ispin).gt.maxocc(ispin)) then
        maxocc(ispin)=e_mo(i,ispin)
        imaxocc(ispin)=i
      endif
    enddo
    do i0=1,n_virt(ispin)
      i=virt(i0,ispin)
      if (e_mo(i,ispin).lt.minvirt(ispin)) then
        minvirt(ispin)=e_mo(i,ispin)
        iminvirt(ispin)=i
      endif
    enddo
  enddo
  double precision :: e_mo_thresh
  e_mo_thresh = 1.d-8
  !these should both just be 2x2 arrays, but performance here doesn't really matter and this is more readable
  !if (maxocc(1).ge.maxocc(2)) then
  if ((maxocc(2)-maxocc(1)).le.e_mo_thresh) then
    spin_homo_lumo(1)=1
  else
    spin_homo_lumo(1)=2
  endif
  if ((minvirt(1)-minvirt(2)).le.e_mo_thresh) then
    spin_homo_lumo(2)=1
  else
    spin_homo_lumo(2)=2
  endif

  idx_homo_lumo(1)=imaxocc(spin_homo_lumo(1))
  idx_homo_lumo(2)=iminvirt(spin_homo_lumo(2))

  deallocate(e_mo,occ,virt,key_virt)

end

subroutine get_list_hp_banned_ab(tmp_det,N_hp,exc_is_banned,spin_hp,sign_hp,idx_hp,nint,all_banned)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! input determinant tmp_det and list of single holes/particles
  ! for each hole/particle, determine whether it is filled/empty in tmp_det
  ! return which are disallowed in exc_is_banned
  ! if all are banned, set all_banned to true
  END_DOC
  integer, intent(in) :: N_hp,nint
  integer, intent(in) :: spin_hp(N_hp), idx_hp(N_hp)
  double precision, intent(in) :: sign_hp(N_hp)
  integer(bit_kind), intent(in) :: tmp_det(nint,2)
  logical, intent(out) :: exc_is_banned(N_hp)
  logical, intent(out) :: all_banned
  
  integer :: i
  logical :: is_filled

  all_banned = .True.
  do i=1,N_hp
    call orb_is_filled(tmp_det,idx_hp(i),spin_hp(i),nint,is_filled)
    if (sign_hp(i).gt.0) then ! particle creation, banned if already filled
      exc_is_banned(i) = is_filled
    else ! hole creation, banned if already empty
      exc_is_banned(i) = (.not.is_filled)
    endif
    all_banned = (all_banned.and.exc_is_banned(i))
  enddo
end

subroutine get_list_hp_banned_single_spin(tmp_spindet,N_hp,exc_is_banned,spin_hp,sign_hp,idx_hp,ispin,nint,all_banned)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! input spindeterminant tmp_spindet and list of single holes/particles
  ! tmp_spindet is only one spin part of a full det, with spin==ispin
  ! for each hole/particle, determine whether it is filled/empty in tmp_det
  ! return which are disallowed in exc_is_banned
  ! if all are banned, set all_banned to true
  END_DOC
  integer, intent(in) :: N_hp, ispin, nint
  integer, intent(in) :: spin_hp(N_hp), idx_hp(N_hp)
  double precision, intent(in) :: sign_hp(N_hp)
  integer(bit_kind), intent(in) :: tmp_spindet(nint)
  logical, intent(out) :: exc_is_banned(N_hp)
  logical, intent(out) :: all_banned
  
  integer :: i
  logical :: is_filled

  all_banned = .True.
  do i=1,N_hp
    if (spin_hp(i).eq.ispin) then
      call orb_is_filled_single_spin(tmp_spindet,idx_hp(i),nint,is_filled)
      if (sign_hp(i).gt.0) then ! particle creation, banned if already filled
        exc_is_banned(i) = is_filled
      else ! hole creation, banned if already empty
        exc_is_banned(i) = (.not.is_filled)
      endif
    else
      exc_is_banned(i) = .False.
    endif
    all_banned = (all_banned.and.exc_is_banned(i))
  enddo
end

subroutine get_list_hp_banned_spin(tmp_det,N_hp,exc_is_banned,spin_hp,sign_hp,idx_hp,ispin,nint,all_banned)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! input determinant tmp_det and list of single holes/particles
  ! for each hole/particle, determine whether it is filled/empty in tmp_det
  ! return which are disallowed in exc_is_banned
  ! if all are banned, set all_banned to true
  ! only consider tmp_det(1:N_int, ispin)
  END_DOC
  integer, intent(in) :: N_hp, ispin, nint
  integer, intent(in) :: spin_hp(N_hp), idx_hp(N_hp)
  double precision, intent(in) :: sign_hp(N_hp)
  integer(bit_kind), intent(in) :: tmp_det(nint,2)
  logical, intent(out) :: exc_is_banned(N_hp)
  logical, intent(out) :: all_banned

  integer(bit_kind) :: spindet(nint)
  
  integer :: i
  logical :: is_filled
  spindet(1:nint) = tmp_det(1:nint,ispin)

  all_banned = .True.
  do i=1,N_hp
    if (spin_hp(i).eq.ispin) then
      call orb_is_filled(tmp_det,idx_hp(i),ispin,nint,is_filled)
      if (sign_hp(i).gt.0) then ! particle creation, banned if already filled
        exc_is_banned(i) = is_filled
      else ! hole creation, banned if already empty
        exc_is_banned(i) = (.not.is_filled)
      endif
    else
      exc_is_banned(i) = .False.
    endif
    all_banned = (all_banned.and.exc_is_banned(i))
  enddo
end

  
subroutine spinorb_is_filled_int_bit(key_ref,iorb_int,iorb_bit,Nint,is_filled)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! determine whether iorb is filled in key_ref
  ! iorb is specified by int and bit locations within the determinant
  END_DOC
  integer, intent(in)            :: iorb_int, iorb_bit, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint)
  logical, intent(out) :: is_filled
  
  ASSERT (iorb_int > 0)
  ASSERT (iorb_bit >= 0)
  ASSERT (Nint > 0)
  is_filled = btest(key_ref(iorb_int),iorb_bit)
end

subroutine orb_is_filled_int_bit(key_ref,iorb_int,iorb_bit,ispin,Nint,is_filled)
  use bitmasks
  implicit none
  BEGIN_DOC
  !  todo: not finished
  ! determine whether iorb is filled in key_ref
  ! iorb is specified by int and bit locations within the determinant
  END_DOC
  integer, intent(in)            :: iorb_int, iorb_bit, ispin, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint,2)
  logical, intent(out) :: is_filled
  
  ASSERT (iorb_int > 0)
  ASSERT (iorb_bit >= 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  is_filled = btest(key_ref(iorb_int,ispin),iorb_bit)
!  call spinorb_is_filled_int_bit(key_ref(1,ispin),iorb_int,iorb_bit,Nint,is_filled)
end

subroutine get_orb_int_bit(iorb,iorb_int,iorb_bit)
  BEGIN_DOC
  ! get int and bit corresponding to orbital index iorb 
  END_DOC
  use bitmasks
  implicit none
  integer, intent(in)            :: iorb
  integer, intent(out)            :: iorb_int, iorb_bit
  ASSERT (iorb > 0)
  iorb_int = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (iorb_int > 0)
  iorb_bit = iorb - ishft(iorb_int-1,bit_kind_shift)-1
  ASSERT (iorb_bit >= 0)
end

subroutine orb_is_filled_single_spin(key_ref,iorb,Nint,is_filled)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! determine whether iorb is filled in key_ref
  ! key_ref is single alpha or beta determinant
  END_DOC
  integer, intent(in)            :: iorb, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint)
  logical, intent(out) :: is_filled
  
  integer                        :: k,l
  
  ASSERT (iorb > 0)
  ASSERT (Nint > 0)
  
  ! k is index of the int where iorb is found
  ! l is index of the bit where iorb is found
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  is_filled = btest(key_ref(k),l)  
end

subroutine orb_is_filled(key_ref,iorb,ispin,Nint,is_filled)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! determine whether iorb, ispin is filled in key_ref
  ! key_ref has alpha and beta parts
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint,2)
  logical, intent(out) :: is_filled
  
  integer                        :: k,l
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  ! k is index of the int where iorb is found
  ! l is index of the bit where iorb is found
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  is_filled = btest(key_ref(k,ispin),l)  
end

subroutine ac_operator_phase(key_new,key_ref,iorb,ispin,Nint,phase)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! apply creation operator to key_ref
  ! add electron with spin ispin to orbital with index iorb
  ! output resulting det and phase in key_new and phase
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint,2)
  integer(bit_kind), intent(out) :: key_new(Nint,2)
  double precision, intent(out) :: phase
  
  integer                        :: k,l,i

  double precision, parameter :: p(0:1) = (/ 1.d0, -1.d0 /)
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  key_new=key_ref

  ! alpha det is list of Nint 64-bit ints
  ! k is index of the int where iorb is found
  ! l is index of the bit where iorb is found
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  key_new(k,ispin) = ibset(key_new(k,ispin),l)
  
  integer(bit_kind) :: parity_filled

  ! I assume here that the ordering is all alpha spinorbs and then all beta spinorbs
  ! If we add an alpha electron, parity is not affected by beta part of determinant
  !   (only need number of alpha occupied orbs below iorb)
 
  ! If we add a beta electron, the parity is affected by alpha part
  !    (need total number of occupied alpha orbs (all of which come before beta)
  !      and total number of beta occupied orbs below iorb)

  if (ispin==1) then
    parity_filled=0_bit_kind
  else
    parity_filled=iand(int(elec_alpha_num,bit_kind),1_bit_kind)
  endif

  ! get parity due to orbs in other ints (with lower indices)
  do i=1,k-1
    parity_filled = iand(int(popcnt(key_ref(i,ispin)),bit_kind),parity_filled)
  enddo
  
  ! get parity due to orbs in same int as iorb
  ! ishft(1_bit_kind,l)-1 has its l rightmost bits set to 1, other bits set to 0
  parity_filled = iand(int(popcnt(iand(ishft(1_bit_kind,l)-1,key_ref(k,ispin))),bit_kind),parity_filled)
  phase = p(iand(1_bit_kind,parity_filled))

end

subroutine a_operator_phase(key_new,key_ref,iorb,ispin,Nint,phase)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! apply annihilation operator to key_ref
  ! remove electron with spin ispin to orbital with index iorb
  ! output resulting det and phase in key_new and phase
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint,2)
  integer(bit_kind), intent(out) :: key_new(Nint,2)
  double precision, intent(out) :: phase
  
  integer                        :: k,l,i

  double precision, parameter :: p(0:1) = (/ 1.d0, -1.d0 /)
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  key_new=key_ref

  ! alpha det is list of Nint 64-bit ints
  ! k is index of the int where iorb is found
  ! l is index of the bit where iorb is found
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  key_new(k,ispin) = ibclr(key_new(k,ispin),l)
  
  integer(bit_kind) :: parity_filled

  ! I assume here that the ordering is all alpha spinorbs and then all beta spinorbs
  ! If we add an alpha electron, parity is not affected by beta part of determinant
  !   (only need number of alpha occupied orbs below iorb)
 
  ! If we add a beta electron, the parity is affected by alpha part
  !    (need total number of occupied alpha orbs (all of which come before beta)
  !      and total number of beta occupied orbs below iorb)

  if (ispin==1) then
    parity_filled=0_bit_kind
  else
    parity_filled=iand(int(elec_alpha_num,bit_kind),1_bit_kind)
  endif

  ! get parity due to orbs in other ints (with lower indices)
  do i=1,k-1
    parity_filled = iand(int(popcnt(key_ref(i,ispin)),bit_kind),parity_filled)
  enddo
  
  ! get parity due to orbs in same int as iorb
  ! ishft(1_bit_kind,l)-1 has its l rightmost bits set to 1, other bits set to 0
  parity_filled = iand(int(popcnt(iand(ishft(1_bit_kind,l)-1,key_ref(k,ispin))),bit_kind),parity_filled)
  phase = p(iand(1_bit_kind,parity_filled))

end
!BEGIN_PROVIDER [ double precision, mo_mono_elec_integral_diag,(mo_num)]
!  implicit none
!  integer                        :: i
!  BEGIN_DOC
!  ! diagonal elements of mo_mono_elec_integral array
!  END_DOC
!  print*,'Providing the mono electronic integrals (diagonal)'
!
!  do i = 1, mo_num
!    mo_mono_elec_integral_diag(i) = real(mo_mono_elec_integral(i,i))
!  enddo
!
!END_PROVIDER
