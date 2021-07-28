
subroutine four_idx_novvvv_old
  use map_module
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Retransform MO integrals for next CAS-SCF step
  END_DOC
  integer(bit_kind)              :: mask_ijkl(N_int,4)
  integer(bit_kind)              :: mask_ijk(N_int,3)

  print*,'Using partial transformation'
  print*,'It will not transform all integrals with at least 3 indices within the virtuals'
    integer                        :: i,j,k,l
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  I I I I !!!!!!!!!!!!!!!!!!!!
    ! (core+inact+act) ^ 4
    ! <ii|ii>
    print*, ''
    print*, '<ii|ii>'
    do i = 1,N_int
      mask_ijkl(i,1) =  core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,2) =  core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,3) =  core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,4) =  core_inact_act_bitmask_4(i,1)
    enddo
    call add_integrals_to_map(mask_ijkl)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  I I V V !!!!!!!!!!!!!!!!!!!!
    ! (core+inact+act) ^ 2  (virt) ^2
    ! <iv|iv>  = J_iv
    print*, ''
    print*, '<iv|iv>'
    do i = 1,N_int
      mask_ijkl(i,1) =  core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,2) =  virt_bitmask(i,1)
      mask_ijkl(i,3) =  core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,4) =  virt_bitmask(i,1)
    enddo
    call add_integrals_to_map(mask_ijkl)

    ! (core+inact+act) ^ 2  (virt) ^2
    ! <ii|vv> = (iv|iv)
    print*, ''
    print*, '<ii|vv>'
    do i = 1,N_int
      mask_ijkl(i,1) = core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,2) = core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,3) = virt_bitmask(i,1)
      mask_ijkl(i,4) = virt_bitmask(i,1)
    enddo
    call add_integrals_to_map(mask_ijkl)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! V V V !!!!!!!!!!!!!!!!!!!!!!!
!    if(.not.no_vvv_integrals)then
      print*, ''
      print*, '<rv|sv> and <rv|vs>'
      do i = 1,N_int
        mask_ijk(i,1) =  virt_bitmask(i,1)
        mask_ijk(i,2) =  virt_bitmask(i,1)
        mask_ijk(i,3) =  virt_bitmask(i,1)
      enddo
      call add_integrals_to_map_three_indices(mask_ijk)
!    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  I I I V !!!!!!!!!!!!!!!!!!!!
    ! (core+inact+act) ^ 3  (virt) ^1
    ! <iv|ii>
    print*, ''
    print*, '<iv|ii>'
    do i = 1,N_int
      mask_ijkl(i,1) =  core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,2) =  core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,3) =  core_inact_act_bitmask_4(i,1)
      mask_ijkl(i,4) =  virt_bitmask(i,1)
    enddo
    call add_integrals_to_map(mask_ijkl)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  I V V V !!!!!!!!!!!!!!!!!!!!
    ! (core+inact+act) ^ 1  (virt) ^3
    ! <iv|vv>
!    if(.not.no_ivvv_integrals)then
      print*, ''
      print*, '<iv|vv>'
      do i = 1,N_int
        mask_ijkl(i,1) =  core_inact_act_bitmask_4(i,1)
        mask_ijkl(i,2) =  virt_bitmask(i,1)
        mask_ijkl(i,3) =  virt_bitmask(i,1)
        mask_ijkl(i,4) =  virt_bitmask(i,1)
      enddo
      call add_integrals_to_map_no_exit_34(mask_ijkl)
end
