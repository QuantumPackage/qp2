program kick_the_mos

  !BEGIN_DOC
  ! To do a small rotation of the MOs
  !END_DOC

  implicit none

  kick_in_mos = .True.
  TOUCH kick_in_mos

  print*, 'Security mo_class:', security_mo_class

  ! The default mo_classes are setted only if the MOs to localize are not specified
  if (security_mo_class .and. (dim_list_act_orb == mo_num .or. &
      dim_list_core_orb + dim_list_act_orb == mo_num)) then

    print*, 'WARNING'
    print*, 'You must set different mo_class with qp set_mo_class'
    print*, 'If you want to kick all the orbital:'
    print*, 'qp set Orbital_optimization security_mo_class false'
    print*, ''
    print*, 'abort'

    call abort
  
  endif
  
  call apply_pre_rotation
  
end
