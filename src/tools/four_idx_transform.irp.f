program four_idx_transform
  implicit none
  BEGIN_DOC
! 4-index transformation of two-electron integrals from |AO| to |MO| integrals.
!
! This program will compute the two-electron integrals on the |MO| basis and store it into the |EZFIO| folder.
!
! This program can be useful if the AO --> MO transformation is an expensive step by itself.
  END_DOC

  io_mo_two_e_integrals = 'Write'
  SOFT_TOUCH io_mo_two_e_integrals
  if (.true.) then
    PROVIDE mo_two_e_integrals_in_map
  endif
end
