subroutine save_erf_two_e_integrals_mo
  BEGIN_DOC
  ! Saves the map of the MO range-separated (erf) two-electron
  ! integrals to the EZFIO file, in the work directory
  ! (work/mo_ints_erf), marks the work directory as not empty, and
  ! sets io_mo_two_e_integrals_erf to 'Read', so that the
  ! integrals are read back from the file at the next run.
  END_DOC
 implicit none
 integer :: i,j,k,l
 PROVIDE mo_two_e_integrals_erf_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_erf',mo_integrals_erf_map)
 call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals_erf('Read')
end


subroutine save_erf_two_e_ints_mo_into_ints_mo
  BEGIN_DOC
  ! Saves the map of the MO range-separated (erf) two-electron
  ! integrals to the EZFIO file, but in the file of the regular MO
  ! integrals (work/mo_ints), marks the work directory as not empty,
  ! and sets io_mo_two_e_integrals to 'Read', so that at the next run
  ! the erf integrals are read as the regular two-electron integrals.
  END_DOC
 implicit none
 integer :: i,j,k,l
 PROVIDE mo_two_e_integrals_erf_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_erf_map)
 call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
end

