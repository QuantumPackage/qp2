subroutine save_erf_two_e_integrals_ao
  BEGIN_DOC
  ! Saves the map of the AO range-separated (erf) two-electron
  ! integrals to the EZFIO file, in the work directory
  ! (work/ao_ints_erf), marks the work directory as not empty, and
  ! sets io_ao_two_e_integrals_erf to 'Read', so that the
  ! integrals are read back from the file at the next run.
  END_DOC
 implicit none
 integer :: i,j,k,l
 PROVIDE ao_two_e_integrals_erf_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_erf',ao_integrals_erf_map)
 call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals_erf('Read')
end

subroutine save_erf_two_e_ints_ao_into_ints_ao
  BEGIN_DOC
  ! Saves the map of the AO range-separated (erf) two-electron
  ! integrals to the EZFIO file, but in the file of the regular AO
  ! integrals (work/ao_ints), marks the work directory as not
  ! empty, and sets io_ao_two_e_integrals and do_ao_cholesky to
  ! 'Read', so that at the next run the erf integrals are read as
  ! the regular two-electron integrals.
  END_DOC
 implicit none
 integer :: i,j,k,l
 PROVIDE ao_two_e_integrals_erf_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_erf_map)
 call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read')
 call ezfio_set_ao_two_e_ints_do_ao_cholesky('Read')
end

