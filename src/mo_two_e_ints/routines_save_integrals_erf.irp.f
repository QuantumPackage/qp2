subroutine save_erf_two_e_integrals_mo
 implicit none
 integer :: i,j,k,l
 PROVIDE mo_two_e_integrals_erf_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_erf',mo_integrals_erf_map)
 call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals_erf('Read')
end


subroutine save_erf_two_e_ints_mo_into_ints_mo
 implicit none
 integer :: i,j,k,l
 PROVIDE mo_two_e_integrals_erf_in_map
 call ezfio_set_work_empty(.False.)
 call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_erf_map)
 call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')
end

