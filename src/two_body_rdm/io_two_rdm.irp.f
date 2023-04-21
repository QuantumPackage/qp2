subroutine write_array_two_rdm(n_orb,nstates,array_tmp,name_file)
 implicit none
 integer, intent(in) :: n_orb,nstates
 character*(128),  intent(in) :: name_file
 double precision, intent(in) :: array_tmp(n_orb,n_orb,n_orb,n_orb,nstates)

 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 PROVIDE ezfio_filename
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'W')
 call lock_io()
 write(i_unit_output)array_tmp
 call unlock_io()
 close(unit=i_unit_output)
end

subroutine read_array_two_rdm(n_orb,nstates,array_tmp,name_file)
 implicit none
 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 integer, intent(in) :: n_orb,nstates
 character*(128),  intent(in)  :: name_file
 double precision, intent(out) :: array_tmp(n_orb,n_orb,n_orb,n_orb,N_states)
 PROVIDE ezfio_filename
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'R')
 call lock_io()
 read(i_unit_output)array_tmp
 call unlock_io()
 close(unit=i_unit_output)
end

