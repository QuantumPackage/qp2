subroutine write_array_two_rdm(n_orb,nstates,array_tmp,name_file)
 BEGIN_DOC
 ! Writes array_tmp as an unformatted file named <name_file> in
 ! the EZFIO work directory, guarded by lock_io()/unlock_io().
 ! array_tmp(i,j,k,l,istate) is the 2-RDM of state istate over
 ! n_orb orbitals, in the MO basis, physicists' notation.
 ! The callers in act_2_rdm.irp.f store the array before the
 ! factor of 2 that they apply to the provided value.
 END_DOC
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
 BEGIN_DOC
 ! Reads array_tmp from the unformatted file named <name_file> in
 ! the EZFIO work directory, guarded by lock_io()/unlock_io().
 ! array_tmp(i,j,k,l,istate) is the 2-RDM of state istate over
 ! n_orb orbitals, in the MO basis, physicists' notation.
 ! The callers in act_2_rdm.irp.f read the array before the
 ! factor of 2 that they apply to the provided value.
 END_DOC
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


subroutine write_array_two_trans_rdm(n_orb,nstates,array_tmp,name_file)
 BEGIN_DOC
 ! Writes array_tmp as an unformatted file named <name_file> in
 ! the EZFIO work directory, guarded by lock_io()/unlock_io().
 ! array_tmp(i,j,k,l,istate,jstate) is the transition 2-RDM
 ! between states istate and jstate over n_orb orbitals, in the
 ! MO basis, physicists' notation, as stored by
 ! act_2_rdm_trans_spin_trace_mo.
 END_DOC
 implicit none
 integer, intent(in) :: n_orb,nstates
 character*(128),  intent(in) :: name_file
 double precision, intent(in) :: array_tmp(n_orb,n_orb,n_orb,n_orb,nstates,nstates)

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

subroutine read_array_two_trans_rdm(n_orb,nstates,array_tmp,name_file)
 BEGIN_DOC
 ! Reads array_tmp from the unformatted file named <name_file> in
 ! the EZFIO work directory, guarded by lock_io()/unlock_io().
 ! array_tmp(i,j,k,l,istate,jstate) is the transition 2-RDM
 ! between states istate and jstate over n_orb orbitals, in the
 ! MO basis, physicists' notation, as stored by
 ! act_2_rdm_trans_spin_trace_mo.
 END_DOC
 implicit none
 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 integer, intent(in) :: n_orb,nstates
 character*(128),  intent(in)  :: name_file
 double precision, intent(out) :: array_tmp(n_orb,n_orb,n_orb,n_orb,N_states,nstates)
 PROVIDE ezfio_filename
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'R')
 call lock_io()
 read(i_unit_output)array_tmp
 call unlock_io()
 close(unit=i_unit_output)
end

