program basis_corr_su_pbe_ot
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf
  no_core_density = .True.
  touch no_core_density
  if(io_mo_two_e_integrals .ne. "Read")then
   provide ao_two_e_integrals_in_map
  endif
  provide mo_two_e_integrals_in_map
  call print_su_pbe_ot
 
end

subroutine print_su_pbe_ot
 implicit none
 integer :: istate 
 do istate = 1, N_states
  write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD PBE-UEG       , state ',istate,' = ',ecmd_pbe_ueg_mu_of_r(istate)
  write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ecmd_pbe_ueg_test  , state ',istate,' = ',ecmd_pbe_ueg_test(istate)
 enddo
! do istate = 1, N_states
!  write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD SU-PBE-OT     , state ',istate,' = ',ecmd_pbe_on_top_su_mu_of_r(istate)
! enddo

end
