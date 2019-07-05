program print_2rdm
 implicit none
 BEGIN_DOC
 ! get the active part of the bielectronic energy on a given wave function.
 !
 ! useful to test the active part of the spin trace 2 rdms
 END_DOC
 no_vvvv_integrals = .True.
 read_wf = .True.
 touch read_wf no_vvvv_integrals
 call routine
end

subroutine routine
 print *,  psi_energy_with_nucl_rep
end
