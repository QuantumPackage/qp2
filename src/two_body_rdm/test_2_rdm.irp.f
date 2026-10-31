program test_2_rdm
 BEGIN_DOC
 ! Test program for the two-electron density matrix (2-RDM)
 ! providers of the two_body_rdm module: it sets read_wf so that
 ! the wave function is read from the EZFIO file, then calls
 ! routine_active_only, routine_full_mos and
 ! routine_active_only_trans, which print the two-electron
 ! repulsion energy obtained from the 2-RDM providers.
 END_DOC
 implicit none
 read_wf = .True.
 touch read_wf 
 call routine_active_only
 call routine_full_mos
 call routine_active_only_trans
end

