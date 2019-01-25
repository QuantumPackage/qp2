BEGIN_PROVIDER [ character*32,h0_type ]
 implicit none
 BEGIN_DOC
 ! Type of zeroth-order Hamiltonian
 END_DOC
 if (s2_eig) then
   h0_type = 'SOP'
 else
   h0_type = 'EN'
 endif
END_PROVIDER

