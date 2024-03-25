program my_program
  implicit none
  BEGIN_DOC
! This program is there essentially to show how one can use providers in programs 
  END_DOC
 integer :: i,j
 double precision :: accu
 print*,'Trace on the AO basis '
 print*,trace_ao_one_e_ints
 print*,'Trace on the AO basis after projection on the MO basis'
 print*,trace_ao_one_e_ints_from_mo
 print*,'Trace of MO integrals '
 print*,trace_mo_one_e_ints
 print*,'ao_num = ',ao_num
 print*,'mo_num = ',mo_num
 if(ao_num .ne. mo_num)then
  print*,'The AO basis and MO basis are different ...'
  print*,'Trace on the AO basis should not be the same as Trace of MO integrals'
  print*,'Only the second one must be equal to the trace on the MO integrals'
 else 
  print*,'The AO basis and MO basis are the same !'
  print*,'All traces should coincide '
 endif
end
