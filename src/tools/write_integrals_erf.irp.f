program write_integrals_erf
 implicit none
 BEGIN_DOC
 ! Saves the two-electron integrals with the $erf(\mu r_{12})/r_{12}$
 ! oprerator into the EZFIO directory.
 END_DOC
 io_mo_two_e_integrals = 'None'
 touch io_mo_two_e_integrals
 io_ao_two_e_integrals = 'None'
 touch io_ao_two_e_integrals
 call routine

end

subroutine routine
 implicit none
 BEGIN_DOC
 ! Saves the maps of the AO and MO range-separated (erf)
 ! two-electron integrals to the |EZFIO| work directory
 ! (work/ao_ints_erf and work/mo_ints_erf), so that they are read
 ! back at the next run.
 END_DOC
 call save_erf_two_e_integrals_ao
 call save_erf_two_e_integrals_mo

end

