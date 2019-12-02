BEGIN_PROVIDER [ double precision, ao_pseudo_integrals, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Pseudo-potential integrals in the |AO| basis set.
  END_DOC

  if (read_ao_integrals_pseudo) then
    call ezfio_get_ao_one_e_ints_ao_integrals_pseudo(ao_pseudo_integrals)
    print *,  'AO pseudopotential integrals read from disk'
  else

    if (do_pseudo) then
      print *, irp_here, 'Not yet implemented'
      stop -1
    endif
  endif

  if (write_ao_integrals_pseudo) then
    call ezfio_set_ao_one_e_ints_ao_integrals_pseudo(ao_pseudo_integrals)
    print *,  'AO pseudopotential integrals written to disk'
  endif

END_PROVIDER


