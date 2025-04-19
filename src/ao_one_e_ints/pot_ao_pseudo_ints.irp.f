BEGIN_PROVIDER [ double precision, ao_pseudo_integrals, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Pseudo-potential integrals in the |AO| basis set.
  END_DOC

  if (read_ao_integrals_pseudo) then
    call ezfio_get_ao_one_e_ints_ao_integrals_pseudo(ao_pseudo_integrals)
    print *,  'AO pseudopotential integrals read from disk'
  else
   call ao_cart_to_ao_basis(ao_cart_pseudo_integrals, ao_cart_num, ao_pseudo_integrals, ao_num)
  endif

  if (write_ao_integrals_pseudo) then
    call ezfio_set_ao_one_e_ints_ao_integrals_pseudo(ao_pseudo_integrals)
    print *,  'AO pseudopotential integrals written to disk'
  endif

END_PROVIDER

