BEGIN_PROVIDER [double precision, mo_pseudo_integrals, (mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  ! Pseudopotential integrals in |MO| basis
  END_DOC

  if (read_mo_integrals_pseudo) then
    call ezfio_get_mo_one_e_ints_mo_integrals_pseudo(mo_pseudo_integrals)
    print *,  'MO pseudopotential integrals read from disk'
  else if (do_pseudo) then
      call ao_to_mo(                                                   &
        ao_pseudo_integrals,                                         &
        size(ao_pseudo_integrals,1),                                 &
        mo_pseudo_integrals,                                         &
        size(mo_pseudo_integrals,1)                                  &
        )
  else
      mo_pseudo_integrals = 0.d0
  endif

  if (write_mo_integrals_pseudo) then
    call ezfio_set_mo_one_e_ints_mo_integrals_pseudo(mo_pseudo_integrals)
    print *,  'MO pseudopotential integrals written to disk'
  endif

END_PROVIDER


