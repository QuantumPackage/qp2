BEGIN_PROVIDER [complex*16, mo_pseudo_integrals_complex, (mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  ! Pseudopotential integrals in |MO| basis
  END_DOC
  integer :: i,j

  if (read_mo_integrals_pseudo) then
    call ezfio_get_mo_one_e_ints_mo_integrals_pseudo_complex(mo_pseudo_integrals_complex)
    print *,  'MO pseudopotential integrals read from disk'
  else if (do_pseudo) then
    call ao_to_mo_complex(                                            &
        ao_pseudo_integrals_complex,                                 &
        size(ao_pseudo_integrals_complex,1),                         &
        mo_pseudo_integrals_complex,                                 &
        size(mo_pseudo_integrals_complex,1)                          &
        )
  else
    mo_pseudo_integrals_complex = (0.d0,0.d0)
  endif
  if (write_mo_integrals_pseudo) then
    call ezfio_set_mo_one_e_ints_mo_integrals_pseudo_complex(mo_pseudo_integrals_complex)
    print *,  'MO pseudopotential integrals written to disk'
  endif

END_PROVIDER


