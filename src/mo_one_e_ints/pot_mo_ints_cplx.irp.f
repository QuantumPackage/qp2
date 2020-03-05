BEGIN_PROVIDER [complex*16, mo_integrals_n_e_complex, (mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  !  Kinetic energy integrals in the MO basis
  END_DOC
  integer :: i,j

  if (read_mo_integrals_e_n) then
    call ezfio_get_mo_one_e_ints_mo_integrals_e_n_complex(mo_integrals_n_e_complex)
    print *,  'MO N-e integrals read from disk'
  else
    call ao_to_mo_complex(                                            &
        ao_integrals_n_e_complex,                                 &
        size(ao_integrals_n_e_complex,1),                         &
        mo_integrals_n_e_complex,                                 &
        size(mo_integrals_n_e_complex,1)                          &
        )
  endif
  if (write_mo_integrals_e_n) then
    call ezfio_set_mo_one_e_ints_mo_integrals_e_n_complex(mo_integrals_n_e_complex)
    print *,  'MO N-e integrals written to disk'
  endif

END_PROVIDER


