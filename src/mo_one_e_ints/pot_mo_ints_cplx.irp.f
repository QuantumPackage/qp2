BEGIN_PROVIDER [complex*16, mo_integrals_n_e_complex, (mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  !  Kinetic energy integrals in the MO basis
  END_DOC
  integer :: i,j
  
  print *, 'Providing MO N-e integrals'
  if (read_mo_integrals_n_e) then
    call ezfio_get_mo_one_e_ints_mo_integrals_n_e_complex(mo_integrals_n_e_complex)
    print *,  'MO N-e integrals read from disk'
  else
  print *, 'Providing MO N-e integrals from AO N-e integrals'
    call ao_to_mo_complex(                                            &
        ao_integrals_n_e_complex,                                 &
        size(ao_integrals_n_e_complex,1),                         &
        mo_integrals_n_e_complex,                                 &
        size(mo_integrals_n_e_complex,1)                          &
        )
  endif
  if (write_mo_integrals_n_e) then
    call ezfio_set_mo_one_e_ints_mo_integrals_n_e_complex(mo_integrals_n_e_complex)
    print *,  'MO N-e integrals written to disk'
  endif

END_PROVIDER

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!

BEGIN_PROVIDER [complex*16, mo_integrals_n_e_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num)]
  implicit none
  BEGIN_DOC
  !  Kinetic energy integrals in the MO basis
  END_DOC
  integer :: i,j
  
  print *, 'Providing MO N-e integrals'
  if (read_mo_integrals_n_e) then
    call ezfio_get_mo_one_e_ints_mo_integrals_n_e_kpts(mo_integrals_n_e_kpts)
    print *,  'MO N-e integrals read from disk'
  else
  print *, 'Providing MO N-e integrals from AO N-e integrals'
    call ao_to_mo_kpts(                                            &
        ao_integrals_n_e_kpts,                                 &
        size(ao_integrals_n_e_kpts,1),                         &
        mo_integrals_n_e_kpts,                                 &
        size(mo_integrals_n_e_kpts,1)                          &
        )
  endif
  if (write_mo_integrals_n_e) then
    call ezfio_set_mo_one_e_ints_mo_integrals_n_e_kpts(mo_integrals_n_e_kpts)
    print *,  'MO N-e integrals written to disk'
  endif

END_PROVIDER
