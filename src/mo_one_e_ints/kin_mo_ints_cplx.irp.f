BEGIN_PROVIDER [complex*16, mo_kinetic_integrals_complex, (mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  !  Kinetic energy integrals in the MO basis
  END_DOC
  integer :: i,j

  print *, 'Providing MO kinetic integrals'
  if (read_mo_integrals_kinetic) then
    call ezfio_get_mo_one_e_ints_mo_integrals_kinetic_complex(mo_kinetic_integrals_complex)
    print *,  'MO kinetic integrals read from disk'
  else
  print *, 'Providing MO kinetic integrals from AO kinetic integrals'
    call ao_to_mo_complex(                                            &
        ao_kinetic_integrals_complex,                                 &
        size(ao_kinetic_integrals_complex,1),                         &
        mo_kinetic_integrals_complex,                                 &
        size(mo_kinetic_integrals_complex,1)                          &
        )
  endif
  if (write_mo_integrals_kinetic) then
    call ezfio_set_mo_one_e_ints_mo_integrals_kinetic_complex(mo_kinetic_integrals_complex)
    print *,  'MO kinetic integrals written to disk'
  endif

END_PROVIDER

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!

BEGIN_PROVIDER [complex*16, mo_kinetic_integrals_kpts, (mo_num_per_kpt,mo_num_per_kpt,kpt_num)]
  implicit none
  BEGIN_DOC
  !  Kinetic energy integrals in the MO basis
  END_DOC
  integer :: i,j

  print *, 'Providing MO kinetic integrals'
  if (read_mo_integrals_kinetic) then
    call ezfio_get_mo_one_e_ints_mo_integrals_kinetic_kpts(mo_kinetic_integrals_kpts)
    print *,  'MO kinetic integrals read from disk'
  else
  print *, 'Providing MO kinetic integrals from AO kinetic integrals'
    call ao_to_mo_kpts(                                            &
        ao_kinetic_integrals_kpts,                                 &
        size(ao_kinetic_integrals_kpts,1),                         &
        mo_kinetic_integrals_kpts,                                 &
        size(mo_kinetic_integrals_kpts,1)                          &
        )
  endif
  if (write_mo_integrals_kinetic) then
    call ezfio_set_mo_one_e_ints_mo_integrals_kinetic_kpts(mo_kinetic_integrals_kpts)
    print *,  'MO kinetic integrals written to disk'
  endif

END_PROVIDER

