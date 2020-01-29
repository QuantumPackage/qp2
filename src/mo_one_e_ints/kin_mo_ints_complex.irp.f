 BEGIN_PROVIDER [double precision, mo_kinetic_integrals_real, (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_kinetic_integrals_imag, (mo_num,mo_num)]
&BEGIN_PROVIDER [complex*16, mo_kinetic_integrals_complex, (mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  !  Kinetic energy integrals in the MO basis
  END_DOC
  integer :: i,j

  if (read_mo_integrals_kinetic) then
    mo_kinetic_integrals_real = 0.d0
    mo_kinetic_integrals_imag = 0.d0
    call ezfio_get_mo_one_e_ints_mo_integrals_kinetic_real(mo_kinetic_integrals_real)
    call ezfio_get_mo_one_e_ints_mo_integrals_kinetic_imag(mo_kinetic_integrals_imag)
    print *,  'MO kinetic integrals read from disk'
    do i=1,mo_num
      do j=1,mo_num
        mo_kinetic_integrals_complex(j,i) = dcmplx(mo_kinetic_integrals_real(j,i), &
                                                   mo_kinetic_integrals_imag(j,i))
      enddo
    enddo
  else
    call ao_to_mo_complex(                                            &
        ao_kinetic_integrals_complex,                                 &
        size(ao_kinetic_integrals_complex,1),                         &
        mo_kinetic_integrals_complex,                                 &
        size(mo_kinetic_integrals_complex,1)                          &
        )
  endif
  if (write_mo_integrals_kinetic) then
    !mo_kinetic_integrals_real = 0.d0
    !mo_kinetic_integrals_imag = 0.d0
    do i=1,mo_num
      do j=1,mo_num
        mo_kinetic_integrals_real(j,i)=dble(mo_kinetic_integrals_complex(j,i))
        mo_kinetic_integrals_imag(j,i)=dimag(mo_kinetic_integrals_complex(j,i))
      enddo
    enddo
    call ezfio_set_mo_one_e_ints_mo_integrals_kinetic_real(mo_kinetic_integrals_real)
    call ezfio_set_mo_one_e_ints_mo_integrals_kinetic_imag(mo_kinetic_integrals_imag)
    print *,  'MO kinetic integrals written to disk'
  endif

END_PROVIDER

