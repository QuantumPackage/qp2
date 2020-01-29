 BEGIN_PROVIDER [double precision, mo_integrals_n_e_real, (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_integrals_n_e_imag, (mo_num,mo_num)]
&BEGIN_PROVIDER [complex*16, mo_integrals_n_e_complex, (mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  !  Kinetic energy integrals in the MO basis
  END_DOC
  integer :: i,j

  if (read_mo_integrals_e_n) then
    mo_integrals_n_e_real = 0.d0
    mo_integrals_n_e_imag = 0.d0
    call ezfio_get_mo_one_e_ints_mo_integrals_e_n_real(mo_integrals_n_e_real)
    call ezfio_get_mo_one_e_ints_mo_integrals_e_n_imag(mo_integrals_n_e_imag)
    print *,  'MO N-e integrals read from disk'
    do i=1,mo_num
      do j=1,mo_num
        mo_integrals_n_e_complex(j,i) = dcmplx(mo_integrals_n_e_real(j,i), &
                                                   mo_integrals_n_e_imag(j,i))
      enddo
    enddo
  else
    call ao_to_mo_complex(                                            &
        ao_integrals_n_e_complex,                                 &
        size(ao_integrals_n_e_complex,1),                         &
        mo_integrals_n_e_complex,                                 &
        size(mo_integrals_n_e_complex,1)                          &
        )
  endif
  if (write_mo_integrals_e_n) then
    !mo_integrals_n_e_real = 0.d0
    !mo_integrals_n_e_imag = 0.d0
    do i=1,mo_num
      do j=1,mo_num
        mo_integrals_n_e_real(j,i)=dble(mo_integrals_n_e_complex(j,i))
        mo_integrals_n_e_imag(j,i)=dimag(mo_integrals_n_e_complex(j,i))
      enddo
    enddo
    call ezfio_set_mo_one_e_ints_mo_integrals_e_n_real(mo_integrals_n_e_real)
    call ezfio_set_mo_one_e_ints_mo_integrals_e_n_imag(mo_integrals_n_e_imag)
    print *,  'MO N-e integrals written to disk'
  endif

END_PROVIDER


