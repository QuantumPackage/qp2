 BEGIN_PROVIDER [double precision, mo_pseudo_integrals_real, (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_pseudo_integrals_imag, (mo_num,mo_num)]
&BEGIN_PROVIDER [complex*16, mo_pseudo_integrals_complex, (mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  ! Pseudopotential integrals in |MO| basis
  END_DOC
  integer :: i,j

  if (read_mo_integrals_pseudo) then
    mo_pseudo_integrals_real = 0.d0
    mo_pseudo_integrals_imag = 0.d0
    call ezfio_get_mo_one_e_ints_mo_integrals_pseudo_real(mo_pseudo_integrals_real)
    call ezfio_get_mo_one_e_ints_mo_integrals_pseudo_imag(mo_pseudo_integrals_imag)
    print *,  'MO pseudopotential integrals read from disk'
    do i=1,mo_num
      do j=1,mo_num
        mo_pseudo_integrals_complex(j,i) = dcmplx(mo_pseudo_integrals_real(j,i), &
                                                  mo_pseudo_integrals_imag(j,i))
      enddo
    enddo
  else if (do_pseudo) then
    call ao_to_mo_complex(                                            &
        ao_pseudo_integrals_complex,                                 &
        size(ao_pseudo_integrals_complex,1),                         &
        mo_pseudo_integrals_complex,                                 &
        size(mo_pseudo_integrals_complex,1)                          &
        )
    do i=1,mo_num
      do j=1,mo_num
        mo_pseudo_integrals_real(j,i)=dble(mo_pseudo_integrals_complex(j,i))
        mo_pseudo_integrals_imag(j,i)=dimag(mo_pseudo_integrals_complex(j,i))
      enddo
    enddo
  else
    mo_pseudo_integrals_real = 0.d0
    mo_pseudo_integrals_imag = 0.d0
    mo_pseudo_integrals_complex = (0.d0,0.d0)
  endif
  if (write_mo_integrals_pseudo) then
    !mo_pseudo_integrals_real = 0.d0
    !mo_pseudo_integrals_imag = 0.d0
    do i=1,mo_num
      do j=1,mo_num
        mo_pseudo_integrals_real(j,i)=dble(mo_pseudo_integrals_complex(j,i))
        mo_pseudo_integrals_imag(j,i)=dimag(mo_pseudo_integrals_complex(j,i))
      enddo
    enddo
    call ezfio_set_mo_one_e_ints_mo_integrals_pseudo_real(mo_pseudo_integrals_real)
    call ezfio_set_mo_one_e_ints_mo_integrals_pseudo_imag(mo_pseudo_integrals_imag)
    print *,  'MO pseudopotential integrals written to disk'
  endif

END_PROVIDER


