 BEGIN_PROVIDER [double precision, mo_integrals_n_e, (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_integrals_n_e_imag, (mo_num,mo_num)]
&BEGIN_PROVIDER [complex*16, mo_integrals_n_e_complex, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
! Nucleus-electron interaction on the |MO| basis
 END_DOC

  if (is_periodic) then
    integer :: i,j
    if (read_mo_integrals_e_n) then
      call ezfio_get_mo_one_e_ints_mo_integrals_e_n(mo_integrals_n_e)
      call ezfio_get_mo_one_e_ints_mo_integrals_e_n_imag(mo_integrals_n_e_imag)
      print *,  'MO N-e integrals read from disk'
      do i=1,mo_num
        do j=1,mo_num
          mo_integrals_n_e_complex(j,i) = dcmplx(mo_integrals_n_e(j,i), &
                                                 mo_integrals_n_e_imag(j,i))
        enddo
      enddo
    else
      call ao_to_mo_complex(                                            &
               ao_integrals_n_e_complex,                            &
          size(ao_integrals_n_e_complex,1),                         &
               mo_integrals_n_e_complex,                            &
          size(mo_integrals_n_e_complex,1)                          &
          )
      do i=1,mo_num
        do j=1,mo_num
          mo_integrals_n_e(j,i)=dble(mo_integrals_n_e_complex(j,i))
          mo_integrals_n_e_imag(j,i)=dimag(mo_integrals_n_e_complex(j,i))
        enddo
      enddo
    endif
    if (write_mo_integrals_e_n) then
      call ezfio_set_mo_one_e_ints_mo_integrals_e_n(mo_integrals_n_e)
      call ezfio_set_mo_one_e_ints_mo_integrals_e_n_imag(mo_integrals_n_e_imag)
      print *,  'MO N-e integrals written to disk'
    endif
  else
    if (read_mo_integrals_e_n) then
       call ezfio_get_mo_one_e_ints_mo_integrals_e_n(mo_integrals_n_e)
      print *,  'MO N-e integrals read from disk'
    else
      call ao_to_mo(                                                   &
          ao_integrals_n_e,                                       &
          size(ao_integrals_n_e,1),                               &
          mo_integrals_n_e,                                       &
          size(mo_integrals_n_e,1)                                &
          )
    endif
    if (write_mo_integrals_e_n) then
       call ezfio_set_mo_one_e_ints_mo_integrals_e_n(mo_integrals_n_e)
      print *,  'MO N-e integrals written to disk'
    endif
  endif

END_PROVIDER


BEGIN_PROVIDER [double precision, mo_integrals_n_e_per_atom, (mo_num,mo_num,nucl_num)]
 implicit none
 BEGIN_DOC
! mo_integrals_n_e_per_atom(i,j,k) =
! $\langle \phi_i| -\frac{1}{|r-R_k|} | \phi_j \rangle$.
! where R_k is the coordinate of the k-th nucleus.
 END_DOC

 integer :: k
 mo_integrals_n_e_per_atom = 0.d0
 do k = 1, nucl_num
   call ao_to_mo(                                                    &
       ao_integrals_n_e_per_atom(1,1,k),                        &
       size(ao_integrals_n_e_per_atom,1),                       &
       mo_integrals_n_e_per_atom(1,1,k),                        &
       size(mo_integrals_n_e_per_atom,1)                        &
       )
 enddo

END_PROVIDER

