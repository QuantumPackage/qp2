BEGIN_PROVIDER [double precision, mo_integrals_n_e, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
! Nucleus-electron interaction on the |MO| basis
 END_DOC

  if (read_mo_integrals_n_e) then
     call ezfio_get_mo_one_e_ints_mo_integrals_n_e(mo_integrals_n_e)
    print *,  'MO N-e integrals read from disk'
  else
    call ao_to_mo(                                                   &
        ao_integrals_n_e,                                       &
        size(ao_integrals_n_e,1),                               &
        mo_integrals_n_e,                                       &
        size(mo_integrals_n_e,1)                                &
        )
  endif
  if (write_mo_integrals_n_e) then
     call ezfio_set_mo_one_e_ints_mo_integrals_n_e(mo_integrals_n_e)
    print *,  'MO N-e integrals written to disk'
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

BEGIN_PROVIDER [ double precision, mo_integrals_n_e_diag,(mo_num)]
  implicit none
  integer                        :: i
  BEGIN_DOC
  ! diagonal elements of mo_integrals_n_e or mo_integrals_n_e_complex
  END_DOC
  
  if (is_complex) then
    integer :: k,i_shft
    PROVIDE mo_integrals_n_e_kpts
    do k=1,kpt_num
      i_shft = (k-1)*mo_num_per_kpt
      do i=1,mo_num_per_kpt
        mo_integrals_n_e_diag(i+i_shft) = dble(mo_integrals_n_e_kpts(i,i,k))
      enddo
    enddo
  else
    PROVIDE mo_integrals_n_e
    do i=1,mo_num
      mo_integrals_n_e_diag(i) = mo_integrals_n_e(i,i)
    enddo
  endif
END_PROVIDER
