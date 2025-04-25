
! ---

BEGIN_PROVIDER [ double precision, ao_integrals_n_e, (ao_num,ao_num)]

  BEGIN_DOC
  !  Nucleus-electron interaction, in the |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
  !
  !  These integrals also contain the pseudopotential integrals.
  END_DOC

  implicit none
  integer          :: num_A, num_B, power_A(3), power_B(3)
  integer          :: i, j, k, l, n_pt_in, m
  double precision :: alpha, beta
  double precision :: A_center(3),B_center(3),C_center(3)
  double precision :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult

  ao_integrals_n_e = 0.d0

  if (read_ao_integrals_n_e) then

    call ezfio_get_ao_one_e_ints_ao_integrals_n_e(ao_integrals_n_e)
    print *,  'AO N-e integrals read from disk'

  else

   call ao_cart_to_ao_basis(ao_cart_integrals_n_e, ao_cart_num, ao_integrals_n_e, ao_num)

    IF(do_pseudo) THEN
       ao_integrals_n_e += ao_pseudo_integrals
    ENDIF
    IF(point_charges) THEN
       ao_integrals_n_e += ao_integrals_pt_chrg
    ENDIF

  endif


  if (write_ao_integrals_n_e) then
    call ezfio_set_ao_one_e_ints_ao_integrals_n_e(ao_integrals_n_e)
    print *,  'AO N-e integrals written to disk'
  endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_integrals_n_e_imag, (ao_num,ao_num)]
  BEGIN_DOC
  !  Nucleus-electron interaction, in the |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
  END_DOC
  implicit none
  double precision               :: alpha, beta
  integer                        :: num_A,num_B
  double precision               :: A_center(3),B_center(3),C_center(3)
  integer                        :: power_A(3),power_B(3)
  integer                        :: i,j,k,l,n_pt_in,m
  double precision               :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult

  if (read_ao_integrals_n_e) then
    call ezfio_get_ao_one_e_ints_ao_integrals_n_e_imag(ao_integrals_n_e_imag)
    print *,  'AO N-e integrals read from disk'
  else
   print *,  irp_here, ': Not yet implemented'
  endif
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_integrals_n_e_per_atom, (ao_num,ao_num,nucl_num)]
  BEGIN_DOC
! Nucleus-electron interaction in the |AO| basis set, per atom A.
!
! :math:`\langle \chi_i | -\frac{1}{|r-R_A|} | \chi_j \rangle`
  END_DOC
  implicit none
  integer :: i
  do i = 1, nucl_num
   call ao_cart_to_ao_basis(ao_cart_integrals_n_e_per_atom(1,1,i), ao_cart_num,ao_integrals_n_e_per_atom(1,1,i), ao_num)
  enddo
END_PROVIDER

