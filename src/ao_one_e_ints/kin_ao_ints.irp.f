
! ---

 BEGIN_PROVIDER [ double precision, ao_deriv2_x, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_deriv2_y, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_deriv2_z, (ao_num, ao_num) ]

  BEGIN_DOC
  ! Second derivative matrix elements in the |AO| basis.
  !
  ! .. math::
  !
  !   {\tt ao\_deriv2\_x} =
  !   \langle \chi_i(x,y,z) | \frac{\partial^2}{\partial x^2} |\chi_j (x,y,z) \rangle
  !
  END_DOC

  implicit none
  integer          :: i, j

  call ao_cart_to_ao_basis(ao_cart_deriv2_x, ao_cart_num, ao_deriv2_x, ao_num)
  call ao_cart_to_ao_basis(ao_cart_deriv2_y, ao_cart_num, ao_deriv2_y, ao_num)
  call ao_cart_to_ao_basis(ao_cart_deriv2_z, ao_cart_num, ao_deriv2_z, ao_num)

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_kinetic_integrals, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Kinetic energy integrals in the |AO| basis.
  !
  ! $\langle \chi_i |\hat{T}| \chi_j \rangle$
  !
  END_DOC
  integer                        :: i,j,k,l

  if (read_ao_integrals_kinetic) then
    call ezfio_get_ao_one_e_ints_ao_integrals_kinetic(ao_kinetic_integrals)
    print *,  'AO kinetic integrals read from disk'
  else
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP  PRIVATE(i,j) &
    !$OMP  SHARED(ao_num, ao_kinetic_integrals,ao_deriv2_x,ao_deriv2_y,ao_deriv2_z)
    do j = 1, ao_num
      do i = 1, ao_num
      ao_kinetic_integrals(i,j) = -0.5d0 * (ao_deriv2_x(i,j) + ao_deriv2_y(i,j) + ao_deriv2_z(i,j) )
      enddo
    enddo
    !$OMP END PARALLEL DO
  endif
  if (write_ao_integrals_kinetic) then
    call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(ao_kinetic_integrals)
    print *,  'AO kinetic integrals written to disk'
  endif
END_PROVIDER

BEGIN_PROVIDER [double precision, ao_kinetic_integrals_imag, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Kinetic energy integrals in the |AO| basis.
  !
  ! $\langle \chi_i |\hat{T}| \chi_j \rangle$
  !
  END_DOC
  integer                        :: i,j,k,l

  if (read_ao_integrals_kinetic) then
    call ezfio_get_ao_one_e_ints_ao_integrals_kinetic(ao_kinetic_integrals_imag)
    print *,  'AO kinetic integrals read from disk'
  else
    print *,  irp_here, ': Not yet implemented'
  endif
  if (write_ao_integrals_kinetic) then
    call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(ao_kinetic_integrals_imag)
    print *,  'AO kinetic integrals written to disk'
  endif
END_PROVIDER

