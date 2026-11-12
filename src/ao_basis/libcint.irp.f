
BEGIN_PROVIDER [logical, use_libcint]

  implicit none

  BEGIN_DOC
  ! If true, compute the AO integrals with libcint.
  !
  ! The request is only a preference: it is ignored when the basis or the
  ! Hamiltonian cannot be represented in libcint, see
  ! :c:data:`use_libcint_effective`.
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  use_libcint = .True.
  if (mpi_master) then
    call ezfio_has_ao_basis_use_libcint(has)
    if (has) then
      call ezfio_get_ao_basis_use_libcint(use_libcint)
    else
      call ezfio_set_ao_basis_use_libcint(use_libcint)
    endif
  endif

END_PROVIDER
