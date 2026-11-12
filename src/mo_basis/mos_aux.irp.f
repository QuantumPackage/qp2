
! ---

BEGIN_PROVIDER [double precision, mo_coef_aux, (ao_num,mo_num)]

  implicit none
  BEGIN_DOC
  ! Auxiliary MO coefficients, same (ao_num,mo_num) layout as
  ! mo_coef: coefficient of the i-th AO on the j-th aux MO.
  ! Read from the EZFIO database if the entry exists (master
  ! reads, then MPI broadcast), otherwise set to the $S^{-1/2}$
  ! canonical orthonormalized AOs (ao_ortho_canonical_coef).
  ! Auxiliary basis, probably to go beyond the basis set for
  ! integral evaluation (author's reading, 2026-10-31); introduced
  ! by AbdAmmar (commit 868d5a162, 2023-05-06, "jast 201 added").
  END_DOC
  integer                        :: i, j
  logical                        :: exists
  double precision, allocatable  :: buffer(:,:)

  PROVIDE ezfio_filename

  if (mpi_master) then
    ! Coefs
    call ezfio_has_mo_basis_mo_coef_aux(exists)
  endif



  if (exists) then
    if (mpi_master) then
      call ezfio_get_mo_basis_mo_coef_aux(mo_coef_aux)
      write(*,*) 'Read  mo_coef_aux'
    endif
  else
    ! Orthonormalized AO basis
    do i = 1, mo_num
      do j = 1, ao_num
        mo_coef_aux(j,i) = ao_ortho_canonical_coef(j,i)
      enddo
    enddo
  endif

END_PROVIDER

