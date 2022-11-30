use trexio

BEGIN_PROVIDER [ logical, use_trexio ]
 implicit none
 BEGIN_DOC
 ! Returns the content of the QP_USE_TREXIO variable
 END_DOC
 character(32) :: buffer

 call getenv('QP_USE_TREXIO', buffer)
 if (trim(buffer) == '0') then
   print *, 'Using EZFIO storage'
   use_trexio = .False.
 else
   print *, 'Using TREXIO storage'
   use_trexio = .True.
 endif

END_PROVIDER

BEGIN_PROVIDER [ character*(1024), trexio_filename ]
  implicit none
  BEGIN_DOC
  ! Name of the TREXIO file.
  END_DOC
  trexio_filename = trim(ezfio_work_dir)//'/trexio.h5'
END_PROVIDER

BEGIN_PROVIDER [ integer(trexio_back_end_t), trexio_backend ]
  implicit none
  BEGIN_DOC
  ! Name of the TREXIO file.
  END_DOC
  trexio_backend = TREXIO_HDF5
END_PROVIDER

BEGIN_PROVIDER [ integer(trexio_t), trexio_file ]
  implicit none
  BEGIN_DOC
  ! Name of the TREXIO file.
  END_DOC
  integer (trexio_exit_code) :: rc

  trexio_file = 0_trexio_t
  if (use_trexio) then
    trexio_file = trexio_open(trexio_filename, 'u', trexio_backend, rc)
    call trexio_assert(rc, TREXIO_SUCCESS)
  endif

END_PROVIDER

