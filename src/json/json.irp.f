BEGIN_PROVIDER [ character*(2048), json_filename ]
 implicit none
 BEGIN_DOC
 ! Path of the JSON result file of this run: <EZFIO>/json/NNNNN.json
 !
 ! The EZFIO setter below creates the json/ directory and marks the EZFIO
 ! file as carrying JSON output; it has to run before the file is opened.
 END_DOC
 integer                      :: counter
 character*(2048)             :: prefix
 logical                      :: exists

 prefix = trim(ezfio_filename)//'/json/'

 call lock_io
 call ezfio_set_json_empty(.False.)
 exists = .True.
 counter = 0
 do while (exists)
   counter += 1
   write(json_filename, '(A,I5.5,A)') trim(prefix), counter, '.json'
   INQUIRE(FILE=trim(json_filename), EXIST=exists)
 enddo
 call unlock_io

END_PROVIDER

subroutine json_init(kind)
  use iso_c_binding, only: c_int
  use qp_json, only: qp_json_begin, json_add
  implicit none
  BEGIN_DOC
  ! Start the JSON logging of this run and write the first version of the
  ! file. Called once per program, before any json_add.
  END_DOC
  character(len=*), intent(in) :: kind

  PROVIDE json_filename
  call qp_json_begin(trim(kind), trim(json_filename), trim(ezfio_filename), &
                     merge(1_c_int, 0_c_int, mpi_master))
  if (mpi_master) then
    ! ezfio_files is excluded from the automatic input capture: a NEED
    ! entry on json would close a dependency cycle. Its one input that
    ! the calculation really inspects (ao_basis checks it against a
    ! format version) is logged here by hand.
    call json_add('/input/ezfio_files/ezfio_convention', ezfio_convention)
  endif
end
