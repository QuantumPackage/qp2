 BEGIN_PROVIDER [ character*(128), qp_stop_filename ]
&BEGIN_PROVIDER [ character*(128), qp_kill_filename ]
&BEGIN_PROVIDER [ integer, qp_stop_variable ]
 implicit none
 BEGIN_DOC
 ! Name of the file to check for qp stop
 END_DOC
 qp_stop_filename = trim(ezfio_filename)//'/work/qpstop'
 qp_kill_filename = trim(ezfio_filename)//'/work/qpkill'
 qp_stop_variable = 0
END_PROVIDER

logical function qp_stop()
  implicit none
  BEGIN_DOC
! Checks if the qp_stop command was invoked for the clean termination of the program
  END_DOC
  integer                        :: iunit
  integer, external              :: getUnitAndOpen

  if (qp_stop_variable == 0) then

    INQUIRE(FILE=trim(qp_kill_filename), EXIST=qp_stop)
    if (qp_stop) then
      qp_stop_variable = 1
      ! qp_stop is true
      return
    endif

    INQUIRE(FILE=trim(qp_stop_filename), EXIST=qp_stop)
    if (qp_stop) then
      iunit = getUnitAndOpen(trim(qp_stop_filename),'r')
      close(iunit, STATUS='DELETE')
    endif

  else 

    qp_stop = .True.

  endif
end



