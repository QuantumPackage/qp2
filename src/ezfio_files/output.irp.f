 BEGIN_PROVIDER [ double precision, output_wall_time_0 ]
&BEGIN_PROVIDER [ double precision, output_cpu_time_0 ]
  implicit none
  BEGIN_DOC
  ! Initial CPU and wall times when printing in the output files
  END_DOC
  call cpu_time(output_wall_time_0)
  call wall_time(output_wall_time_0)
END_PROVIDER


subroutine write_time(iunit)
  implicit none
  BEGIN_DOC
  ! Write a time stamp in the output for chronological reconstruction
  END_DOC
  integer, intent(in)            :: iunit
  double precision               :: wt, ct
  if (.not.mpi_master) then
    return
  endif
  write(6,*)
  call print_memory_usage()
  call cpu_time(ct)
  ct = ct - output_cpu_time_0
  call wall_time(wt)
  wt = wt - output_wall_time_0
  write(6,'(A,F14.6,A,F14.6,A)') &
    '.. >>>>> [ WALL TIME: ', wt, '  s ] [ CPU  TIME: ', ct, '  s ] <<<<< ..'
  write(6,*)
end

subroutine write_double(iunit,value,label)
  implicit none
  BEGIN_DOC
  ! Write a double precision value in output
  END_DOC
  if (.not.mpi_master) then
    return
  endif
  integer, intent(in)            :: iunit
  double precision               :: value
  character*(*)                  :: label
  character*(64), parameter      :: f = '(A50,G24.16)'
  character*(50)                 :: newlabel
  write(newlabel,'(A,A)') '* ',trim(label)
  write(6,f) newlabel, value
end


subroutine write_int(iunit,value,label)
  implicit none
  BEGIN_DOC
  ! Write an integer value in output
  END_DOC
  if (.not.mpi_master) then
    return
  endif
  integer, intent(in)            :: iunit
  integer                        :: value
  character*(*)                  :: label
  character*(64), parameter      :: f = '(A50,I16)'
  character*(50)                 :: newlabel
  write(newlabel,'(A,A)') '* ',trim(label)
  write(6,f) newlabel, value
end


subroutine write_bool(iunit,value,label)
  implicit none
  BEGIN_DOC
  ! Write an logical value in output
  END_DOC
  if (.not.mpi_master) then
    return
  endif
  integer, intent(in)            :: iunit
  logical                        :: value
  character*(*)                  :: label
  character*(64), parameter      :: f = '(A50,L1)'
  character*(50)                 :: newlabel
  write(newlabel,'(A,A)') '* ',trim(label)
  write(6,f) newlabel, value
end


