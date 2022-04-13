subroutine set_multiple_levels_omp(activate)

  BEGIN_DOC
! If true, activate OpenMP nested parallelism. If false, deactivate.
  END_DOC

  implicit none
  logical, intent(in) :: activate

  if (activate) then
    call omp_set_max_active_levels(3)

    IRP_IF SET_NESTED
      call omp_set_nested(.True.)
    IRP_ENDIF

  else

    call omp_set_max_active_levels(1)

    IRP_IF SET_NESTED
      call omp_set_nested(.False.)
    IRP_ENDIF
  end if

end
