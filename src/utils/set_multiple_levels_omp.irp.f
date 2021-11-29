subroutine set_multiple_levels_omp()

! Doc : idk

  implicit none

  IRP_IF SET_MAX_ACT
    !print*,'SET_MAX_ACT: True, call omp_set_max_active_levels(5)'
    call omp_set_max_active_levels(5)
  IRP_ENDIF
  IRP_IF SET_NESTED
    !print*,'SET_NESTED: True, call omp_set_nested(.True.)'
    call omp_set_nested(.True.)
  IRP_ENDIF

end
