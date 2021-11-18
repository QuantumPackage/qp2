subroutine intel_check_omp()

! Doc : idk

  implicit none

  IRP_IF INTEL_CHECK_OMP
    call omp_set_max_active_levels(5)
    print*,'INTEL_CHECK_OMP: true'
  IRP_ELSE
    call omp_set_nested(.True.)
    !call omp_set_nested(.False.)
    print*,'INTEL_CHECK_OMP: false'
  IRP_ENDIF

end
