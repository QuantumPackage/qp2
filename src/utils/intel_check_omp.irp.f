subroutine intel_check_omp()

! Doc : idk

  implicit none

  IRP_IF INTEL2021_CHECK_OMP
    call omp_set_max_active_levels(5)
    print*,'INTEL2021_CHECK_OMP: true'
  IRP_ENDIF
  IRP_IF INTEL2019_CHECK_OMP
    call omp_set_nested(.True.)
    print*,'INTEL2019_CHECK_OMP: true'
  IRP_ENDIF
  IRP_IF GNU_CHECK_OMP
    call omp_set_nested(.True.)
    print*,'GNU_CHECK_OMP: true'
  IRP_ENDIF

end
