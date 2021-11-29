program check_omp_v2
 
  use omp_lib

  implicit none

  integer :: accu, accu2
  integer :: s, n_setting
  logical :: verbose, test_versions
  logical, allocatable :: is_working(:)

  verbose = .False.
  test_versions = .True.
  n_setting = 4

  allocate(is_working(n_setting))

  is_working = .False.

  ! set the number of threads
  call omp_set_num_threads(2) 

  do s = 1, n_setting

    accu = 0
    accu2 = 0

    call omp_set_max_active_levels(1)
    call omp_set_nested(.False.)

    if (s==1) then
      !call set_multiple_levels_omp()
      cycle
    elseif (s==2) then
      call omp_set_max_active_levels(5)
    elseif (s==3) then
      call omp_set_nested(.True.)
    else
      call omp_set_nested(.True.)
      call omp_set_max_active_levels(5)
    endif

    ! Level 1
    !$OMP PARALLEL 
      if (verbose) then
        print*,'Num threads level 1:',omp_get_num_threads()
      endif

      ! Level 2
      !$OMP PARALLEL 
        if (verbose) then
          print*,'Num threads level 2:',omp_get_num_threads()
        endif

        ! Level 3
        !$OMP PARALLEL 
          if (verbose) then
            print*,'Num threads level 3:',omp_get_num_threads()
          endif

          call check_omp_in_subroutine(accu2)
     
          ! Level 4
          !$OMP PARALLEL 

            if (verbose) then
              print*,'Num threads level 4:',omp_get_num_threads()
            endif

           !$OMP ATOMIC    
           accu = accu + 1
           !$OMP END ATOMIC 
       
          !$OMP END PARALLEL
     
     
        !$OMP END PARALLEL


      !$OMP END PARALLEL


    !$OMP END PARALLEL

  if (verbose) then
    print*,'Setting:',s,'accu=',accu
    print*,'Setting:',s,'accu2=',accu2
  endif

  if (accu == 16 .and. accu2 == 16) then
    is_working(s) = .True.
  endif

  enddo

  if (verbose) then
    if (is_working(2)) then
      print*,'The parallelization works on 4 levels with:'
      print*,'call omp_set_max_active_levels(5)'
      print*,''
      print*,'Please use the irpf90 flags -DSET_MAX_ACT in qp2/config/${compiler_name}.cfg'
    elseif (is_working(3)) then
      print*,'The parallelization works on 4 levels with:'
      print*,'call omp_set_nested(.True.)'
      print*,''
      print*,'Please use the irpf90 flag -DSET_NESTED in qp2/config/${compiler_name}.cfg'
    elseif (is_working(4)) then
      print*,'The parallelization works on 4 levels with:'
      print*,'call omp_set_nested(.True.)'
      print*,'+'
      print*,'call omp_set_max_active_levels(5)'
      print*,''
      print*,'Please use the irpf90 flags -DSET_NESTED -DSET_MAX_ACT in qp2/config/${compiler_name}.cfg'
    else
      print*,'The parallelization on multiple levels does not work with:'
      print*,'call omp_set_max_active_levels(5)'
      print*,'or'
      print*,'call omp_set_nested(.True.)'
      print*,'or'
      print*,'call omp_set_nested(.True.)'
      print*,'+'
      print*,'call omp_set_max_active_levels(5)'
      print*,''
      print*,'Try an other compiler and good luck...'
    endif

   ! if (is_working(1)) then
   !   print*,''
   !   print*,'=========================================================='
   !   print*,'Your actual set up works for parallelization with 4 levels'
   !   print*,'=========================================================='
   !   print*,''
   ! else
   !   print*,''
   !   print*,'==================================================================='
   !   print*,'Your actual set up does not work for parallelization with 4 levels'
   !   print*,'Please look at the previous messages to understand the requirements'
   !   print*,'==================================================================='
   !   print*,''
   ! endif
  endif

  ! List of working flags
  if (test_versions) then
    print*,'Tests:',is_working(2:4)
  endif

  ! IRPF90_FLAGS
  if (is_working(2)) then
    print*,'-DSET_MAX_ACT'
  elseif (is_working(3)) then
    print*,'-DSET_NESTED'
  elseif (is_working(4)) then
    print*,'-DSET_MAX_ACT -DSET_NESTED'
  else
    print*,'ERROR'
  endif

end

subroutine check_omp_in_subroutine(accu2)

  implicit none
 
  integer, intent(inout) :: accu2

  !$OMP PARALLEL

  !$OMP ATOMIC
  accu2 = accu2 + 1
  !$OMP END ATOMIC

  !$OMP END PARALLEL

end
