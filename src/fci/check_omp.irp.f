program check_omp
  
  use omp_lib

  implicit none

  integer :: i,j,k,l,m,n,x,z,setting
  double precision :: w1,w2,c1,c2
  double precision, allocatable :: accu(:,:,:,:) 
  logical :: must_exit, verbose, is_working

  x = 4
  allocate(accu(x,x,x,x))

  verbose = .False.

  accu = 0d0
  must_exit = .False.

  !$OMP PARALLEL
  if (OMP_GET_NUM_THREADS() == 1) then
    print*,''
    print*,'1 thread, no parallelization possible'
    print*,''
    must_exit=.True.
  endif
  !$OMP END PARALLEL
  if (must_exit) then
    call abort
  endif

  ! reset the number of max active levels
  !call omp_set_max_active_levels(1)

  !print*,'omp_get_max_active_levels:',omp_get_max_active_levels()
  !call intel_check_omp()
  !print*,'omp_get_max_active_levels:',omp_get_max_active_levels()

  ! set the number of threads
  call omp_set_num_threads(2) 

  do z = 1, 4

    if (must_exit) then
      exit
    endif 

    call omp_set_max_active_levels(1)
    call omp_set_nested(.False.)

    if (z==1) then
      call test_set_multiple_levels_omp()
      !call test_set_multiple_levels_omp
    elseif (z==2) then
      call omp_set_max_active_levels(5)
    elseif (z==3) then
      call omp_set_nested(.True.)
    else
      call omp_set_nested(.True.) 
      call omp_set_max_active_levels(5)
    endif

    setting = z-1

    !$OMP PARALLEL &
      !$OMP PRIVATE(i,j,k,l,m,n) &
      !$OMP SHARED(accu) 
  
      if (verbose) then
        print*,'Nb threads level 1:', omp_get_num_threads()
      endif
  
      !$OMP MASTER
      if (omp_get_num_threads()==1) then
        print*,'Setting',setting,"error at level 1"
        setting = -1
      endif
      !$OMP END MASTER
      
  !    !$OMP DO
  !    do l = 1, x
  !      do k = 1, x
  !        do j = 1, x
  !          do i = 1, x
  !            accu(i,j,k,l) = accu(i,j,k,l) + 1d0
  !          enddo
  !        enddo
  !      enddo
  !    enddo
  !    !$OMP END DO
  
      !$OMP PARALLEL & 
        !$OMP PRIVATE(i,j,k,l,m,n) &
        !$OMP SHARED(accu)
  
        if (verbose) then
          print*,'Nb threads level 2:', omp_get_num_threads()
        endif
  
        !$OMP MASTER 
        if (omp_get_num_threads()==1 .and. setting >= 0) then
          print*,'Setting',setting,"error at level 2"
          setting = -1
        endif
        !$OMP END MASTER
        
  !      !$OMP DO
  !      do l = 1, x
  !        do k = 1, x
  !          do j = 1, x
  !            do i = 1, x
  !              accu(i,j,k,l) = accu(i,j,k,l)+ 1d0
  !            enddo
  !          enddo
  !        enddo
  !      enddo
  !      !$OMP END DO
  
        !$OMP PARALLEL & 
          !$OMP PRIVATE(i,j,k,l,m,n) &
          !$OMP SHARED(accu)
          
          if (verbose) then 
            print*,'Nb threads level 3:', omp_get_num_threads()
          endif
  
          !$OMP MASTER
          if (omp_get_num_threads()==1 .and. setting >= 0) then
           print*,'Setting',setting,"error at level 3"
           setting = -1
          endif
          !$OMP END MASTER
  
  !        !$OMP DO
  !        do l = 1, x
  !          do k = 1, x
  !            do j = 1, x
  !              do i = 1, x
  !                accu(i,j,k,l) = accu(i,j,k,l)+ 1d0
  !              enddo
  !            enddo
  !          enddo
  !        enddo
  !        !$OMP END DO
  
          !$OMP PARALLEL &
            !$OMP PRIVATE(i,j,k,l,m,n) &
            !$OMP SHARED(accu)
            
            if (verbose) then
              print*,'Nb threads level 4:', omp_get_num_threads()
            endif

            !$OMP MASTER
            if (omp_get_num_threads()==1 .and. setting >= 0) then
              print*,'Setting',setting,"error at level 4"
            elseif(omp_get_num_threads()==1 .or. setting == 0) then
            else
              must_exit = .True.
            endif

            if ( z == 1 .and. setting == 0) then
              is_working = .True.
            elseif (z == 1 .and. setting == -1) then
              is_working = .False.
            else
            endif
            !$OMP END MASTER
  
  !          !$OMP DO
  !           do l = 1, x
  !             do k = 1, x
  !               do j = 1, x
  !                 do i = 1, x
  !                   accu(i,j,k,l) = accu(i,j,k,l)+ 1d0
  !                 enddo
  !               enddo
  !             enddo
  !           enddo
  !           !$OMP END DO       
   
          !$OMP END PARALLEL
  
        !$OMP END PARALLEL
    
      !$OMP END PARALLEL
      
    !$OMP END PARALLEL

  enddo
  
  print*,''

  if (setting == 1) then
    print*,'The parallelization works on 4 levels with:'
    print*,'call omp_set_max_active_levels(5)'
    print*,''
    print*,'Please use the irpf90 flags -DSET_MAX_ACT in qp2/config/${compiler_name}.cfg'
  elseif (setting == 2) then
    print*,'The parallelization works on 4 levels with:'
    print*,'call omp_set_nested(.True.)'
    print*,'' 
    print*,'Please use the irpf90 flag -DSET_NESTED in qp2/config/${compiler_name}.cfg'
  elseif (setting == 3) then
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
    print*,'Good luck...'
  endif

  if (is_working) then
    print*,''
    print*,'=========================================================='
    print*,'Your actual set up works for parallelization with 4 levels'
    print*,'=========================================================='
    print*,''
  else
    print*,''
    print*,'==================================================================='
    print*,'Your actual set up works for parallelization with 4 levels'
    print*,'Please look at the previous messages to understand the requirements'
    print*,'If it does not work even with the right irpf90 flags, clean and'
    print*,'recompile your code at ${QP_ROOT}'
    print*,'==================================================================='
    print*,''
  endif

end

