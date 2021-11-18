program test_intel_check_omp
  
  use omp_lib

  implicit none

  integer :: i,j,k,l,m,n,x
  double precision :: w1,w2,c1,c2
  double precision, allocatable :: accu(:,:,:,:) 

  x = 4
  allocate(accu(x,x,x,x))

  accu = 0d0

  !$OMP PARALLEL
  print*, 'Hello1 from:', OMP_GET_THREAD_NUM()
  !$OMP END PARALLEL

  print*,'omp_get_max_active_levels:',omp_get_max_active_levels()
  call intel_check_omp()
  print*,'omp_get_max_active_levels:',omp_get_max_active_levels()

  !call omp_set_max_active_levels(20000)

  !$OMP PARALLEL
  print*, 'Hello2 from:', OMP_GET_THREAD_NUM()
  !$OMP END PARALLEL

  call wall_time(w1)
  call cpu_time(c1) 
  !$OMP PARALLEL &
    !$OMP PRIVATE(i,j,k,l,m,n) &
    !$OMP SHARED(accu) 

    print*,'level 1',omp_get_num_threads()
    !$OMP DO
    do l = 1, x
      do k = 1, x
        do j = 1, x
          do i = 1, x
            accu(i,j,k,l) = accu(i,j,k,l) + 1d0
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP PARALLEL & 
      !$OMP PRIVATE(i,j,k,l,m,n) &
      !$OMP SHARED(accu)
        
      print*,'level 2',omp_get_num_threads()
      !$OMP DO
      do l = 1, x
        do k = 1, x
          do j = 1, x
            do i = 1, x
              accu(i,j,k,l) = accu(i,j,k,l)+ 1d0
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO

      !$OMP PARALLEL & 
        !$OMP PRIVATE(i,j,k,l,m,n) &
        !$OMP SHARED(accu)

        print*,'level 3',omp_get_num_threads()
        !$OMP DO
        do l = 1, x
          do k = 1, x
            do j = 1, x
              do i = 1, x
                accu(i,j,k,l) = accu(i,j,k,l)+ 1d0
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO

        !$OMP PARALLEL &
            !$OMP PRIVATE(i,j,k,l,m,n) &
            !$OMP SHARED(accu)

            print*,'level 4',omp_get_num_threads()
            !$OMP DO
             do l = 1, x
               do k = 1, x
                 do j = 1, x
                   do i = 1, x
                     accu(i,j,k,l) = accu(i,j,k,l)+ 1d0
                   enddo
                 enddo
               enddo
             enddo
             !$OMP END DO       
 
        !$OMP END PARALLEL

      !$OMP END PARALLEL
  
    !$OMP END PARALLEL
    
  !$OMP END PARALLEL

  call wall_time(w2)
  call cpu_time(c2)

  print*,accu(1,1,1,1)
  print*,'wall time:', w2-w1
  print*,'cpu time:', c2-c1
  print*,'ration',(c2-c1)/(w2-w1)
end
