
BEGIN_PROVIDER [ logical, ao_two_e_integrals_erf_in_map ]
  implicit none
  use f77_zmq
  use map_module
  BEGIN_DOC
  !  Map of Atomic integrals
  !     i(r1) j(r2) 1/r12 k(r1) l(r2)
  END_DOC

  integer                        :: i,j,k,l
  double precision               :: ao_two_e_integral_erf,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  include 'utils/constants.include.F'

  ! For integrals file
  integer(key_kind),allocatable  :: buffer_i(:)
  integer,parameter              :: size_buffer = 1024*64
  real(integral_kind),allocatable :: buffer_value(:)

  integer                        :: n_integrals, rc
  integer                        :: kk, m, j1, i1, lmax
  character*(64)                 :: fmt

  integral = ao_two_e_integral_erf(1,1,1,1)

  double precision               :: map_mb
  PROVIDE read_ao_two_e_integrals_erf io_ao_two_e_integrals_erf
  if (read_ao_two_e_integrals_erf) then
    print*,'Reading the AO ERF integrals'
      call map_load_from_disk(trim(ezfio_filename)//'/work/ao_ints_erf',ao_integrals_erf_map)
      print*, 'AO ERF integrals provided'
      ao_two_e_integrals_erf_in_map = .True.
      return
  endif

  print*, 'Providing the AO ERF integrals'
  call wall_time(wall_0)
  call wall_time(wall_1)
  call cpu_time(cpu_1)

  integer(ZMQ_PTR) :: zmq_to_qp_run_socket, zmq_socket_pull
  call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,'ao_integrals_erf')

  character(len=:), allocatable :: task
  allocate(character(len=ao_num*12) :: task)
  write(fmt,*) '(', ao_num, '(I5,X,I5,''|''))'
  do l=1,ao_num
    write(task,fmt) (i,l, i=1,l)
    integer, external :: add_task_to_taskserver
    if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task)) == -1) then
      stop 'Unable to add task to server'
    endif
  enddo
  deallocate(task)

  integer, external :: zmq_set_running
  if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
    print *,  irp_here, ': Failed in zmq_set_running'
  endif

  PROVIDE nproc
  !$OMP PARALLEL DEFAULT(shared) private(i) num_threads(nproc+1)
      i = omp_get_thread_num()
      if (i==0) then
        call ao_two_e_integrals_erf_in_map_collector(zmq_socket_pull)
      else
        call ao_two_e_integrals_erf_in_map_slave_inproc(i)
      endif
  !$OMP END PARALLEL

  call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'ao_integrals_erf')


  print*, 'Sorting the map'
  call map_sort(ao_integrals_erf_map)
  call cpu_time(cpu_2)
  call wall_time(wall_2)
  integer(map_size_kind)         :: get_ao_erf_map_size, ao_erf_map_size
  ao_erf_map_size = get_ao_erf_map_size()

  print*, 'AO ERF integrals provided:'
  print*, ' Size of AO ERF map :         ', map_mb(ao_integrals_erf_map) ,'MB'
  print*, ' Number of AO ERF integrals :', ao_erf_map_size
  print*, ' cpu  time :',cpu_2 - cpu_1, 's'
  print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'

  ao_two_e_integrals_erf_in_map = .True.

  if (write_ao_two_e_integrals_erf) then
    call ezfio_set_work_empty(.False.)
    call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_erf',ao_integrals_erf_map)
    call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals_erf('Read')
  endif

END_PROVIDER




BEGIN_PROVIDER [ double precision, ao_two_e_integral_erf_schwartz,(ao_num,ao_num)  ]
  implicit none
  BEGIN_DOC
  !  Needed to compute Schwartz inequalities
  END_DOC

  integer                        :: i,k
  double precision               :: ao_two_e_integral_erf,cpu_1,cpu_2, wall_1, wall_2

  ao_two_e_integral_erf_schwartz(1,1) = ao_two_e_integral_erf(1,1,1,1)
  !$OMP PARALLEL DO PRIVATE(i,k)                                     &
      !$OMP DEFAULT(NONE)                                            &
      !$OMP SHARED (ao_num,ao_two_e_integral_erf_schwartz)              &
      !$OMP SCHEDULE(dynamic)
  do i=1,ao_num
    do k=1,i
      ao_two_e_integral_erf_schwartz(i,k) = dsqrt(ao_two_e_integral_erf(i,k,i,k))
      ao_two_e_integral_erf_schwartz(k,i) = ao_two_e_integral_erf_schwartz(i,k)
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER


