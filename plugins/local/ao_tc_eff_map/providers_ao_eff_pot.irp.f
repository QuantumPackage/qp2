
BEGIN_PROVIDER [ logical, ao_tc_sym_two_e_pot_in_map ]
  implicit none
  use omp_lib
  use map_module
  BEGIN_DOC
  !  Map of Atomic integrals
  !     i(r1) j(r2) 1/r12 k(r1) l(r2)
  !
  ! In-process port (2026-11-12, the distributed layer was removed):
  ! the former task-server
  ! dispatch (one task per column l = 1..ao_num, nproc workers each
  ! computing the first pair of the task string and pushing the result
  ! to a collector that appended to the map) is replaced by one OpenMP
  ! team of nproc threads running a parallel-do over l, each thread
  ! computing with its own buffer and appending to the map under a
  ! named critical region (the map is not thread-safe; the legacy
  ! serialized the same append through the collector). The computed
  ! set is exactly the legacy one: for each l, compute_ao_tc_sym_two_e_
  ! pot_jl is called on the first pair of the column task, i.e.
  ! (j,l) = (1,l) [legacy quirk, kept faithfully: the worker read
  ! only the first pair of each task string].
  END_DOC

  integer                        :: i,j,k,l
  double precision               :: ao_tc_sym_two_e_pot,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  include 'utils/constants.include.F'

  ! For integrals file
  integer(key_kind),allocatable  :: buffer_i(:)
  integer,parameter              :: size_buffer = 1024*64
  real(integral_kind),allocatable :: buffer_value(:)

  integer                        :: n_integrals, rc
  integer                        :: kk, m, j1, i1, lmax
  character*(64)                 :: fmt

  integral = ao_tc_sym_two_e_pot(1,1,1,1)

  double precision               :: map_mb

  print*, 'Providing the ao_tc_sym_two_e_pot_map integrals'
  call wall_time(wall_0)
  call wall_time(wall_1)
  call cpu_time(cpu_1)

  PROVIDE nproc

  ! One buffer per team member; each thread computes one (1,l) column
  ! at a time and appends the result to the shared map under a critical
  ! region (the legacy collector did the same appends, serialized by
  ! the arrival order of the old distributed workers).
  integer*8 :: accu
  accu = 0_8

  ! Serial warm-up on the first column: IRPF90 providers are lazy and
  ! not thread-safe, and compute_ao_tc_sym_two_e_pot_jl reads
  ! providers (ao_two_e_integral_erf and friends) that are built on
  ! first access.  The column is deterministic and is skipped in the
  ! parallel-do because its result is already in the map (the
  ! pt2_omp pattern).
  allocate ( buffer_i(ao_num*ao_num), buffer_value(ao_num*ao_num) )
  call compute_ao_tc_sym_two_e_pot_jl(1,1,n_integrals,buffer_i,buffer_value)
  call insert_into_ao_tc_sym_two_e_pot_map(n_integrals,buffer_i,buffer_value)
  accu = accu + n_integrals
  deallocate ( buffer_i, buffer_value )

  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nproc) PRIVATE(i, n_integrals)
     i = omp_get_thread_num()
     integer(key_kind), allocatable :: buffer_i_t(:)
     real(integral_kind), allocatable :: buffer_value_t(:)
     allocate ( buffer_i_t(ao_num*ao_num), buffer_value_t(ao_num*ao_num) )
     !$OMP DO SCHEDULE(DYNAMIC,1)
     do l=2,ao_num
       call compute_ao_tc_sym_two_e_pot_jl(1,l,n_integrals,buffer_i_t,buffer_value_t)
       !$OMP CRITICAL
         call insert_into_ao_tc_sym_two_e_pot_map(n_integrals,buffer_i_t,buffer_value_t)
         accu = accu + n_integrals
       !$OMP END CRITICAL
     enddo
     !$OMP END DO
     deallocate ( buffer_i_t, buffer_value_t )
  !$OMP END PARALLEL

  integer(map_size_kind)         :: get_ao_tc_sym_two_e_pot_map_size, ao_eff_pot_map_size
  ao_eff_pot_map_size = get_ao_tc_sym_two_e_pot_map_size()

  if (ao_eff_pot_map_size /= accu) then
      print *, ''
      print *, irp_here
      print *, 'Control : ', ao_eff_pot_map_size
      print *, 'Accu    : ', accu
      print *, 'Some integrals were lost during the parallel computation.'
      print *, 'Try to reduce the number of threads.'
      stop
  endif

  print*, 'Sorting the map'
  call map_sort(ao_tc_sym_two_e_pot_map)
  call cpu_time(cpu_2)
  call wall_time(wall_2)

  print*, 'AO eff_pot integrals provided:'
  print*, ' Size of AO eff_pot map :         ', map_mb(ao_tc_sym_two_e_pot_map) ,'MB'
  print*, ' Number of AO eff_pot integrals :', ao_eff_pot_map_size
  print*, ' cpu  time :',cpu_2 - cpu_1, 's'
  print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'

  ao_tc_sym_two_e_pot_in_map = .True.


END_PROVIDER
