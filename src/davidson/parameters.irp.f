BEGIN_PROVIDER [ character(64), davidson_criterion ]
 implicit none
 BEGIN_DOC
 ! Can be : [  energy  | residual | both | wall_time | cpu_time | iterations ]
 END_DOC
 davidson_criterion = 'residual'
END_PROVIDER

subroutine davidson_converged(energy,residual,wall,iterations,cpu,N_st,converged)
  implicit none
  BEGIN_DOC
! True if the Davidson algorithm is converged
  END_DOC
  integer, intent(in) :: N_st, iterations
  logical, intent(out) :: converged
  double precision, intent(in) :: energy(N_st), residual(N_st)
  double precision, intent(in) :: wall, cpu
  double precision :: E(N_st), time
  double precision, allocatable, save :: energy_old(:)

  if (iterations < 2) then
    converged = .False.
    return
  endif

  if (.not.allocated(energy_old)) then
    allocate(energy_old(N_st))
    energy_old = 0.d0
  endif

  E = energy - energy_old
  energy_old = energy
  if (davidson_criterion == 'energy') then
    converged = dabs(maxval(E(1:N_st))) < threshold_davidson
  else if (davidson_criterion == 'residual') then
    converged = dabs(maxval(residual(1:N_st))) < threshold_davidson
  else if (davidson_criterion == 'both') then
    converged = dabs(maxval(residual(1:N_st))) + dabs(maxval(E(1:N_st)) ) &
       < threshold_davidson
  else if (davidson_criterion == 'wall_time') then
    call wall_time(time)
    converged = time - wall > threshold_davidson
  else if (davidson_criterion == 'cpu_time') then
    call cpu_time(time)
    converged = time - cpu > threshold_davidson
  else if (davidson_criterion == 'iterations') then
    converged = iterations >= int(threshold_davidson)
  endif
end
