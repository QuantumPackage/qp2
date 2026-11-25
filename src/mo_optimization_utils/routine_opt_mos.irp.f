
subroutine run_optimization_mos_CIPSI
  BEGIN_DOC
  ! Orbital optimization with a wave function that grows as the
  ! optimization proceeds: each cycle runs a stochastic CIPSI step
  ! (run_stochastic_cipsi) then an orbital optimization
  ! (run_orb_opt_trust_v2), and doubles N_det_max. Each cycle writes
  ! n_det, the energy before and after the optimization and their
  ! difference to /results/mo_opt of the JSON result document.
  END_DOC

  use qp_json, only: json_add
  implicit none

  double precision :: e_cipsi, e_opt, delta_e
  double precision, allocatable :: Ev(:),PT2(:)
  integer :: nb_iter,i
  logical :: not_converged

  call json_init('cipsi_orb_opt')

  PROVIDE psi_det psi_coef all_mo_integrals ao_pseudo_integrals
  allocate(Ev(N_states),PT2(N_states))

  not_converged = .True.
  nb_iter = 0

  ! To start from the wf
  N_det_max = max(n_det,5)
  TOUCH N_det_max

  ! State average energy of the wave function we start from
  call state_average_energy(e_cipsi)

  do while (not_converged)
      print*,''
      print*,'======================'
      print*,' Cipsi step:', nb_iter
      print*,'======================'
      print*,''
      print*,'********** cipsi step **********'
      ! cispi calculation
      call run_stochastic_cipsi(Ev,PT2)

      ! State average energy after the cipsi step
      call state_average_energy(e_cipsi)

      print*,''
      print*,'********** optimization step **********'
      ! orbital optimization
      call run_orb_opt_trust_v2

      ! State average energy after the orbital optimization
      call state_average_energy(e_opt)

      print*,''
      print*,'********** diff step **********'
      ! Gain in energy
      delta_e = e_opt - e_cipsi
      print*, 'Gain in energy during the orbital optimization:', delta_e

      call json_add('/results/mo_opt/iteration', nb_iter)
      call json_add('/results/mo_opt/n_det', n_det)
      call json_add('/results/mo_opt/energy_cipsi', e_cipsi)
      call json_add('/results/mo_opt/energy_opt', e_opt)
      call json_add('/results/mo_opt/delta_energy', delta_e)

      ! Exit
      if (delta_e > 1d-12) then
          print*, 'WARNING, something wrong happened'
          print*, 'The gain (delta_e) in energy during the optimization process'
          print*, 'is > 0, but it must be < 0'
          print*, 'The program will exit'
          exit
      endif

      if (n_det > n_det_max_opt) then
          print*, 'The number of determinants in the wf > n_det_max_opt'
          print*, 'The program will exit'
          exit
      endif

      ! To double the number of determinants in the wf
      N_det_max = int(dble(n_det * 2)*0.9)
      TOUCH N_det_max

      nb_iter = nb_iter + 1
  enddo

end
