
subroutine run_optimization_mos_CIPSI

  implicit none

  double precision :: e_cipsi, e_opt, delta_e
  double precision, allocatable :: Ev(:),PT2(:)
  integer :: nb_iter,i
  logical :: not_converged
  character (len=100) :: filename

  PROVIDE psi_det psi_coef mo_two_e_integrals_in_map ao_pseudo_integrals
  allocate(Ev(N_states),PT2(N_states))

  not_converged = .True.
  nb_iter = 0

  ! To start from the wf
  N_det_max = max(n_det,5)
  TOUCH N_det_max

  open(unit=10, file=trim(ezfio_filename)//'/mo_optimization/result_opt')
  write(10,*) "   Ndet        E_cipsi         E_opt          Delta_e"
  call state_average_energy(e_cipsi)
  write(10,'(I10, 3F15.7)') n_det, e_cipsi, e_cipsi, 0d0
  close(10)

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

      open(unit=10, file=trim(ezfio_filename)//'/mo_optimization/result_opt', position='append')
      write(10,'(I10, 3F15.7)') n_det, e_cipsi, e_opt, delta_e
      close(10)

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
