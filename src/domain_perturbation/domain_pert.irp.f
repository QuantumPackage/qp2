program domain_pert
  implicit none
  BEGIN_DOC
  ! Perturbation theory derived from the domain Dyson equation, for all
  ! the states stored in the |EZFIO| directory.
  !
  ! The reference space is the set of determinants of the wave
  ! function. The orders 2 to :option:`domain_perturbation n_pert` of
  ! the energy are computed by a stochastic algorithm with
  ! annihilation, or exactly if
  ! :option:`domain_perturbation n_walkers` is 0.
  END_DOC
  read_wf = .True.
  SOFT_TOUCH read_wf
  PROVIDE all_mo_integrals
  call run_domain_pert
end

subroutine run_domain_pert
  implicit none
  BEGIN_DOC
  ! Runs the calculation for all the states and prints the results.
  END_DOC
  double precision, allocatable :: e(:), err(:), errt(:), ecb(:)
  double precision :: etot(N_states), eerr(N_states), s, t0, t1
  double precision :: eg, a, b
  integer :: k, p
  logical :: auto, positive

  if (n_pert < 2) then
    stop 'n_pert should be at least 2'
  endif
  allocate(e(2:n_pert),err(2:n_pert),errt(2:n_pert),ecb(2:n_pert))

  call write_int(6,N_det,'Number of determinants in the reference')
  call write_int(6,n_pert,'Maximum order')
  call write_int(6,n_det_domain,'Maximum size of the external domains')
  call write_int(6,n_walkers,'Number of walkers (0 : exact)')

  do k=1,N_states
    if (domain_residual(k) > 1.d-4) then
      print *, 'Warning: the reference state', k, ' is not converged.'
      print *, 'Norm of the residual : ', domain_residual(k)
    endif
    call wall_time(t0)
    if (domain_chebyshev) then
      call domain_state_cheb(k,eg,a,b,auto,positive,e,err,ecb,errt)
    else
      call domain_state(k,e,err,errt)
    endif
    call wall_time(t1)
    print *, ''
    print '(A,I4)', ' State ', k
    print '(A)', ' ============='
    print '(A,F20.12)', ' E_0         = ', domain_e0(k)+nuclear_repulsion
    if (domain_chebyshev) then
      print '(A,F20.12)', ' E_0 + E_2   = ', eg+nuclear_repulsion
      print '(A,2F10.5,A,L2)', ' Interval    = ', a, b,               &
          '   automatic:', auto
      if (.not.positive) then
        print *, 'Warning: the moments are not those of a positive'
        print *, 'measure (complex spectrum or noise). The Chebyshev'
        print *, 'series can be unreliable for this state.'
      endif
      print '(A4,A22,A14,A22,A14)', 'p', 'Brillouin-Wigner',         &
          'error', 'Chebyshev', 'error'
      do p=2,n_pert
        print '(I4,F22.12,E14.4,F22.12,E14.4)', p,                   &
            e(p)+nuclear_repulsion, err(p),                          &
            ecb(p)+nuclear_repulsion, errt(p)
      enddo
      etot(k) = ecb(n_pert)+nuclear_repulsion
    else
      print '(A4,A22,A14,A22,A14)', 'p', 'E_p', 'error',             &
          'E_0+...+E_p', 'error'
      s = domain_e0(k)+nuclear_repulsion
      do p=2,n_pert
        s = s + e(p)
        print '(I4,F22.12,E14.4,F22.12,E14.4)', p, e(p), err(p),     &
            s, errt(p)
      enddo
      etot(k) = s
    endif
    eerr(k) = errt(n_pert)
    call write_double(6,t1-t0,'Wall time (s)')
  enddo
  call write_double(6,domain_min_gap_seen,'Smallest |H_jj - E_0|')
  call ezfio_set_domain_perturbation_domain_energy(etot)
  call ezfio_set_domain_perturbation_domain_energy_error(eerr)
  deallocate(e,err,errt,ecb)
end
