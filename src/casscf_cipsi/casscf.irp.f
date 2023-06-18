program casscf
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  call reorder_orbitals_for_casscf
!  no_vvvv_integrals = .True.
!  touch no_vvvv_integrals
  n_det_max_full = 500
  touch n_det_max_full
  pt2_relative_error = 0.04
  touch pt2_relative_error 
!  call run_stochastic_cipsi
  call run
end

subroutine run
  implicit none
  double precision               :: energy_old, energy, pt2_max_before, ept2_before,delta_E
  logical                        :: converged,state_following_casscf_save
  integer                        :: iteration
  converged = .False.

  energy = 0.d0
  mo_label = "MCSCF"
  iteration = 1
  state_following_casscf_save = state_following_casscf
  state_following_casscf = .True.
  touch state_following_casscf
  ept2_before = 0.d0
  if(adaptive_pt2_max)then
   pt2_max = 0.005
   SOFT_TOUCH pt2_max
  endif
  do while (.not.converged)
    print*,'pt2_max = ',pt2_max
    call run_stochastic_cipsi
    energy_old = energy
    energy = eone+etwo+ecore
    pt2_max_before = pt2_max

    call write_time(6)
    call write_int(6,iteration,'CAS-SCF iteration = ')
    call write_double(6,energy,'CAS-SCF energy = ')
    if(n_states == 1)then
     double precision :: E_PT2, PT2
     call ezfio_get_casscf_energy_pt2(E_PT2)
     call ezfio_get_casscf_energy(PT2)
     PT2 -= E_PT2
     call write_double(6,E_PT2,'E + PT2 energy = ')
     call write_double(6,PT2,'  PT2          = ')
     call write_double(6,pt2_max,' PT2_MAX       = ')
    endif

    print*,''
    call write_double(6,norm_grad_vec2,'Norm of gradients = ')
    call write_double(6,norm_grad_vec2_tab(1), '     Core-active  gradients = ')
    call write_double(6,norm_grad_vec2_tab(2), '     Core-virtual gradients = ')
    call write_double(6,norm_grad_vec2_tab(3), '   Active-virtual gradients = ')
    print*,''
    call write_double(6,energy_improvement, 'Predicted energy improvement = ')

    if(criterion_casscf == "energy")then
     converged = dabs(energy_improvement) < thresh_scf
    else if (criterion_casscf == "gradients")then
     converged = norm_grad_vec2 < thresh_scf
    else if (criterion_casscf == "e_pt2")then
     delta_E = dabs(E_PT2 - ept2_before)
     converged = dabs(delta_E) < thresh_casscf
    endif
    ept2_before = E_PT2
    if(adaptive_pt2_max)then
     pt2_max = dabs(energy_improvement / (pt2_relative_error))
     pt2_max = min(pt2_max, pt2_max_before)
     if(n_act_orb.ge.n_big_act_orb)then
      pt2_max = max(pt2_max,pt2_min_casscf)
     endif
    endif
    print*,''
    call write_double(6,pt2_max, 'PT2_MAX for next iteration = ')

    mo_coef = NewOrbs
    mo_occ  = occnum
    call save_mos
    if(.not.converged)then
     iteration += 1
     if(norm_grad_vec2.gt.0.01d0)then
      N_det = N_states
     else
      N_det = max(N_det/8 ,N_states)
     endif
     psi_det = psi_det_sorted
     psi_coef = psi_coef_sorted
     read_wf = .True.
     call clear_mo_map
     SOFT_TOUCH mo_coef N_det psi_det psi_coef
     if(adaptive_pt2_max)then
       SOFT_TOUCH pt2_max  
     endif
     if(iteration .gt. 3)then
      state_following_casscf = state_following_casscf_save
      soft_touch state_following_casscf
     endif
    endif

  enddo

end


