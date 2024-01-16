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
  if(small_active_space)then
   pt2_relative_error = 0.00001
  else
   thresh_scf = 1.d-4
   pt2_relative_error = 0.04
  endif
  touch pt2_relative_error 
  call run
end

subroutine run
  implicit none
  double precision               :: energy_old, energy, pt2_max_before,delta_E
  logical                        :: converged,state_following_casscf_cipsi_save
  integer                        :: iteration,istate
  double precision, allocatable :: E_PT2(:), PT2(:), Ev(:), ept2_before(:)
  allocate(E_PT2(N_states), PT2(N_states), Ev(N_states), ept2_before(N_states))
  converged = .False.

  energy = 0.d0
  mo_label = "MCSCF"
  iteration = 1
  state_following_casscf_cipsi_save = state_following_casscf
  state_following_casscf = .True.
  touch state_following_casscf
  ept2_before = 0.d0
  if(small_active_space)then
   pt2_max = 1.d-10
   SOFT_TOUCH pt2_max
  else
   if(adaptive_pt2_max)then
    pt2_max = 0.005
    SOFT_TOUCH pt2_max
   endif
  endif
  do while (.not.converged)
    print*,'pt2_max = ',pt2_max
    call run_stochastic_cipsi(Ev,PT2)
    print*,'Ev,PT2',Ev(1),PT2(1)
    E_PT2(1:N_states) = Ev(1:N_states) + PT2(1:N_states)
    energy_old = energy
    energy = eone+etwo+ecore
    pt2_max_before = pt2_max

    call write_time(6)
    call write_int(6,iteration,'CAS-SCF iteration = ')
    call write_double(6,energy,'CAS-SCF energy = ')
!    if(n_states == 1)then
!     call ezfio_get_casscf_cipsi_energy_pt2(E_PT2)
!     call ezfio_get_casscf_cipsi_energy(PT2)
     do istate=1,N_states
     call write_double(6,E_PT2(istate),'E + PT2 energy = ')
     call write_double(6,PT2(istate),'  PT2          = ')
     enddo
     call write_double(6,pt2_max,' PT2_MAX       = ')
!    endif

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
     delta_E = 0.d0
     do istate = 1, N_states
      delta_E += dabs(E_PT2(istate) - ept2_before(istate))
     enddo
     converged = dabs(delta_E) < thresh_casscf
    endif
    ept2_before = E_PT2
    if(.not.small_active_space)then
     if(adaptive_pt2_max)then
      pt2_max = dabs(energy_improvement / (pt2_relative_error))
      pt2_max = min(pt2_max, pt2_max_before)
      if(n_act_orb.ge.n_big_act_orb)then
       pt2_max = max(pt2_max,pt2_min_casscf)
      endif
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
     if(.not.small_active_space)then
      if(adaptive_pt2_max)then
        SOFT_TOUCH pt2_max  
      endif
     endif
     if(iteration .gt. 3)then
      state_following_casscf = state_following_casscf_cipsi_save
      soft_touch state_following_casscf
     endif
    endif

  enddo
     integer :: i
    print*,'Converged CASSCF '
    print*,'--------------------------'
    write(6,*) ' occupation numbers of orbitals '
    do i=1,mo_num
      write(6,*) i,occnum(i)
    end do
    print*,'--------------'
!
!     write(6,*)
!     write(6,*) ' the diagonal of the inactive effective Fock matrix '
!     write(6,'(5(i3,F12.5))') (i,Fipq(i,i),i=1,mo_num)
!     write(6,*)
  print*,'Fock MCSCF'
  do i = 1, mo_num
   write(*,*)i,mcscf_fock_diag_mo(i)
!   write(*,*)mcscf_fock_alpha_mo(i,i)
  enddo


end


