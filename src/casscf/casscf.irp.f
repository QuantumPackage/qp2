program casscf
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  no_vvvv_integrals = .True.
  SOFT_TOUCH no_vvvv_integrals 
  threshold_davidson = 1.d-7
  touch threshold_davidson 
  if(cisd_guess)then
    logical :: converged 
    integer :: iteration
    double precision :: energy
    print*,'*******************************'
    print*,'*******************************'
    print*,'*******************************'
    print*,'USING A CISD WAVE FUNCTION AS GUESS FOR THE MCSCF WF'
    print*,'*******************************'
    print*,'*******************************'
    converged = .False.
    iteration = 0
    generators_type = "HF"
    touch generators_type
    read_wf = .False.
    touch read_wf 
    logical :: do_cisdtq
    do_cisdtq  = .True.
    double precision :: thr
    thr = 5.d-3
    do while (.not.converged)
     call cisd_scf_iteration(converged,iteration,energy,thr)
     if(HF_index.ne.1.and.iteration.gt.0)then
      print*,'*******************************'
      print*,'*******************************'
      print*,'The HF determinant is not the dominant determinant in the CISD WF ...'
      print*,'Therefore we skip the CISD WF ..'
      print*,'*******************************'
      print*,'*******************************'
      do_cisdtq = .False.
      exit
     endif
     if(iteration.gt.15.and..not.converged)then
      print*,'It seems that the orbital optimization for the CISD WAVE FUNCTION CANNOT CONVERGE ...' 
      print*,'Passing to CISDTQ WAVE FUNCTION'
      exit
     endif
    enddo
    if(do_cisdtq)then
     print*,'*******************************'
     print*,'*******************************'
     print*,'*******************************'
     print*,'SWITCHING WITH A CISDTQ WAVE FUNCTION AS GUESS FOR THE MCSCF WF'
     print*,'*******************************'
     print*,'*******************************'
     converged = .False.
     iteration = 0
     read_wf = .False.
     touch read_wf 
     pt2_max = 0.01d0
     touch pt2_max 
     energy = 0.d0
     do while (.not.converged)
      call cisdtq_scf_iteration(converged,iteration,energy,thr)
      if(HF_index.ne.1.and.iteration.gt.0)then
       print*,'*******************************'
       print*,'*******************************'
       print*,'The HF determinant is not the dominant determinant in the CISDTQ WF ...'
       print*,'Therefore we skip the CISDTQ WF ..'
       print*,'*******************************'
       print*,'*******************************'
       exit
      endif
      if(iteration.gt.15.and..not.converged)then
       print*,'It seems that the orbital optimization for the CISDTQ WAVE FUNCTION CANNOT CONVERGE ...' 
       print*,'Passing to CISDTQ WAVE FUNCTION'
       exit
      endif
     enddo
    endif
  endif
  read_wf = .False.
  touch read_wf 
  pt2_max = 0.0d0
  touch pt2_max 
! call run_cipsi_scf
  call run 
end

subroutine run
  implicit none
  double precision               :: energy_old, energy
  logical                        :: converged
  integer                        :: iteration
  converged = .False.

  energy = 0.d0
  mo_label = "MCSCF"
  iteration = 1
  do while (.not.converged)
    call run_stochastic_cipsi
    energy_old = energy
    energy = eone+etwo+ecore

    call write_time(6)
    call write_int(6,iteration,'CAS-SCF iteration')
    call write_double(6,energy,'CAS-SCF energy')
    call write_double(6,energy_improvement, 'Predicted energy improvement')

    converged = dabs(energy_improvement) < thresh_scf
!   pt2_max = dabs(energy_improvement / pt2_relative_error)

    mo_coef = NewOrbs
    call save_mos
    iteration += 1
    N_det = N_det/2
    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    read_wf = .True.
    call clear_mo_map
    SOFT_TOUCH mo_coef N_det pt2_max  psi_det psi_coef 

  enddo

end
