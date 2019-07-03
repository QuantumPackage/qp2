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
  generators_type = "CAS"
  touch generators_type
  read_wf = .False.
  touch read_wf 
  pt2_max = 0.015d0
  touch pt2_max 
  call run_cipsi_scf
end

