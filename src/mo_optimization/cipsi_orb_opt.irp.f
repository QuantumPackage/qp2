program cipsi_orb_opt
  
  BEGIN_DOC
  ! Full orbital optimization in which the wave function grows as the
  ! optimization proceeds. Sets read_wf so that the wave function is
  ! read from the EZFIO file, then calls run_optimization_mos_CIPSI,
  ! which alternates stochastic CIPSI steps and orbital optimizations,
  ! doubling N_det_max after each cycle, until the wave function
  ! contains more than n_det_max_opt determinants.
  END_DOC
  read_wf = .true. ! must be True for the orbital optimization !!!
  TOUCH read_wf 
  call run_optimization_mos_CIPSI

end

