program cipsi_orb_opt
  
  read_wf = .true. ! must be True for the orbital optimization !!!
  TOUCH read_wf 
  call run_optimization_mos_CIPSI

end

