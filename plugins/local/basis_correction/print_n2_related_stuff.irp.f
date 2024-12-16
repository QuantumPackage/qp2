program print_n2_stuffs
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf
  no_core_density = .True.
  touch no_core_density
  call routine
end

subroutine routine
 implicit none
 print*,'average_extrapolated_on_top = ',average_extrapolated_on_top
 print*,'average_on_top              = ',average_on_top
 print*,'mu_average_prov             = ',mu_average_prov
end
