! ---

program tc_som

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, ' starting ...'
  print *, ' do not forget to do tc-scf first'

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
!  my_n_pt_r_grid = 10 ! small grid for quick debug
!  my_n_pt_a_grid = 26 ! small grid for quick debug
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  PROVIDE mu_erf 
  print *, ' mu = ', mu_erf
  PROVIDE j1b_type
  print *, ' j1b_type = ', j1b_type
  print *, j1b_pen

  read_wf = .true.
  touch read_wf

  call main()

end

! ---

subroutine main()

  implicit none
  integer          :: i, i_HF, degree
  double precision :: hmono_1, htwoe_1, hthree_1, htot_1
  double precision :: hmono_2, htwoe_2, hthree_2, htot_2
  double precision :: U_SOM

  PROVIDE N_int N_det

  do i = 1, N_det
    call get_excitation_degree(HF_bitmask, psi_det(1,1,i), degree, N_int)
    if(degree == 0) then
      i_HF = i
      exit
    endif
  enddo
  print *, ' HF determinants:', i_HF
  print *, '          N_det :', N_det

  U_SOM = 0.d0 
  do i = 1, N_det
    if(i == i_HF) cycle
    call htilde_mu_mat_bi_ortho(psi_det(1,1,i_HF), psi_det(1,1,i), N_int, hmono_1, htwoe_1, hthree_1, htot_1)
    call htilde_mu_mat_bi_ortho(psi_det(1,1,i), psi_det(1,1,i_HF), N_int, hmono_2, htwoe_2, hthree_2, htot_2)
    U_SOM += htot_1 * htot_2
  enddo
  U_SOM = 0.5d0 * U_SOM
  print *, ' U_SOM = ', U_SOM
  
  return
end subroutine main

! ---

