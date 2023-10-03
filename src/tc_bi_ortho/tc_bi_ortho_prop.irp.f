
! ---

program tc_bi_ortho_prop

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'Hello world'

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  !call routine_diag
  call test

end

! ---

subroutine test
 implicit none
 integer :: i
 print*,'TC Dipole components'
 do i= 1, 3
  print*,tc_bi_ortho_dipole(i,1)
 enddo
 integer, allocatable           :: occ(:,:)
 integer                        :: n_occ_ab(2)
 allocate(occ(N_int*bit_kind_size,2))
 call bitstring_to_list_ab(ref_bitmask, occ, n_occ_ab, N_int)
 integer :: ispin,j,jorb
 double precision :: accu
 accu = 0.d0
 do ispin=1, 2
  do i = 1, n_occ_ab(ispin)
   jorb = occ(i,ispin)
   accu += mo_bi_orth_bipole_z(jorb,jorb)
  enddo
 enddo
 print*,'accu = ',accu
 
end
