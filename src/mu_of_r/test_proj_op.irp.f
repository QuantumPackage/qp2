program projected_operators
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
  ! You specify that you want to avoid any contribution from 
  ! orbitals coming from core 
  no_core_density = .True.
  touch no_core_density
  mu_of_r_potential = "cas_full"
  touch mu_of_r_potential 
  print*,'Using Valence Only functions'
!  call test_f_HF_valence_ab
!  call routine_full_mos
!   call test_f_ii_valence_ab
!   call test_f_ia_valence_ab
!  call test_f_ii_ia_aa_valence_ab
 call test
end


subroutine test
 implicit none
 integer :: i_point
 double precision :: ref, new, accu, weight
 accu = 0.d0
 do i_point = 1, n_points_final_grid
  ref = f_hf_cholesky_sparse(i_point)
  new = f_hf_cholesky_sparse_bis(i_point)
  weight = final_weight_at_r_vector(i_point)
  accu += dabs(ref - new) * weight
 enddo
 print*,'accu = ',accu

end
