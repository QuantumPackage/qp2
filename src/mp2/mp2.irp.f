program mp2
  call run
end

subroutine run
  implicit none
  double precision, allocatable  :: pt2(:), norm_pert(:)
  double precision               :: H_pert_diag, E_old
  integer                        :: N_st, iter
  PROVIDE all_mo_integrals  Fock_matrix_diag_mo H_apply_buffer_allocated
  N_st = N_states
  allocate (pt2(N_st), norm_pert(N_st))
  E_old = HF_energy
  call H_apply_mp2(pt2, norm_pert, H_pert_diag,  N_st)
  print *,  'N_det    = ', N_det
  print *,  'N_states = ', N_states
  print *,  'MP2      = ', pt2
  print *,  'E        = ', E_old
  print *,  'E+MP2    = ', E_old+pt2
  deallocate(pt2,norm_pert)
end
