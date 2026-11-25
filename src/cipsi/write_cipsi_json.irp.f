subroutine write_cipsi_json(pt2_data, pt2_data_err)
    use selection_types
    use qp_json, only: json_add
    implicit none
    BEGIN_DOC
!  Logs the CIPSI results of the current iteration.
    END_DOC
    type(pt2_type), intent(in) :: pt2_data, pt2_data_err
    integer :: n
    double precision, allocatable :: buf(:)

    n = min(N_states,N_det)
    if (n < 1) return

    call json_add('/results/cipsi/n_det', N_det)
    if (s2_eig) call json_add('/results/cipsi/n_cfg', N_Configuration)
    call json_add('/results/cipsi/states/energy', psi_energy_with_nucl_rep(1:n))
    call json_add('/results/cipsi/states/s2', psi_s2(1:n))
    call json_add('/results/cipsi/states/pt2', pt2_data % pt2(1:n))
    call json_add('/results/cipsi/states/pt2_err', pt2_data_err % pt2(1:n))
    call json_add('/results/cipsi/states/rpt2', pt2_data % rpt2(1:n))
    call json_add('/results/cipsi/states/rpt2_err', pt2_data_err % rpt2(1:n))
    call json_add('/results/cipsi/states/variance', pt2_data % variance(1:n))
    call json_add('/results/cipsi/states/variance_err', &
                  pt2_data_err % variance(1:n))
    if (N_iter > 0) then
      call json_add('/results/cipsi/history/n_iter', N_iter)
      call json_add('/results/cipsi/history/n_iter_max', N_iter_max)
      call json_add('/results/cipsi/states/extrapolated_energy', &
                    extrapolated_energy(N_iter,1:n))
      call json_add('/results/cipsi/history/shape', (/ n, N_iter /))
      allocate(buf(n*N_iter))
      buf = reshape(energy_iterations(1:n,1:N_iter), (/ n*N_iter /))
      call json_add('/results/cipsi/history/energy', buf)
      buf = reshape(pt2_iterations(1:n,1:N_iter), (/ n*N_iter /))
      call json_add('/results/cipsi/history/rpt2', buf)
      deallocate(buf)
    endif
end
