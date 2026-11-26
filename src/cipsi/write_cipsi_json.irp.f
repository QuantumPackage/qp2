subroutine write_cipsi_json(pt2_data, pt2_data_err)
    use selection_types
    use qp_json, only: json_add
    use c_functions, only: qp_epoch_c
    implicit none
    BEGIN_DOC
!  Logs the CIPSI results of the current iteration.
    END_DOC
    ! One row per call, appended to SAVE arrays; the JSON document is
    ! rewritten from scratch each call (json_add replaces).
    type(pt2_type), intent(in) :: pt2_data, pt2_data_err

    integer, save :: nst = 0, nh = 0
    double precision, allocatable, save :: h_energy(:), h_pt2(:), h_rpt2(:), &
        h_var(:), h_perr(:), h_samp(:), h_epoch(:), h_wall(:), h_mem(:), &
        h_ptime(:), h_pstart(:), h_pend(:)
    integer, allocatable, save :: h_ndet(:), h_ncfg(:), h_win(:), h_task(:)

    integer :: n
    double precision :: wt, mem

    n = min(N_states,N_det)
    if (n < 1) return

    if (nst /= n) then                 ! a different CI space: new series
        if (nh > 0) deallocate(h_energy, h_pt2, h_rpt2, h_var, h_perr, &
             h_samp, h_epoch, h_wall, h_mem, h_ptime, h_pstart, h_pend, &
             h_ndet, h_ncfg, h_win, h_task)
        nh = 0
        nst = n
    endif
    if (nh == 0) allocate(h_energy(0), h_pt2(0), h_rpt2(0), h_var(0), &
        h_perr(0), h_samp(0), h_epoch(0), h_wall(0), h_mem(0), h_ptime(0), &
        h_pstart(0), h_pend(0), h_ndet(0), h_ncfg(0), h_win(0), h_task(0))

    nh = nh + 1
    h_energy = [h_energy, psi_energy_with_nucl_rep(1:n)]
    h_pt2    = [h_pt2,    pt2_data % pt2(1:n)]
    h_rpt2   = [h_rpt2,   pt2_data % rpt2(1:n)]
    h_var    = [h_var,    pt2_data % variance(1:n)]
    h_perr   = [h_perr,   pt2_data_err % pt2(1:n)]
    h_samp   = [h_samp,   pt2_data % samples]      ! scalar, broadcast by shape
    h_ndet   = [h_ndet,   N_det]
    h_ncfg   = [h_ncfg,   merge(N_Configuration, 0, s2_eig)]
    h_win    = [h_win,    pt2_data % windows]
    h_task   = [h_task,   pt2_data % tasks]
    h_ptime  = [h_ptime,  pt2_data % time]
    h_pstart = [h_pstart, pt2_data % epoch_start]
    h_pend   = [h_pend,   pt2_data % epoch_end]
    h_epoch  = [h_epoch,  qp_epoch_c()]
    call wall_time(wt)
    h_wall   = [h_wall,   wt - human_t0]
    call resident_memory(mem)
    h_mem    = [h_mem,    mem]

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
    endif
    call json_add('/results/cipsi/history/shape', (/ n, nh /))
    call json_add('/results/cipsi/history/energy', h_energy)
    call json_add('/results/cipsi/history/pt2', h_pt2)
    call json_add('/results/cipsi/history/rpt2', h_rpt2)
    call json_add('/results/cipsi/history/variance', h_var)
    call json_add('/results/cipsi/history/pt2_err', h_perr)
    call json_add('/results/cipsi/history/n_det', h_ndet)
    if (s2_eig) call json_add('/results/cipsi/history/n_cfg', h_ncfg)
    call json_add('/results/cipsi/history/pt2_windows', h_win)
    call json_add('/results/cipsi/history/pt2_tasks', h_task)
    call json_add('/results/cipsi/history/samples', h_samp)
    call json_add('/results/cipsi/history/pt2_time', h_ptime)
    call json_add('/results/cipsi/history/pt2_epoch_start', h_pstart)
    call json_add('/results/cipsi/history/pt2_epoch_end', h_pend)
    call json_add('/results/cipsi/history/epoch', h_epoch)
    call json_add('/results/cipsi/history/wall_time', h_wall)
    call json_add('/results/cipsi/history/memory', h_mem)
end
