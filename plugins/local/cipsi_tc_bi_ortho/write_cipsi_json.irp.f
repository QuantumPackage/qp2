subroutine write_cipsi_json(pt2_data, pt2_data_err)
    use selection_types
    implicit none
    BEGIN_DOC
!  Writes JSON data for CIPSI runs
    END_DOC
    type(pt2_type), intent(in) :: pt2_data, pt2_data_err
    integer :: i,j,k

    call lock_io
    character*(64), allocatable    :: fmtk(:)
    integer :: N_states_p, N_iter_p
    N_states_p = min(N_states,N_det)
    N_iter_p = min(N_iter,8)
    allocate(fmtk(0:N_iter_p))
    fmtk(:) = '(''   '',E22.15,'','')'
    fmtk(N_iter_p) = '(''   '',E22.15)'

    write(json_unit, json_dict_uopen_fmt)
    write(json_unit, json_int_fmt) 'n_det', N_det
    if (s2_eig) then
      write(json_unit, json_int_fmt) 'n_cfg', N_configuration
      if (only_expected_s2) then
        write(json_unit, json_int_fmt) 'n_csf', N_csf
      endif
    endif
    write(json_unit, json_array_open_fmt) 'states'
    do k=1,N_states_p
      write(json_unit, json_dict_uopen_fmt)
      write(json_unit, json_real_fmt) 'energy', psi_energy_with_nucl_rep(k)
      write(json_unit, json_real_fmt) 's2', psi_s2(k)
      write(json_unit, json_real_fmt) 'pt2', pt2_data % pt2(k)
      write(json_unit, json_real_fmt) 'pt2_err', pt2_data_err % pt2(k)
      write(json_unit, json_real_fmt) 'rpt2', pt2_data % rpt2(k)
      write(json_unit, json_real_fmt) 'rpt2_err', pt2_data_err % rpt2(k)
      write(json_unit, json_real_fmt) 'variance', pt2_data % variance(k)
      write(json_unit, json_real_fmt) 'variance_err', pt2_data_err % variance(k)
      write(json_unit, json_array_open_fmt) 'ex_energy'
      do i=2,N_iter_p
          write(json_unit, fmtk(i)) extrapolated_energy(i,k)
      enddo
      write(json_unit, json_array_close_fmtx)
      if (k < N_states_p) then
        write(json_unit, json_dict_close_fmt)
      else
        write(json_unit, json_dict_close_fmtx)
      endif
    enddo
    write(json_unit, json_array_close_fmtx)
    write(json_unit, json_dict_close_fmt)
    deallocate(fmtk)
    call unlock_io
end
