use sizet
logical function ao_one_e_integral_zero(i,k)
 use sizet
  BEGIN_DOC
  ! Returns .True. if the AO integrals between the |AO|s $i$ and
  ! $k$ can be taken as zero. No screening when
  ! $io\_ao\_integrals\_overlap$ is not 'None' (overlap integrals
  ! read or written to the EZFIO file), when $is\_periodic$, or
  ! for complex AOs ($use\_cgtos$). Otherwise, zero if
  ! $ao\_overlap\_abs(i,k) < ao\_one\_e\_integrals\_threshold$.
  END_DOC
  implicit none
  integer(size_t), intent(in) :: i,k

  ao_one_e_integral_zero = .False.
  if (.not.((io_ao_integrals_overlap/='None').or.is_periodic.or.use_cgtos)) then
    if (ao_overlap_abs(i,k) < ao_one_e_integrals_threshold) then
        ao_one_e_integral_zero = .True.
        return
    endif
  endif
end

