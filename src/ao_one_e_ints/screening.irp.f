logical function ao_one_e_integral_zero(i,k)
  implicit none
  integer, intent(in) :: i,k

  ao_one_e_integral_zero = .False.
  if (.not.((io_ao_integrals_overlap/='None').or.is_complex)) then
    if (ao_overlap_abs(i,k) < ao_integrals_threshold) then
        ao_one_e_integral_zero = .True.
        return
    endif
  endif
end

