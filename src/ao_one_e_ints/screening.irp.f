logical function ao_one_e_integral_zero(i,k)
  implicit none
  integer, intent(in) :: i,k

  ao_one_e_integral_zero = .False.
  if (.not.(read_ao_one_e_integrals.or.is_periodic)) then
    if (ao_overlap_abs(i,k) < ao_integrals_threshold) then
        ao_one_e_integral_zero = .True.
        return
    endif
  endif
  if (ao_two_e_integral_schwartz(i,k) < ao_integrals_threshold) then
      ao_one_e_integral_zero = .True.
  endif
end

