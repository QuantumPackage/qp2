logical function ao_one_e_integral_zero(i,k)
  implicit none
  integer, intent(in) :: i,k

  ao_one_e_integral_zero = .False.
  if((io_ao_integrals_overlap=='None').and.(.not.use_cgtos)) then
    if (ao_overlap_abs(i,k) < ao_one_e_integrals_threshold) then
        ao_one_e_integral_zero = .True.
        return
    endif
  endif
end

