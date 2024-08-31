logical function ao_two_e_integral_zero(i,j,k,l)
  implicit none
  integer, intent(in) :: i,j,k,l

  ao_two_e_integral_zero = .False.
  if(do_torus) return

  if (.not.(read_ao_two_e_integrals.or.is_periodic)) then
    if (ao_overlap_abs(j,l)*ao_overlap_abs(i,k)  < ao_integrals_threshold) then
        ao_two_e_integral_zero = .True.
        return
    endif
    if (ao_two_e_integral_schwartz(j,l)*ao_two_e_integral_schwartz(i,k)  < ao_integrals_threshold) then
        ao_two_e_integral_zero = .True.
    endif
  endif
end
