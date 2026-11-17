logical function ao_two_e_integral_zero(i,j,k,l)
  BEGIN_DOC
  ! Returns .True. if the AO integral $(ij|kl)$ can be taken as
  ! zero. No screening when the integrals are read from the EZFIO
  ! file or when $is\_periodic$.
  ! Otherwise, zero if $ao\_overlap\_abs(i,k)\,ao\_overlap\_abs(j,l)$
  ! or the Schwartz bound $\sqrt{(ii|kk)\,(jj|ll)}$ is below
  ! $ao\_integrals\_threshold$.
  ! For cGTOs ($use\_cgtos$) only the Schwartz bound is used: the AOs are
  ! real, so the bound is rigorous, but the overlap of two oscillating
  ! functions can vanish while their integrals do not.
  END_DOC
  implicit none
  integer, intent(in) :: i,j,k,l

  ao_two_e_integral_zero = .False.
  if (.not.(read_ao_two_e_integrals.or.is_periodic)) then
    if (.not.use_cgtos) then
      if (ao_overlap_abs(j,l)*ao_overlap_abs(i,k)  < ao_integrals_threshold) then
        ao_two_e_integral_zero = .True.
        return
      endif
    endif
    if (ao_two_e_integral_schwartz(j,l)*ao_two_e_integral_schwartz(i,k)  < ao_integrals_threshold) then
        ao_two_e_integral_zero = .True.
    endif
  endif
end
