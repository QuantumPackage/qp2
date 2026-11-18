logical function ao_two_e_integral_zero(i,j,k,l)
  BEGIN_DOC
  ! Returns .True. if the AO integral $(ij|kl)$ can be taken as
  ! zero. No screening when the integrals are read from the EZFIO
  ! file or when $is\_periodic$.
  ! Otherwise, zero if $ao\_overlap\_abs(i,k)\,ao\_overlap\_abs(j,l)$
  ! or the Schwartz bound $\sqrt{(ii|kk)\,(jj|ll)}$ is below
  ! $ao\_integrals\_threshold$.
  ! For cGTOs ($use\_cgtos$) the overlap test is not used (the overlap of
  ! two oscillating functions can vanish while their integrals do not);
  ! the test is the Cauchy-Schwarz bound $\sqrt{(ik|ik)\,(jl|jl)}$,
  ! rigorous for the real cGTO AOs ($ao\_two\_e\_integral\_schwartz\_cgtos$).
  END_DOC
  implicit none
  integer, intent(in) :: i,j,k,l

  ao_two_e_integral_zero = .False.
  if (read_ao_two_e_integrals.or.is_periodic) return

  if (use_cgtos) then
    ! Cauchy-Schwarz bound |(ik|jl)| <= sqrt( (ik|ik) (jl|jl) )
    if (ao_two_e_integral_schwartz_cgtos(j,l) *                     &
        ao_two_e_integral_schwartz_cgtos(i,k) < ao_integrals_threshold) then
      ao_two_e_integral_zero = .True.
    endif
  else
    if (ao_overlap_abs(j,l)*ao_overlap_abs(i,k)  < ao_integrals_threshold) then
        ao_two_e_integral_zero = .True.
        return
    endif
    if (ao_two_e_integral_schwartz(j,l)*ao_two_e_integral_schwartz(i,k)  < ao_integrals_threshold) then
        ao_two_e_integral_zero = .True.
    endif
  endif
end
