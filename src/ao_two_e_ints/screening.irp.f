logical function ao_two_e_integral_zero(i,j,k,l)
  BEGIN_DOC
  ! Returns .True. if the AO integral $(ik|jl)$ can be taken as
  ! zero: the 1st and 3rd arguments form the pair on electron 1,
  ! the 2nd and 4th the pair on electron 2, i.e. the integral is
  ! $ao\_two\_e\_integral(i,k,j,l)$. No screening when the
  ! integrals are read from the EZFIO file or when $is\_periodic$.
  ! For real GTOs the test is the rigorous Cauchy-Schwarz bound
  ! $\sqrt{(ik|ik)\,(jl|jl)}$, with
  ! $ao\_two\_e\_integral\_schwartz(p,q) = \sqrt{(pq|pq)}$ and the
  ! diagonals summed over all primitive quartets; for cGTOs
  ! ($use\_cgtos$) the same bound via
  ! $ao\_two\_e\_integral\_schwartz\_cgtos$ (the overlap of two
  ! oscillating functions can vanish while their integrals do
  ! not). The bound is rigorous, so no integral above
  ! $ao\_integrals\_threshold$ is screened.
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
    ! The former overlap test ao_overlap_abs(j,l)*ao_overlap_abs(i,k)
    ! was removed (2026-11, issue #419): the product of overlaps is
    ! not a bound on (ik|jl) and dropped integrals up to ~10x the
    ! threshold (S_ik*S_jl ~ 1e-15 with (ik|jl) = 1.1e-15); the CS
    ! bound above covers it at zero cost.
    if (ao_two_e_integral_schwartz(j,l)*ao_two_e_integral_schwartz(i,k)  < ao_integrals_threshold) then
        ao_two_e_integral_zero = .True.
    endif
  endif
end
