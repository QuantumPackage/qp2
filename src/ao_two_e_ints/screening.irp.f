logical function ao_two_e_integral_zero(i,j,k,l)
  BEGIN_DOC
  ! Returns .True. if the AO integral $(ik|jl)$ can be taken as
  ! zero: the 1st+3rd arguments form the first electron pair and
  ! the 2nd+4th arguments the second one, i.e. the integral
  ! $(ik|jl)$ of ao_two_e_integral(i,k,j,l). No screening when the
  ! integrals are read from the EZFIO file, when $is\_periodic$, or
  ! for complex AOs ($use\_cgtos$). Otherwise, zero if
  ! $ao\_overlap\_abs(i,k)\,ao\_overlap\_abs(j,l)$ or the
  ! Cauchy-Schwarz bound $\sqrt{(ik|ik)\,(jl|jl)}$ (with
  ! $ao\_two\_e\_integral\_schwartz(p,q)=\sqrt{(pq|pq)}$) is below
  ! $ao\_integrals\_threshold$.
  END_DOC
  implicit none
  integer, intent(in) :: i,j,k,l

  ao_two_e_integral_zero = .False.
  if (.not.(read_ao_two_e_integrals.or.is_periodic.or.use_cgtos)) then
    if (ao_overlap_abs(j,l)*ao_overlap_abs(i,k)  < ao_integrals_threshold) then
        ao_two_e_integral_zero = .True.
        return
    endif
    if (ao_two_e_integral_schwartz(j,l)*ao_two_e_integral_schwartz(i,k)  < ao_integrals_threshold) then
        ao_two_e_integral_zero = .True.
    endif
  endif
end
