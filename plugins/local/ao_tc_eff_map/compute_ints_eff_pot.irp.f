

subroutine compute_ao_tc_sym_two_e_pot_jl(j, l, n_integrals, buffer_i, buffer_value)

  use map_module

  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC

  implicit none

  integer, intent(in)             :: j, l
  integer,intent(out)             :: n_integrals
  integer(key_kind),intent(out)   :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)

  integer                         :: i, k
  integer                         :: kk, m, j1, i1
  double precision                :: cpu_1, cpu_2, wall_1, wall_2
  double precision                :: integral, wall_0, integral_pot, integral_erf
  double precision                :: thr

  logical, external               :: ao_two_e_integral_zero
  double precision                :: ao_tc_sym_two_e_pot, ao_two_e_integral_erf
  double precision                :: env_gauss_2e_j1, env_gauss_2e_j2



  thr = ao_integrals_threshold

  n_integrals = 0

  j1 = j+ishft(l*l-l,-1)
  do k = 1, ao_num           ! r1
    i1 = ishft(k*k-k,-1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif

      if (ao_two_e_integral_erf_schwartz(i,k)*ao_two_e_integral_erf_schwartz(j,l) < thr) then
        cycle
      endif

      !DIR$ FORCEINLINE
      integral_pot = ao_tc_sym_two_e_pot  (i, k, j, l)  ! i,k : r1    j,l : r2
      integral_erf = ao_two_e_integral_erf(i, k, j, l)
      integral     = integral_erf + integral_pot

      if(abs(integral) < thr) then
        cycle
      endif

      n_integrals += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i, j, k, l, buffer_i(n_integrals))
      buffer_value(n_integrals) = integral
    enddo
  enddo

end subroutine compute_ao_tc_sym_two_e_pot_jl

