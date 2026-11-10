!> Two-electron |AO| integral :math:`\langle ij | kl \rangle` from libcint.
!>
!> The integral is taken from the block of the whole shell quartet, which is
!> computed once and kept in a one-entry per-thread cache
!> (:f:func:`libcint_two_e_block`): the callers walk the shell quartets in
!> an order that reuses a block `ncart(l)` times in a row, so a single entry
!> captures most of the locality.
!>
!> `lr` selects the long-range operator :math:`\mathrm{erf}(\mu r_{12})/r_{12}`
!> instead of :math:`1/r_{12}`; it is passed by the dispatcher so that this
!> module does not need to depend on the :file:`ao_two_e_ints` module. An
!> `lr` request with `mu_erf = 0` returns zero, since libcint reads a zero
!> `PTR_RANGE_OMEGA` as the plain Coulomb operator.
!>
!> The result is multiplied by :f:func:`libcint_scale` on the four indices,
!> which carries the per-|AO| normalization that libcint cannot represent.
double precision function libcint_ao_two_e_integral(i, j, k, l, lr)

  use libcint_block_cache, only: libcint_two_e_block
  implicit none

  BEGIN_DOC
  ! :math:`\langle ij | kl \rangle` with libcint, plain or long-range
  ! Coulomb operator.
  END_DOC

  integer, intent(in) :: i, j, k, l
  logical, intent(in) :: lr

  integer        :: si, sj, sk, sl, ni, nj, nk
  integer        :: dims(4)
  integer        :: ip, jp, kp, lp
  double precision, pointer :: blk(:)

  PROVIDE ao_shell ao_first_of_shell
  PROVIDE libcint_atm libcint_bas libcint_env libcint_env_lr
  PROVIDE libcint_scale libcint_natm libcint_nbas
  PROVIDE mu_erf

  ! libcint reads PTR_RANGE_OMEGA = 0 as the plain Coulomb operator, while
  ! the QP long-range kernel erf(mu*r12)/r12 is exactly zero at mu = 0.
  if (lr .and. mu_erf == 0.d0) then
    libcint_ao_two_e_integral = 0.d0
    return
  endif

  si = ao_shell(i)
  sj = ao_shell(j)
  sk = ao_shell(k)
  sl = ao_shell(l)

  call libcint_two_e_block(si, sj, sk, sl, lr, libcint_atm, libcint_natm, &
      libcint_bas, libcint_nbas, libcint_env, libcint_env_lr, blk, dims)

  ni = dims(1)
  nj = dims(2)
  nk = dims(3)

  ! libcint stores out[ip + jp*ni + kp*ni*nj + lp*ni*nj*nk], ip fastest
  ip = i - ao_first_of_shell(si)
  jp = (j - ao_first_of_shell(sj)) * ni
  kp = (k - ao_first_of_shell(sk)) * ni*nj
  lp = (l - ao_first_of_shell(sl)) * ni*nj*nk

  libcint_ao_two_e_integral = blk(1 + ip + jp + kp + lp) &
      * libcint_scale(i) * libcint_scale(j) &
      * libcint_scale(k) * libcint_scale(l)

end function libcint_ao_two_e_integral
