use iso_c_binding
use libcint_interface

!> Fill a matrix with one-electron integrals computed by libcint.
!>
!> `op` is one of `LIB_OP_OVLP`, `LIB_OP_KIN`, `LIB_OP_NUC`. The matrix is
!> symmetric and both triangles are written, so that the caller does not
!> have to symmetrize.
!>
!> The block of a shell pair is corrected by `libcint_scale` on both
!> indices, which carries the per-|AO| normalization that libcint cannot
!> represent (see :f:func:`libcint_scale`).
!>
!> The operator is selected once, before any loop, and each branch owns its
!> copy of the shell-pair loop: no branch on `op` is executed inside the
!> loops.
subroutine libcint_ao_one_e(op, mat)

  use iso_c_binding
  use libcint_interface
  implicit none

  BEGIN_DOC
  ! One-electron |AO| integrals with libcint: overlap (`LIB_OP_OVLP`),
  ! kinetic energy (`LIB_OP_KIN`) or nucleus-electron attraction
  ! (`LIB_OP_NUC`).
  END_DOC

  integer, intent(in)           :: op
  double precision, intent(out) :: mat(ao_num, ao_num)

  integer                       :: si, sj, ni, nj, ip, jp, i, j
  integer                       :: shls(2)
  integer(c_int)                :: info
  double precision, allocatable :: blk(:)

  PROVIDE ao_first_of_shell shell_ang_mom ao_l_max shell_num
  PROVIDE libcint_atm libcint_bas libcint_env libcint_scale
  PROVIDE libcint_natm libcint_nbas

  select case (op)

  case (LIB_OP_OVLP)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(si,sj,ni,nj,ip,jp,i,j,shls,blk,info)
    allocate(blk(libcint_ncart(ao_l_max)**2))
    !$OMP DO SCHEDULE(DYNAMIC)
    do sj = 1, shell_num
      nj = libcint_ncart(shell_ang_mom(sj))
      do si = 1, sj
        ni = libcint_ncart(shell_ang_mom(si))
        shls(1) = si - 1
        shls(2) = sj - 1
        ! the "block not empty" flag returned by libcint is unused
        info = cint1e_ovlp_cart(blk, shls, libcint_atm, libcint_natm, &
            libcint_bas, libcint_nbas, libcint_env, c_null_ptr)
        do jp = 1, nj
          j = ao_first_of_shell(sj) + jp - 1
          do ip = 1, ni
            i = ao_first_of_shell(si) + ip - 1
            mat(i,j) = blk(ip + (jp-1)*ni) &
                * libcint_scale(i) * libcint_scale(j)
            if (i /= j) then
              mat(j,i) = mat(i,j)
            endif
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    deallocate(blk)
    !$OMP END PARALLEL

  case (LIB_OP_KIN)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(si,sj,ni,nj,ip,jp,i,j,shls,blk,info)
    allocate(blk(libcint_ncart(ao_l_max)**2))
    !$OMP DO SCHEDULE(DYNAMIC)
    do sj = 1, shell_num
      nj = libcint_ncart(shell_ang_mom(sj))
      do si = 1, sj
        ni = libcint_ncart(shell_ang_mom(si))
        shls(1) = si - 1
        shls(2) = sj - 1
        ! the "block not empty" flag returned by libcint is unused
        info = cint1e_kin_cart(blk, shls, libcint_atm, libcint_natm, &
            libcint_bas, libcint_nbas, libcint_env, c_null_ptr)
        do jp = 1, nj
          j = ao_first_of_shell(sj) + jp - 1
          do ip = 1, ni
            i = ao_first_of_shell(si) + ip - 1
            mat(i,j) = blk(ip + (jp-1)*ni) &
                * libcint_scale(i) * libcint_scale(j)
            if (i /= j) then
              mat(j,i) = mat(i,j)
            endif
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    deallocate(blk)
    !$OMP END PARALLEL

  case (LIB_OP_NUC)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(si,sj,ni,nj,ip,jp,i,j,shls,blk,info)
    allocate(blk(libcint_ncart(ao_l_max)**2))
    !$OMP DO SCHEDULE(DYNAMIC)
    do sj = 1, shell_num
      nj = libcint_ncart(shell_ang_mom(sj))
      do si = 1, sj
        ni = libcint_ncart(shell_ang_mom(si))
        shls(1) = si - 1
        shls(2) = sj - 1
        ! the "block not empty" flag returned by libcint is unused
        info = cint1e_nuc_cart(blk, shls, libcint_atm, libcint_natm, &
            libcint_bas, libcint_nbas, libcint_env, c_null_ptr)
        do jp = 1, nj
          j = ao_first_of_shell(sj) + jp - 1
          do ip = 1, ni
            i = ao_first_of_shell(si) + ip - 1
            mat(i,j) = blk(ip + (jp-1)*ni) &
                * libcint_scale(i) * libcint_scale(j)
            if (i /= j) then
              mat(j,i) = mat(i,j)
            endif
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    deallocate(blk)
    !$OMP END PARALLEL

  case default

    stop 'libcint_ao_one_e: unknown operator'

  end select

end subroutine libcint_ao_one_e

!> Overlap matrix in the |AO| basis.
subroutine libcint_ao_overlap(mat)

  use libcint_interface, only: LIB_OP_OVLP
  implicit none

  BEGIN_DOC
  ! Overlap matrix in the |AO| basis, from libcint.
  END_DOC

  double precision, intent(out) :: mat(ao_num, ao_num)

  call libcint_ao_one_e(LIB_OP_OVLP, mat)

end subroutine libcint_ao_overlap

!> Kinetic energy matrix in the |AO| basis.
subroutine libcint_ao_kinetic(mat)

  use libcint_interface, only: LIB_OP_KIN
  implicit none

  BEGIN_DOC
  ! Kinetic energy matrix in the |AO| basis, from libcint.
  END_DOC

  double precision, intent(out) :: mat(ao_num, ao_num)

  call libcint_ao_one_e(LIB_OP_KIN, mat)

end subroutine libcint_ao_kinetic

!> Nucleus-electron attraction matrix in the |AO| basis.
!>
!> Only valid when every nuclear charge is an integer: see
!> :f:func:`libcint_nuc_available`.
subroutine libcint_ao_nuc(mat)

  use libcint_interface, only: LIB_OP_NUC
  implicit none

  BEGIN_DOC
  ! Nucleus-electron attraction matrix in the |AO| basis, from libcint.
  END_DOC

  double precision, intent(out) :: mat(ao_num, ao_num)

  call libcint_ao_one_e(LIB_OP_NUC, mat)

end subroutine libcint_ao_nuc
