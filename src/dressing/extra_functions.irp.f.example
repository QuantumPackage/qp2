 BEGIN_PROVIDER [ integer, N_dress_int_buffer ]
&BEGIN_PROVIDER [ integer, N_dress_double_buffer ]
&BEGIN_PROVIDER [ integer, N_dress_det_buffer ]
  implicit none
  N_dress_int_buffer = 1
  N_dress_double_buffer = 1
  N_dress_det_buffer = 1
END_PROVIDER



subroutine delta_ij_done()
  BEGIN_DOC
  ! This subroutine is executed on the master when the dressing has been computed,
  ! before the diagonalization.
  END_DOC
end

subroutine dress_pulled(ind, int_buf, double_buf, det_buf, N_buf)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Dress the contributions pulled from the slave.
  END_DOC

  integer, intent(in) :: ind, N_buf(3)
  integer, intent(in) :: int_buf(*)
  double precision, intent(in) :: double_buf(*)
  integer(bit_kind), intent(in) :: det_buf(N_int,2,*)
end

subroutine generator_start(i_gen, iproc)
  implicit none
  BEGIN_DOC
  ! This subroutine is executed on the slave before computing the contribution of a generator.
  END_DOC

  integer, intent(in) :: i_gen, iproc
  integer :: i
end

subroutine generator_done(i_gen, int_buf, double_buf, det_buf, N_buf, iproc)
  implicit none
  BEGIN_DOC
  ! This subroutine is executed on the slave after computing the contribution of a generator.
  END_DOC
  integer, intent(in) :: i_gen, iproc
  integer, intent(out) :: int_buf(N_dress_int_buffer), N_buf(3)
  double precision, intent(out) :: double_buf(N_dress_double_buffer)
  integer(bit_kind), intent(out) :: det_buf(N_int, 2, N_dress_det_buffer)
  N_buf(:) = 1
  int_buf(:) = 0
  double_buf(:) = 0.d0
  det_buf(:,:,:) = 0
end



