
! ---

subroutine LTxSxR(n, m, L, S, R, C)

  implicit none
  integer,          intent(in)  :: n, m
  double precision, intent(in)  :: L(n,m), S(n,n), R(n,m)
  double precision, intent(out) :: C(m,m)
  integer                       :: i, j
  double precision              :: accu_d, accu_nd
  double precision, allocatable :: tmp(:,:)

  ! L.T x S x R
  allocate(tmp(m,n))
  call dgemm( 'T', 'N', m, n, n, 1.d0      &
            , L, size(L, 1), S, size(S, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  call dgemm( 'N', 'N', m, m, n, 1.d0          &
            , tmp, size(tmp, 1), R, size(R, 1) &
            , 0.d0, C, size(C, 1) )
  deallocate(tmp)

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(j.eq.i) then
        accu_d += dabs(C(j,i))
      else
        accu_nd += C(j,i) * C(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  print*, ' accu_d  = ', accu_d
  print*, ' accu_nd = ', accu_nd

end subroutine LTxR

! ---

