
! ---

subroutine jBH_elem_fct_grad(alpha, r1, r2, fct, grad1_fct)

  implicit none
  double precision, intent(in)  :: alpha, r1(3), r2(3)
  double precision, intent(out) :: fct, grad1_fct(3)
  double precision              :: dist, tmp1, tmp2

  dist = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
              + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
              + (r1(3) - r2(3)) * (r1(3) - r2(3)) )


  if(dist .ge. 1d-10) then
    tmp1 = 1.d0 / (1.d0 + alpha * dist)

    fct = alpha * dist * tmp1
    tmp2 = alpha * tmp1 * tmp1 / dist
    grad1_fct(1) = tmp2 * (r1(1) - r2(1))
    grad1_fct(2) = tmp2 * (r1(2) - r2(2))
    grad1_fct(3) = tmp2 * (r1(3) - r2(3))
  else
    grad1_fct(1) = 0.d0
    grad1_fct(2) = 0.d0
    grad1_fct(3) = 0.d0
    fct = 0.d0
  endif

  return
end 

! ---

