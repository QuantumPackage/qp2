
! ---

subroutine get_grad1_u12_r1_2e(r1, n_grid2, gradx, grady, gradz)

  BEGIN_DOC
  !
  !  d/dx1 j_2e(1,2)
  !  d/dy1 j_2e(1,2)
  !  d/dz1 j_2e(1,2)
  !
  END_DOC

  include 'constants.include.F'

  implicit none
  integer         , intent(in)  :: n_grid2
  double precision, intent(in)  :: r1(3)
  double precision, intent(out) :: gradx(n_grid2)
  double precision, intent(out) :: grady(n_grid2)
  double precision, intent(out) :: gradz(n_grid2)

  integer                       :: jpoint
  integer                       :: i_nucl, p, mpA, npA, opA
  integer                       :: powmax1, powmax, powmax2
  double precision              :: r2(3)
  double precision              :: tmp, tmp1, tmp2
  double precision              :: rn(3), f1A, grad1_f1A(3), f2A, grad2_f2A(3), g12, grad1_g12(3)
  double precision, allocatable :: f1A_power(:), f2A_power(:), double_p(:), g12_power(:)


  powmax1 = max(maxval(jBH_m), maxval(jBH_n))
  powmax2 = maxval(jBH_o)
  powmax  = max(powmax1, powmax2)

  allocate(f1A_power(-1:powmax), f2A_power(-1:powmax), g12_power(-1:powmax), double_p(0:powmax))

  do p = 0, powmax
    double_p(p) = dble(p)
  enddo

  f1A_power(-1) = 0.d0
  f2A_power(-1) = 0.d0
  g12_power(-1) = 0.d0

  f1A_power(0) = 1.d0
  f2A_power(0) = 1.d0
  g12_power(0) = 1.d0

  do jpoint = 1, n_points_extra_final_grid ! r2

    r2(1) = final_grid_points_extra(1,jpoint)
    r2(2) = final_grid_points_extra(2,jpoint)
    r2(3) = final_grid_points_extra(3,jpoint)

    gradx(jpoint) = 0.d0 
    grady(jpoint) = 0.d0 
    gradz(jpoint) = 0.d0 
    do i_nucl = 1, nucl_num 

      rn(1) = nucl_coord(i_nucl,1)
      rn(2) = nucl_coord(i_nucl,2)
      rn(3) = nucl_coord(i_nucl,3)

      call jBH_elem_fct_grad(jBH_en(i_nucl), r1, rn, f1A, grad1_f1A)
      call jBH_elem_fct_grad(jBH_en(i_nucl), r2, rn, f2A, grad2_f2A)
      call jBH_elem_fct_grad(jBH_ee(i_nucl), r1, r2, g12, grad1_g12)

      ! Compute powers of f1A and f2A
      do p = 1, powmax1
        f1A_power(p) = f1A_power(p-1) * f1A
        f2A_power(p) = f2A_power(p-1) * f2A
      enddo
      do p = 1, powmax2
        g12_power(p) = g12_power(p-1) * g12
      enddo

      do p = 1, jBH_size
        mpA = jBH_m(p,i_nucl)
        npA = jBH_n(p,i_nucl)
        opA = jBH_o(p,i_nucl)
        tmp = jBH_c(p,i_nucl)
        if(mpA .eq. npA) then
          tmp = tmp * 0.5d0
        endif

        tmp1 = double_p(mpA) * f1A_power(mpA-1) * f2A_power(npA) + double_p(npA) * f1A_power(npA-1) * f2A_power(mpA)
        tmp1 = tmp1 * g12_power(opA) * tmp
        tmp2 = double_p(opA) * g12_power(opA-1) * (f1A_power(mpA) * f2A_power(npA) + f1A_power(npA) * f2A_power(mpA)) * tmp

        gradx(jpoint) = gradx(jpoint) + tmp1 * grad1_f1A(1) + tmp2 * grad1_g12(1)
        grady(jpoint) = grady(jpoint) + tmp1 * grad1_f1A(2) + tmp2 * grad1_g12(2)
        gradz(jpoint) = gradz(jpoint) + tmp1 * grad1_f1A(3) + tmp2 * grad1_g12(3)
      enddo ! p
    enddo ! i_nucl
  enddo ! jpoint

  return
end

! ---

