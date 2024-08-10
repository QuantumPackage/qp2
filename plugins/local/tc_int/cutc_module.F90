
module cutc_module

  use, intrinsic :: iso_c_binding

  implicit none

  interface

    ! ---

    subroutine cutc_int_c(nxBlocks, nyBlocks, nzBlocks,           &
                          blockxSize, blockySize, blockzSize,     &
                          n_grid1, n_grid2, n_ao, n_nuc, size_bh, &
                          r1, wr1, r2, wr2, rn,                   &
                          aos_data1, aos_data2,                   &
                          c_bh, m_bh, n_bh, o_bh,                 &
                          int2_grad1_u12_ao, int_2e_ao) bind(C, name = "cutc_int_c")

      import c_int, c_double, c_ptr
      integer(c_int), intent(in), value :: nxBlocks, blockxSize
      integer(c_int), intent(in), value :: nyBlocks, blockySize
      integer(c_int), intent(in), value :: nzBlocks, blockzSize
      integer(c_int), intent(in), value :: n_grid1, n_grid2
      integer(c_int), intent(in), value :: n_ao
      integer(c_int), intent(in), value :: n_nuc
      integer(c_int), intent(in), value :: size_bh
      real(c_double), intent(in)        :: r1(3,n_grid1), wr1(n_grid1)
      real(c_double), intent(in)        :: r2(3,n_grid2), wr2(n_grid2)
      real(c_double), intent(in)        :: rn(3,n_nuc)
      real(c_double), intent(in)        :: aos_data1(n_grid1,n_ao,4)
      real(c_double), intent(in)        :: aos_data2(n_grid2,n_ao,4)
      real(c_double), intent(in)        :: c_bh(size_bh,n_nuc)
      integer(c_int), intent(in)        :: m_bh(size_bh,n_nuc)
      integer(c_int), intent(in)        :: n_bh(size_bh,n_nuc)
      integer(c_int), intent(in)        :: o_bh(size_bh,n_nuc)
      real(c_double), intent(out)       :: int2_grad1_u12_ao(n_ao,n_ao,n_grid1,3)
      real(c_double), intent(out)       :: int_2e_ao(n_ao,n_ao,n_ao,n_ao)

    end subroutine cutc_int_c

    ! ---

    subroutine deb_int_2e_ao(nxBlocks, nyBlocks, nzBlocks,           &
                             blockxSize, blockySize, blockzSize,     &
                             n_grid1, n_grid2, n_ao, n_nuc, size_bh, &
                             r1, wr1, r2, wr2, rn,                   &
                             aos_data1, aos_data2,                   &
                             c_bh, m_bh, n_bh, o_bh,                 &
                             int2_grad1_u12_ao, int_2e_ao) bind(C, name = "deb_int_2e_ao")

      import c_int, c_double, c_ptr

      integer(c_int), intent(in), value :: nxBlocks, blockxSize
      integer(c_int), intent(in), value :: nyBlocks, blockySize
      integer(c_int), intent(in), value :: nzBlocks, blockzSize
      integer(c_int), intent(in), value :: n_grid1, n_grid2
      integer(c_int), intent(in), value :: n_ao
      integer(c_int), intent(in), value :: n_nuc
      integer(c_int), intent(in), value :: size_bh
      real(c_double), intent(in)        :: r1(3,n_grid1), wr1(n_grid1)
      real(c_double), intent(in)        :: r2(3,n_grid2), wr2(n_grid2)
      real(c_double), intent(in)        :: rn(3,n_nuc)
      real(c_double), intent(in)        :: aos_data1(n_grid1,n_ao,4)
      real(c_double), intent(in)        :: aos_data2(n_grid2,n_ao,4)
      real(c_double), intent(in)        :: c_bh(size_bh,n_nuc)
      integer(c_int), intent(in)        :: m_bh(size_bh,n_nuc)
      integer(c_int), intent(in)        :: n_bh(size_bh,n_nuc)
      integer(c_int), intent(in)        :: o_bh(size_bh,n_nuc)
      real(c_double), intent(out)       :: int2_grad1_u12_ao(n_ao,n_ao,n_grid1,3)
      real(c_double), intent(out)       :: int_2e_ao(n_ao,n_ao,n_ao,n_ao)

    end subroutine deb_int_2e_ao

    ! ---

    subroutine cutc_no_2e(n_grid1, n_mo, ne_a, ne_b,                   &
                          wr1, mos_l_in_r, mos_r_in_r, int2_grad1_u12, &
                          no_2e) bind(C, name = "cutc_no_2e")

      import c_int, c_double, c_ptr

      integer(c_int), intent(in), value :: n_grid1
      integer(c_int), intent(in), value :: n_mo
      integer(c_int), intent(in), value :: ne_a
      integer(c_int), intent(in), value :: ne_b
      real(c_double), intent(in)        :: wr1(n_grid1)
      real(c_double), intent(in)        :: mos_l_in_r(n_grid1,n_mo)
      real(c_double), intent(in)        :: mos_r_in_r(n_grid1,n_mo)
      real(c_double), intent(in)        :: int2_grad1_u12(n_grid1,3,n_mo,n_mo)
      real(c_double), intent(out)       :: no_2e(n_mo,n_mo,n_mo,n_mo)

    end subroutine cutc_no_2e

    ! ---

    subroutine deb_no_2e(n_grid1, n_mo, ne_a, ne_b,                   &
                         wr1, mos_l_in_r, mos_r_in_r, int2_grad1_u12, &
                         tmpO, tmpJ, tmpA, tmpB, tmpC, tmpD, tmpE,    &
                         no_2e) bind(C, name = "deb_no_2e")

      import c_int, c_double, c_ptr

      integer(c_int), intent(in), value :: n_grid1
      integer(c_int), intent(in), value :: n_mo
      integer(c_int), intent(in), value :: ne_a
      integer(c_int), intent(in), value :: ne_b
      real(c_double), intent(in)        :: wr1(n_grid1)
      real(c_double), intent(in)        :: mos_l_in_r(n_grid1,n_mo)
      real(c_double), intent(in)        :: mos_r_in_r(n_grid1,n_mo)
      real(c_double), intent(in)        :: int2_grad1_u12(n_grid1,3,n_mo,n_mo)
      real(c_double), intent(out)       :: tmpO(n_grid1), tmpJ(n_grid1,3)
      real(c_double), intent(out)       :: tmpA(n_grid1,3,n_mo), tmpB(n_grid1,3,n_mo)
      real(c_double), intent(out)       :: tmpC(n_grid1,4,n_mo,n_mo), tmpD(n_grid1,4,n_mo,n_mo)
      real(c_double), intent(out)       :: tmpE(n_mo,n_mo,n_mo,n_mo)
      real(c_double), intent(out)       :: no_2e(n_mo,n_mo,n_mo,n_mo)

    end subroutine deb_no_2e

    ! ---

    subroutine cutc_no_1e(n_grid1, n_mo, ne_a, ne_b,                   &
                          wr1, mos_l_in_r, mos_r_in_r, int2_grad1_u12, &
                          no_1e) bind(C, name = "cutc_no_1e")

      import c_int, c_double, c_ptr

      integer(c_int), intent(in), value :: n_grid1
      integer(c_int), intent(in), value :: n_mo
      integer(c_int), intent(in), value :: ne_a
      integer(c_int), intent(in), value :: ne_b
      real(c_double), intent(in)        :: wr1(n_grid1)
      real(c_double), intent(in)        :: mos_l_in_r(n_grid1,n_mo)
      real(c_double), intent(in)        :: mos_r_in_r(n_grid1,n_mo)
      real(c_double), intent(in)        :: int2_grad1_u12(n_grid1,3,n_mo,n_mo)
      real(c_double), intent(out)       :: no_1e(n_mo,n_mo)

    end subroutine cutc_no_1e

    ! ---

    subroutine deb_no_1e(n_grid1, n_mo, ne_a, ne_b,                                  &
                         wr1, mos_l_in_r, mos_r_in_r, int2_grad1_u12,                &
                         tmpO, tmpJ, tmpM, tmpS, tmpC, tmpD, tmpL, tmpR, tmpE, tmpF, &
                         no_1e) bind(C, name = "deb_no_1e")

      import c_int, c_double, c_ptr

      integer(c_int), intent(in), value :: n_grid1
      integer(c_int), intent(in), value :: n_mo
      integer(c_int), intent(in), value :: ne_a
      integer(c_int), intent(in), value :: ne_b
      real(c_double), intent(in)        :: wr1(n_grid1)
      real(c_double), intent(in)        :: mos_l_in_r(n_grid1,n_mo)
      real(c_double), intent(in)        :: mos_r_in_r(n_grid1,n_mo)
      real(c_double), intent(in)        :: int2_grad1_u12(n_grid1,3,n_mo,n_mo)
      real(c_double), intent(out)       :: tmpO(n_grid1)
      real(c_double), intent(out)       :: tmpJ(n_grid1,3)
      real(c_double), intent(out)       :: tmpM(n_grid1,3)
      real(c_double), intent(out)       :: tmpS(n_grid1)
      real(c_double), intent(out)       :: tmpC(n_grid1,4,n_mo,n_mo)
      real(c_double), intent(out)       :: tmpD(n_grid1,4)
      real(c_double), intent(out)       :: tmpL(n_grid1,3,n_mo)
      real(c_double), intent(out)       :: tmpR(n_grid1,3,n_mo)
      real(c_double), intent(out)       :: tmpE(n_grid1,5,n_mo)
      real(c_double), intent(out)       :: tmpF(n_grid1,5,n_mo)
      real(c_double), intent(out)       :: no_1e(n_mo,n_mo)

    end subroutine deb_no_1e

    ! ---

  end interface

end module cutc_module


