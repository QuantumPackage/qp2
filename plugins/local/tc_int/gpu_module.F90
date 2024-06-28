
! ---

module gpu_module

  use iso_c_binding

  implicit none

  interface

    subroutine tc_int_bh(n_grid1, n_grid2, ao_num, n_nuc, &
                         size_bh, m_bh, n_bh, o_bh, c_bh, &
                         r1, r2, rn, wr1, wr2, aos_data1, &
                         aos_data2, int2_grad1_u12, tc_int_2e_ao) bind(C)

      import c_int, c_double

      integer(c_int), intent(in), value :: n_grid1, n_grid2, ao_num, n_nuc, size_bh
      integer(c_int), intent(in)        :: m_bh(size_bh,n_nuc)
      integer(c_int), intent(in)        :: n_bh(size_bh,n_nuc)
      integer(c_int), intent(in)        :: o_bh(size_bh,n_nuc)
      real(c_double), intent(in)        :: c_bh(size_bh,n_nuc)
      real(c_double), intent(in)        :: r1(n_grid1,3), r2(n_grid2,3)
      real(c_double), intent(in)        :: rn(n_nuc,3)
      real(c_double), intent(in)        :: wr1(n_grid1), wr2(n_grid2)
      real(c_double), intent(in)        :: aos_data1(n_grid1,ao_num,4), aos_data2(n_grid2,ao_num,4)
      real(c_double), intent(out)       :: int2_grad1_u12(n_grid1,ao_num,ao_num,4)
      real(c_double), intent(out)       :: tc_int_2e_ao(ao_num,ao_num,ao_num,ao_num)

    end subroutine
    
  end interface

end module

! ---

