module slater_rules_opt
  use iso_c_binding
  implicit none

  interface
    subroutine get_excitation_degree_c(key1, key2, degree, nint) bind(C)
      import :: c_int, c_int64_t
      integer(c_int), value, intent(in)  :: nint
      integer(c_int64_t), intent(in)     :: key1(2*nint)
      integer(c_int64_t), intent(in)     :: key2(2*nint)
      integer(c_int), intent(out)        :: degree
    end subroutine get_excitation_degree_c
  end interface

end module slater_rules_opt
