module slater_rules_opt
  use iso_c_binding
  implicit none

  interface
    subroutine get_excitation_degree(nint, key1, key2, degree) bind(C, name="get_excitation_degree_c")
      import :: c_int, c_int64_t
      integer(c_int), value, intent(in)  :: nint
      integer(c_int64_t), intent(in)     :: key1(*)
      integer(c_int64_t), intent(in)     :: key2(*)
      integer(c_int), intent(out)        :: degree
    end subroutine get_excitation_degree
  end interface

end module slater_rules_opt
