module c_functions
  use iso_c_binding

  interface
    subroutine usleep_c(us) bind (C,name="usleep")
      use iso_c_binding
      integer(c_int), value :: us
    end subroutine usleep_c
  end interface

end module

subroutine usleep(us)
  use c_functions
  use iso_c_binding
  implicit none
  integer, intent(in) :: us
  call usleep_c(int(us,c_int))
end
