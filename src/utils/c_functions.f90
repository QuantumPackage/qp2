module c_functions
  use iso_c_binding

  interface
     subroutine usleep_c(us) bind (C,name="usleep_c")
       use iso_c_binding
       integer(c_int), value :: us
     end subroutine usleep_c
  end interface

  interface
     integer(c_int) function atoi_c(a) bind (C,name="atoi")
       use iso_c_binding
       character(kind=c_char), intent(in) :: a(*)
     end function atoi_c
  end interface

  interface
     subroutine sscanf_ss_c(str,s1, s2) bind (C)
       use iso_c_binding
       character(kind=c_char), intent(in ) :: str(*)
       character(kind=c_char), intent(out) :: s1(*),s2(*)
     end subroutine sscanf_ss_c
  end interface

  interface
     subroutine sscanf_ssds_c(str, s1, s2, i, s3) bind (C)
       use iso_c_binding
       character(kind=c_char), intent(in ) :: str(*)
       character(kind=c_char), intent(out) :: s1(*),s2(*),s3(*)
       integer(kind=c_int)   , intent(out) :: i
     end subroutine sscanf_ssds_c
  end interface

  interface
     subroutine sscanf_dd_c(str, i1, i2) bind (C)
       use iso_c_binding
       character(kind=c_char), intent(in ) :: str(*)
       integer(kind=c_int)   , intent(out) :: i1, i2
     end subroutine sscanf_dd_c
  end interface

  interface
     subroutine sscanf_ddd_c(str, i1, i2, i3) bind (C)
       use iso_c_binding
       character(kind=c_char), intent(in ) :: str(*)
       integer(kind=c_int)   , intent(out) :: i1, i2, i3
     end subroutine sscanf_ddd_c
  end interface

  interface
     subroutine sscanf_sd_c(str,s1, i) bind (C)
       use iso_c_binding
       character(kind=c_char), intent(in ) :: str(*)
       character(kind=c_char), intent(out) :: s1(*)
       integer(kind=c_int)   , intent(out) :: i
     end subroutine sscanf_sd_c
  end interface

  interface
    integer(kind=c_int) function mkl_serv_intel_cpu_true() bind(C)
      use iso_c_binding
    end function
  end interface

contains

  integer function atoi(a)
    implicit none
    character(len=*), intent(in) :: a
    atoi = atoi_c(trim(a)//c_null_char)
  end function atoi

end module c_functions

subroutine sscanf_ss(str, s1,s2)
  use c_functions
  use iso_c_binding
  implicit none
  character(*), intent(in)  :: str
  character(*), intent(out) :: s1,s2
  s1 = ' '
  s2 = ' '
  call sscanf_ss_c(trim(str)//c_null_char, s1, s2)
end subroutine sscanf_ss

subroutine sscanf_sd(str, s1,i)
  use c_functions
  use iso_c_binding
  implicit none
  character(*), intent(in)  :: str
  character(*), intent(out) :: s1
  integer, intent(out)      :: i
  s1 = ' '
  call sscanf_sd_c(trim(str)//c_null_char, s1, i)
end subroutine sscanf_sd

subroutine sscanf_ssds(str, s1,s2,i,s3)
  use c_functions
  use iso_c_binding
  implicit none
  character(*), intent(in)  :: str
  character(*), intent(out) :: s1,s2,s3
  integer, intent(out)      :: i
  s1 = ' '
  s2 = ' '
  s3 = ' '
  call sscanf_ssds_c(trim(str)//c_null_char, s1, s2, i, s3)
end subroutine sscanf_ssds

subroutine sscanf_dd(str, i1,i2)
  use c_functions
  use iso_c_binding
  implicit none
  character(*), intent(in)  :: str
  integer, intent(out)      :: i1, i2
  call sscanf_dd_c(trim(str)//c_null_char, i1, i2)
end subroutine sscanf_dd

subroutine sscanf_ddd(str, i1,i2,i3)
  use c_functions
  use iso_c_binding
  implicit none
  character(*), intent(in)  :: str
  integer, intent(out)      :: i1, i2, i3
  call sscanf_ddd_c(trim(str)//c_null_char, i1, i2, i3)
end subroutine sscanf_ddd


subroutine usleep(us)
  use iso_c_binding
  use c_functions
  implicit none
  integer, intent(in) :: us
  integer(c_int) :: u
  u = us
  call usleep_c(u)
end subroutine usleep

