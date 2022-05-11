subroutine format_w_error(value,error,size_nb,max_nb_digits,format_value,str_error)

  implicit none
 
  BEGIN_DOC
  ! Format for double precision, value(error)
  END_DOC

  ! in
  ! | value         | double precision | value...      |
  ! | error         | double precision | error...      |
  ! | size_nb       | integer          | X in FX.Y     |
  ! | max_nb_digits | integer          | Max Y in FX.Y |

  ! out
  ! | format_value  | character | string FX.Y for the format |
  ! | str_error     | character | string of the error |

  ! internal
  ! | str_size      | character | size in string format                             |
  ! | nb_digits     | integer   | number of digits Y in FX.Y depending of the error |
  ! | str_nb_digits | character | nb_digits in string format                        |
  ! | str_exp       | character | string of the value in exponential format         |

  ! in
  double precision, intent(in)   :: error, value
  integer, intent(in)            :: size_nb, max_nb_digits

  ! out
  character(len=20), intent(out) :: str_error, format_value

  ! internal
  character(len=20)              :: str_size, str_nb_digits, str_exp
  integer                        :: nb_digits

  ! max_nb_digit: Y max
  ! size_nb = Size of the double: X (FX.Y)
  write(str_size,'(I3)') size_nb

  ! Error
  write(str_exp,'(1pE20.0)') error
  str_error = trim(adjustl(str_exp))
  
  ! Number of digit: Y (FX.Y) from the exponent
  str_nb_digits = str_exp(19:20)
  read(str_nb_digits,*) nb_digits
 
  ! If the error is 0d0
  if (error <= 1d-16) then 
    write(str_nb_digits,*) max_nb_digits
  endif

  ! If the error is too small 
  if (nb_digits > max_nb_digits) then
      write(str_nb_digits,*) max_nb_digits
      str_error(1:1) = '0'
  endif

  ! If the error is too big (>= 0.5)
  if (error >= 0.5d0) then
    str_nb_digits = '1'
    str_error(1:1) = '*'
  endif

  ! FX.Y,A1,A1,A1 for value(str_error)
  !string = 'F'//trim(adjustl(str_size))//'.'//trim(adjustl(str_nb_digits))//',A1,A1,A1'

  ! FX.Y just for the value 
  format_value = 'F'//trim(adjustl(str_size))//'.'//trim(adjustl(str_nb_digits))

end
