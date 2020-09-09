subroutine set_bit_to_integer(i_physical,key,Nint)
 use bitmasks
 BEGIN_DOC
 ! set to 1 the bit number i_physical in the bitstring key
 END_DOC
 implicit none
 integer, intent(in) :: i_physical,Nint
 integer(bit_kind), intent(inout) :: key(Nint)
 integer :: k,j,i
 k = ishft(i_physical-1,-bit_kind_shift)+1
 j = i_physical-ishft(k-1,bit_kind_shift)-1
 key(k) = ibset(key(k),j)
end


subroutine clear_bit_to_integer(i_physical,key,Nint)
 use bitmasks
 BEGIN_DOC
 ! set to 0 the bit number i_physical in the bitstring key
 END_DOC
 implicit none
 integer, intent(in) :: i_physical,Nint
 integer(bit_kind), intent(inout) :: key(Nint)
 integer :: k,j,i
 k = ishft(i_physical-1,-bit_kind_shift)+1
 j = i_physical-ishft(k-1,bit_kind_shift)-1
 key(k) = ibclr(key(k),j)
end



subroutine bitstring_to_list( string, list, n_elements, Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Gives the indices(+1) of the bits set to 1 in the bit string
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint)
  integer, intent(out)           :: list(Nint*bit_kind_size)
  integer, intent(out)           :: n_elements

  integer                        :: i, ishift
  integer(bit_kind)              :: l

  n_elements = 0
  ishift = 2
  do i=1,Nint
    l = string(i)
    do while (l /= 0_bit_kind)
      n_elements = n_elements+1
      list(n_elements) = ishift+popcnt(l-1_bit_kind) - popcnt(l)
      l = iand(l,l-1_bit_kind)
    enddo
    ishift = ishift + bit_kind_size
  enddo

end

subroutine list_to_bitstring( string, list, n_elements, Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the physical string "string(N_int,2)" from the array of
  ! occupations "list(N_int*bit_kind_size,2)
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(out) :: string(Nint)
  integer, intent(in)            :: list(Nint*bit_kind_size)
  integer, intent(in)            :: n_elements


  integer                        :: i, j
  integer                        :: ipos, iint

  !
  !                                       <== ipos ==>
  !                                                  |
  !                                                  v
  !string :|------------------------|-------------------------|------------------------|
  !        <==== bit_kind_size ====> <==== bit_kind_size ====> <==== bit_kind_size ====>
  !        {        iint            } {         iint         } {         iint         }
  !

  string = 0_bit_kind

  do i=1,n_elements
    iint = shiftr(list(i)-1,bit_kind_shift) + 1
    ipos = list(i)-shiftl((iint-1),bit_kind_shift)-1
    string(iint) = ibset( string(iint), ipos )
  enddo

end


subroutine bitstring_to_str( output, string, Nint )
  use bitmasks
  implicit none
  BEGIN_DOC
! Transform a bit string to a string for printing
  END_DOC
  character*(*), intent(out)     :: output
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint)

  integer                        :: i, j, ibuf
  integer(bit_kind)              :: itemp

  ibuf = 1
  output = ''
  output(ibuf:ibuf) = '|'
  ibuf = ibuf+1
  do i=1,Nint
    itemp = 1_bit_kind
    do j=1,bit_kind_size
      if (iand(itemp,string(i)) == itemp) then
        output(ibuf:ibuf) = '+'
      else
        output(ibuf:ibuf) = '-'
      endif
      ibuf = ibuf+1
      itemp = shiftl(itemp,1)
    enddo
  enddo
  output(ibuf:ibuf) = '|'
end


subroutine bitstring_to_hexa( output, string, Nint )
  use bitmasks
  implicit none
  BEGIN_DOC
! Transform a bit string to a string in hexadecimal format for printing
  END_DOC
  character*(*), intent(out)     :: output
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint)
  integer                        :: i, j, ibuf
  integer(bit_kind)              :: itemp
  character*(32)                 :: f

  write(f,*) '(Z',bit_kind_size/4,'.',bit_kind_size/4,')'
  ibuf = 1
  output = ''
  do i=Nint,1,-1
    write(output(ibuf:ibuf+bit_kind_size/4),f) string(i)
    ibuf = ibuf+bit_kind_size/4
  enddo
end

subroutine debug_det(string,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Subroutine to print the content of a determinant in '+-' notation and
  ! hexadecimal representation.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint,2)
  character*(2048)                :: output(2)
  call bitstring_to_hexa( output(1), string(1,1), Nint )
  call bitstring_to_hexa( output(2), string(1,2), Nint )
  print *,  trim(output(1)) , '|', trim(output(2))

  call print_det(string,Nint)

end

subroutine print_det(string,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Subroutine to print the content of a determinant using the '+-' notation
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint,2)
  character*(2048)                :: output(2)

  call bitstring_to_str( output(1), string(1,1), Nint )
  call bitstring_to_str( output(2), string(1,2), Nint )
  print *,  trim(output(1))
  print *,  trim(output(2))

end

subroutine debug_spindet(string,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Subroutine to print the content of a determinant in '+-' notation and
  ! hexadecimal representation.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint,2)
  character*(2048)                :: output(1)
  call bitstring_to_hexa( output(1), string(1,1), Nint )
  print *,  trim(output(1))
  call print_spindet(string,Nint)

end

subroutine print_spindet(string,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Subroutine to print the content of a determinant using the '+-' notation
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint,2)
  character*(2048)                :: output(1)

  call bitstring_to_str( output(1), string(1,1), Nint )
  print *,  trim(output(1))

end

subroutine debug_single_spindet(string,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Subroutine to print the content of a determinant in '+-' notation and
  ! hexadecimal representation.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint)
  character*(2048)                :: output(1)
  call bitstring_to_hexa( output(1), string(1), Nint )
  print *,  trim(output(1))
  call print_single_spindet(string,Nint)

end

subroutine print_single_spindet(string,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Subroutine to print the content of a determinant using the '+-' notation
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint)
  character*(2048)                :: output(1)

  call bitstring_to_str( output(1), string(1), Nint )
  print *,  trim(output(1))

end

logical function is_integer_in_string(bite,string,Nint)
  use bitmasks
 implicit none
 integer, intent(in) :: bite,Nint
 integer(bit_kind), intent(in) :: string(Nint)
 integer(bit_kind) :: string_bite(Nint)
 integer   :: i,itot,itot_and
 character*(2048)                :: output(1)
 string_bite = 0_bit_kind
 call set_bit_to_integer(bite,string_bite,Nint)
 itot = 0
 itot_and = 0
 is_integer_in_string = .False.
!print*,''
!print*,''
!print*,'bite = ',bite
!call bitstring_to_str( output(1), string_bite, Nint )
! print *,  trim(output(1))
!call bitstring_to_str( output(1), string, Nint )
! print *,  trim(output(1))
 do i = 1, Nint
  itot += popcnt(string(i))
  itot_and += popcnt(ior(string(i),string_bite(i)))
 enddo
!print*,'itot,itot_and',itot,itot_and
 if(itot == itot_and)then
  is_integer_in_string = .True. 
 endif
!pause 
end
