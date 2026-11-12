integer function dress_find_sample(v, w)
  BEGIN_DOC
  ! Binary search on the cumulative weight array
  ! w(0:N_det_generators), assumed non-decreasing: returns the
  ! index i such that w(i-1) < v <= w(i), i.e. the generator
  ! selected by the sample value v. If several consecutive
  ! cumulative values are equal (zero-width bin), the last index
  ! of the run of equal values is returned.
  ! Probably to avoid a division by zero or a zero-width bin
  ! (author, 2026-10-31); the code shows the returned index is a
  ! lookup key into the generator arrays (edI(i), pt2_w(i)) and
  ! the selected weight pt2_w(i) is the divisor of the estimator
  ! x += edI(i) * pt2_W_T / pt2_w(i).
  !
  ! Moved verbatim from the removed dressing module (2026-11-12);
  ! the name is kept so that the external declaration in
  ! pt2_stoch_routines stays untouched.
  END_DOC
  implicit none
  double precision, intent(in) :: v, w(0:N_det_generators)
  integer :: i,l,r

  l = 0
  r = N_det_generators

  do while(r-l > 1)
    i = shiftr(r+l,1)
    if(w(i) < v) then
      l = i
    else
      r = i
    end if
  end do
  i = r
  do r=i+1,N_det_generators
    if (w(r) /= w(i)) then
      exit
    endif
  enddo
  dress_find_sample = r-1
end function
