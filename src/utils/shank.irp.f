double precision function shank_general(array,n,nmax)
 implicit none                                                                                                                                                         
 integer, intent(in) :: n,nmax
 double precision, intent(in) :: array(0:nmax) ! array of the partial sums
 integer :: ntmp,i
 double precision :: sum(0:nmax),shank1(0:nmax)
 if(n.lt.3)then
  print*,'You asked to Shank a sum but the order is smaller than 3 ...'
  print*,'n = ',n
  print*,'stopping ....'
  stop
 endif
 ntmp = n
 sum = array
 i = 0
 do while(ntmp.ge.2)
  i += 1
!  print*,'i = ',i
  call shank(sum,ntmp,nmax,shank1)
  ntmp = ntmp - 2
  sum = shank1
  shank_general = shank1(ntmp)
 enddo
end


subroutine shank(array,n,nmax,shank1)
 implicit none
 integer, intent(in) :: n,nmax
 double precision, intent(in)  :: array(0:nmax)
 double precision, intent(out) :: shank1(0:nmax)
 integer :: i,j
 double precision :: shank_function
 do i = 1, n-1
  shank1(i-1) = shank_function(array,i,nmax)
 enddo
end

double precision function shank_function(array,i,n)
 implicit none
 integer, intent(in) :: i,n
 double precision, intent(in) :: array(0:n)
 double precision :: b_n, b_n1
 b_n = array(i) - array(i-1)
 b_n1 = array(i+1) - array(i)
 if(dabs(b_n1-b_n).lt.1.d-12)then
  shank_function = array(i+1)
 else
  shank_function = array(i+1) - b_n1*b_n1/(b_n1-b_n)
 endif

end


