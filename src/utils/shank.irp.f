double precision function shank3_f(array,n,nmax)
 implicit none                                                                                                                                                         
 integer, intent(in) :: n,nmax
 double precision, intent(in) :: array(0:nmax) ! array of the partial sums
 integer :: ntmp
 double precision :: shank1(0:nmax),shank2(0:nmax),shank3(0:nmax)
 ntmp = n
 call shank(array,ntmp,nmax,shank1)
 ntmp = ntmp - 2
 call shank(shank1,ntmp,nmax,shank2)
 ntmp = ntmp - 2
 call shank(shank2,ntmp,nmax,shank3)
 ntmp = ntmp - 2
 shank3_f = shank3(ntmp)
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


