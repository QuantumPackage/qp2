double precision function plgndr(l,m,x)
 integer, intent(in) :: l,m
 double precision, intent(in) :: x
 BEGIN_DOC
 ! associated Legenre polynom P_l,m(x). Used for the Y_lm(theta,phi)
 ! Taken from https://iate.oac.uncor.edu/~mario/materia/nr/numrec/f6-8.pdf
 END_DOC
 integer :: i,ll
 double precision :: fact,pll,pmm,pmmp1,somx2
 if(m.lt.0.or.m.gt.l.or.dabs(x).gt.1.d0)then 
  print*,'bad arguments in plgndr'
  pause
 endif
 pmm=1.d0
 if(m.gt.0) then
  somx2=dsqrt((1.d0-x)*(1.d0+x))
  fact=1.d0
  do i=1,m
   pmm=-pmm*fact*somx2
   fact=fact+2.d0
  enddo 
 endif ! m > 0
 if(l.eq.m) then
  plgndr=pmm 
 else
  pmmp1=x*(2*m+1)*pmm ! Compute P_m+1^m
   if(l.eq.m+1) then
    plgndr=pmmp1
   else ! Compute P_l^m, l> m+1
    do ll=m+2,l
     pll=(x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm)/(ll-m)
     pmm=pmmp1
     pmmp1=pll
    enddo
    plgndr=pll
   endif ! l.eq.m+1
 endif ! l.eq.m
 return
end

double precision function ortho_assoc_gaus_pol(l1,m1,l2)
 implicit none
 integer, intent(in) :: l1,m1,l2
 double precision :: fact
 if(l1.ne.l2)then
  ortho_assoc_gaus_pol= 0.d0
 else
  ortho_assoc_gaus_pol = 2.d0*fact(l1+m1) / (dble(2*l1+1)*fact(l1-m1))
 endif
end
