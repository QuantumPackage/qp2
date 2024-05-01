 BEGIN_PROVIDER [integer ,  m_max_sm_7]
&BEGIN_PROVIDER [integer ,  n_max_sm_7]
&BEGIN_PROVIDER [integer ,  o_max_sm_7]
 implicit none
 BEGIN_DOC
! maximum value of the "m", "n" and "o" integer in the Jastrow function as in Eq. (4) 
! of Schmidt,Moskowitz, JCP, 93, 4172 (1990)  for the SM_7 version of Table IV
 END_DOC
 m_max_sm_7 = 4
 n_max_sm_7 = 0
 o_max_sm_7 = 4
END_PROVIDER 

 BEGIN_PROVIDER [integer ,  m_max_sm_9]
&BEGIN_PROVIDER [integer ,  n_max_sm_9]
&BEGIN_PROVIDER [integer ,  o_max_sm_9]
 implicit none
 BEGIN_DOC
! maximum value of the "m", "n" and "o" integer in the Jastrow function as in Eq. (4) 
! of Schmidt,Moskowitz, JCP, 93, 4172 (1990)  for the SM_9 version of Table IV
 END_DOC
 m_max_sm_9 = 4
 n_max_sm_9 = 2
 o_max_sm_9 = 4
END_PROVIDER 


 BEGIN_PROVIDER [integer ,  m_max_sm_17]
&BEGIN_PROVIDER [integer ,  n_max_sm_17]
&BEGIN_PROVIDER [integer ,  o_max_sm_17]
 implicit none
 BEGIN_DOC
! maximum value of the "m", "n" and "o" integer in the Jastrow function as in Eq. (4) 
! of Schmidt,Moskowitz, JCP, 93, 4172 (1990)  for the SM_17 version of Table IV
 END_DOC
 m_max_sm_17 = 6
 n_max_sm_17 = 2
 o_max_sm_17 = 6
END_PROVIDER 


BEGIN_PROVIDER [ double precision, c_mn_o_sm_7, (0:m_max_sm_7,0:n_max_sm_7,0:o_max_sm_7,2:10)]
 implicit none
 BEGIN_DOC
 !
 !c_mn_o_7(0:4,0:4,2:10) = coefficient for the SM_7 correlation factor as given is Table IV of 
 !                         Schmidt,Moskowitz, JCP, 93, 4172 (1990)
 !                         the first index (0:4) is the "m"  integer for the 1e part 
 !                         the second index(0:0) is the "n"  integer for the 1e part WHICH IS ALWAYS SET TO 0 FOR SM_7
 !                         the third index (0:4) is the "o" integer for the 2e part
 !                         the fourth index (2:10) is the nuclear charge of the atom 
 END_DOC
 c_mn_o_sm_7 = 0.d0
 integer :: i
 do i = 2, 10 ! loop over nuclear charge
  c_mn_o_sm_7(0,0,1,i) = 0.5d0 ! all the linear terms are set to 1/2 to satisfy the anti-parallel spin condition
 enddo
 ! He atom 
 ! two electron terms 
 c_mn_o_sm_7(0,0,2,2) =  0.50516d0 
 c_mn_o_sm_7(0,0,3,2) = -0.19313d0 
 c_mn_o_sm_7(0,0,4,2) =  0.30276d0 
 ! one-electron terms 
 c_mn_o_sm_7(2,0,0,2) = -0.16995d0
 c_mn_o_sm_7(3,0,0,2) = -0.34505d0
 c_mn_o_sm_7(4,0,0,2) = -0.54777d0
 ! Ne atom 
 ! two electron terms 
 c_mn_o_sm_7(0,0,2,10) = -0.792d0
 c_mn_o_sm_7(0,0,3,10) =  1.05232d0 
 c_mn_o_sm_7(0,0,4,10) = -0.65615d0 
 ! one-electron terms 
 c_mn_o_sm_7(2,0,0,10) = -0.13312d0
 c_mn_o_sm_7(3,0,0,10) = -0.00131d0
 c_mn_o_sm_7(4,0,0,10) =  0.09083d0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, c_mn_o_sm_9, (0:m_max_sm_9,0:n_max_sm_9,0:o_max_sm_9,2:10)]
 implicit none
 BEGIN_DOC
 !
 !c_mn_o_9(0:4,0:4,2:10) = coefficient for the SM_9 correlation factor as given is Table IV of 
 !                         Schmidt,Moskowitz, JCP, 93, 4172 (1990)
 !                         the first index (0:4) is the "m"  integer for the 1e part 
 !                         the second index(0:0) is the "n"  integer for the 1e part WHICH IS ALWAYS SET TO 0 FOR SM_9
 !                         the third index (0:4) is the "o" integer for the 2e part
 !                         the fourth index (2:10) is the nuclear charge of the atom 
 END_DOC
 c_mn_o_sm_9 = 0.d0
 integer :: i
 do i = 2, 10 ! loop over nuclear charge
  c_mn_o_sm_9(0,0,1,i) = 0.5d0 ! all the linear terms are set to 1/2 to satisfy the anti-parallel spin condition
 enddo
 ! He atom 
 ! two electron terms 
 c_mn_o_sm_9(0,0,2,2) =  0.50516d0 
 c_mn_o_sm_9(0,0,3,2) = -0.19313d0 
 c_mn_o_sm_9(0,0,4,2) =  0.30276d0 
 ! one-electron terms 
 c_mn_o_sm_9(2,0,0,2) = -0.16995d0
 c_mn_o_sm_9(3,0,0,2) = -0.34505d0
 c_mn_o_sm_9(4,0,0,2) = -0.54777d0
 ! Ne atom 
 ! two electron terms 
 c_mn_o_sm_9(0,0,2,10) = -0.792d0
 c_mn_o_sm_9(0,0,3,10) =  1.05232d0 
 c_mn_o_sm_9(0,0,4,10) = -0.65615d0 
 ! one-electron terms 
 c_mn_o_sm_9(2,0,0,10) = -0.13312d0
 c_mn_o_sm_9(3,0,0,10) = -0.00131d0
 c_mn_o_sm_9(4,0,0,10) =  0.09083d0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, c_mn_o_sm_17, (0:m_max_sm_17,0:n_max_sm_17,0:o_max_sm_17,2:10)]
 implicit none
 BEGIN_DOC
 !
 !c_mn_o_17(0:4,0:4,2:10) = coefficient for the SM_17 correlation factor as given is Table IV of 
 !                         Schmidt,Moskowitz, JCP, 93, 4172 (1990)
 !                         the first index (0:4) is the "m"  integer for the 1e part 
 !                         the second index(0:0) is the "n"  integer for the 1e part WHICH IS ALWAYS SET TO 0 FOR SM_17
 !                         the third index (0:4) is the "o" integer for the 2e part
 !                         the fourth index (2:10) is the nuclear charge of the atom 
 END_DOC
 c_mn_o_sm_17 = 0.d0
 integer :: i
 do i = 2, 10 ! loop over nuclear charge
  c_mn_o_sm_17(0,0,1,i) = 0.5d0 ! all the linear terms are set to 1/2 to satisfy the anti-parallel spin condition
 enddo
 ! He atom 
 ! two electron terms 
 c_mn_o_sm_17(0,0,2,2) =  0.09239d0 
 c_mn_o_sm_17(0,0,3,2) = -0.38664d0 
 c_mn_o_sm_17(0,0,4,2) =  0.95764d0 
 ! one-electron terms 
 c_mn_o_sm_17(2,0,0,2) =  0.23208d0
 c_mn_o_sm_17(3,0,0,2) = -0.45032d0
 c_mn_o_sm_17(4,0,0,2) =  0.82777d0
 c_mn_o_sm_17(2,2,0,2) = -4.15388d0
 ! ee-n terms 
 c_mn_o_sm_17(2,0,2,2) =  0.80622d0
 c_mn_o_sm_17(2,2,2,2) = 10.19704d0
 c_mn_o_sm_17(4,0,2,2) = -4.96259d0
 c_mn_o_sm_17(2,0,4,2) = -1.35647d0
 c_mn_o_sm_17(4,2,2,2) = -5.90907d0
 c_mn_o_sm_17(6,0,2,2) =  0.90343d0
 c_mn_o_sm_17(4,0,4,2) =  5.50739d0
 c_mn_o_sm_17(2,2,4,2) = -0.03154d0
 c_mn_o_sm_17(2,0,6,2) = -1.1051860


 ! Ne atom 
 ! two electron terms 
 c_mn_o_sm_17(0,0,2,10) = -0.80909d0 
 c_mn_o_sm_17(0,0,3,10) = -0.00219d0 
 c_mn_o_sm_17(0,0,4,10) =  0.59188d0
 ! one-electron terms 
 c_mn_o_sm_17(2,0,0,10) = -0.00567d0
 c_mn_o_sm_17(3,0,0,10) =  0.14011d0
 c_mn_o_sm_17(4,0,0,10) = -0.05671d0
 c_mn_o_sm_17(2,2,0,10) = -3.33767d0
 ! ee-n terms 
 c_mn_o_sm_17(2,0,2,10) =  1.95067d0
 c_mn_o_sm_17(2,2,2,10) =  6.83340d0
 c_mn_o_sm_17(4,0,2,10) = -3.29231d0
 c_mn_o_sm_17(2,0,4,10) = -2.44998d0
 c_mn_o_sm_17(4,2,2,10) = -2.13029d0
 c_mn_o_sm_17(6,0,2,10) =  2.25768d0
 c_mn_o_sm_17(4,0,4,10) =  1.97951d0
 c_mn_o_sm_17(2,2,4,10) = -2.0924160
 c_mn_o_sm_17(2,0,6,10) =  0.35493d0

END_PROVIDER 

 BEGIN_PROVIDER [ double precision, b_I_sm_90,(2:10)]
&BEGIN_PROVIDER [ double precision, d_I_sm_90,(2:10)]
 implicit none
 BEGIN_DOC
! "b_I" and "d_I" parameters of Eqs. (4) and (5) of Schmidt,Moskowitz, JCP, 93, 4172 (1990)
 END_DOC
 b_I_sm_90 = 1.d0
 d_I_sm_90 = 1.d0
 
END_PROVIDER 

subroutine get_full_sm_90_jastrow(r1,r2,rI,sm_j,i_charge, j_1e,j_2e,j_een,j_tot)
 implicit none 
 double precision, intent(in) :: r1(3),r2(3),rI(3)
 integer, intent(in)          :: sm_j, i_charge 
 double precision, intent(out):: j_1e,j_2e,j_een,j_tot
 BEGIN_DOC
 ! Jastrow function as in Eq. (4) of Schmidt,Moskowitz, JCP, 93, 4172 (1990)
 ! the i_charge variable is the integer specifying the charge of the atom for the Jastrow 
 ! the sm_j integer variable represents the "quality" of the jastrow : sm_j = 7, 9, 17
 END_DOC 
 double precision :: r_inucl,r_jnucl,r_ij,b_I, d_I
 b_I = b_I_sm_90(i_charge)
 d_I = d_I_sm_90(i_charge)
 call get_rescaled_variables_j_sm_90(r1,r2,rI,b_I,d_I,r_inucl,r_jnucl,r_ij)
 call jastrow_func_sm_90(r_inucl,r_jnucl,r_ij,sm_j,i_charge, j_1e,j_2e,j_een,j_tot)
end

subroutine get_rescaled_variables_j_sm_90(r1,r2,rI,b_I,d_I,r_inucl,r_jnucl,r_ij)
 implicit none
 BEGIN_DOC
 ! rescaled variables of Eq. (5) and (6) of Schmidt,Moskowitz, JCP, 93, 4172 (1990)
 ! the "b_I" and "d_I" parameters are the same as in Eqs. (5) and (6) 
 END_DOC
 double precision, intent(in) :: r1(3),r2(3),rI(3)
 double precision, intent(in) :: b_I, d_I
 double precision, intent(out):: r_inucl,r_jnucl,r_ij
 double precision :: rin, rjn, rij
 integer :: i
 rin = 0.d0
 rjn = 0.d0
 rij = 0.d0
 do i = 1,3
  rin += (r1(i) - rI(i)) * (r1(i) - rI(i))
  rjn += (r2(i) - rI(i)) * (r2(i) - rI(i))
  rij += (r2(i) - r1(i)) * (r2(i) - r1(i))
 enddo
 rin = dsqrt(rin)
 rjn = dsqrt(rjn)
 rij = dsqrt(rij)
 r_inucl = b_I * rin/(1.d0 + b_I * rin)
 r_jnucl = b_I * rjn/(1.d0 + b_I * rjn)
 r_ij    = d_I * rij/(1.d0 + b_I * rij)
end

subroutine jastrow_func_sm_90(r_inucl,r_jnucl,r_ij,sm_j,i_charge, j_1e,j_2e,j_een,j_tot)
 implicit none
 BEGIN_DOC
 ! Jastrow function as in Eq. (4) of Schmidt,Moskowitz, JCP, 93, 4172 (1990)
 ! Here the r_inucl, r_jnucl are the rescaled variables as defined in Eq. (5) with "b_I" 
 !          r_ij is the rescaled variable as defined in Eq. (6) with "d_I" 
 ! the i_charge variable is the integer specifying the charge of the atom for the Jastrow 
 ! the sm_j integer variable represents the "quality" of the jastrow : sm_j = 7, 9, 17
 !
 ! it returns the j_1e  : sum of terms with "o" = "n" = 0, "m" /= 0,
 !                j_2e  : sum of terms with "m" = "n" = 0, "o" /= 0,
 !                j_een : sum of terms with "m" /=0, "n" /= 0, "o" /= 0,
 !                j_tot : the total sum 
 END_DOC
 double precision, intent(in) :: r_inucl,r_jnucl,r_ij
 integer, intent(in)          :: sm_j,i_charge
 double precision, intent(out):: j_1e,j_2e,j_een,j_tot
 j_1e  = 0.D0
 j_2e  = 0.D0
 j_een = 0.D0
 double precision :: delta_mn,jastrow_sm_90_atomic
 integer :: m,n,o
BEGIN_TEMPLATE
 ! pure 2e part 
 n = 0
 m = 0 
 if(sm_j == $X )then
  do o = 1, o_max_sm_$X
   if(dabs(c_mn_o_sm_$X(m,n,o,i_charge)).lt.1.d-10)cycle 
   j_2e += c_mn_o_sm_$X(m,n,o,i_charge) * jastrow_sm_90_atomic(m,n,o,i_charge,r_inucl,r_jnucl,r_ij)
  enddo
! else 
!  print*,'sm_j = ',sm_j
!  print*,'not implemented, stop'
!  stop
 endif
 ! pure one-e part 
 o = 0 
 if(sm_j == $X)then
  do n = 2, n_max_sm_$X
   do m = 2, m_max_sm_$X
    j_1e += c_mn_o_sm_$X(m,n,o,i_charge) * jastrow_sm_90_atomic(m,n,o,i_charge,r_inucl,r_jnucl,r_ij)
   enddo
  enddo
! else 
!  print*,'sm_j = ',sm_j
!  print*,'not implemented, stop'
!  stop
 endif
 ! e-e-n part 
 if(sm_j == $X)then
  do o = 1, o_max_sm_$X
   do m = 2, m_max_sm_$X
    do n = 2, n_max_sm_$X
     j_een += c_mn_o_sm_$X(m,n,o,i_charge) * jastrow_sm_90_atomic(m,n,o,i_charge,r_inucl,r_jnucl,r_ij)
    enddo
   enddo
  enddo
 else 
!  print*,'sm_j = ',sm_j
!  print*,'not implemented, stop'
!  stop
 endif
 j_tot = j_1e + j_2e + j_een
SUBST [ X]
  7 ;; 
  9 ;; 
  17 ;;
END_TEMPLATE 
end

double precision function jastrow_sm_90_atomic(m,n,o,i_charge,r_inucl,r_jnucl,r_ij)
 implicit none
 BEGIN_DOC
! contribution to the function of Eq. (4) of Schmidt,Moskowitz, JCP, 93, 4172 (1990)
! for a given m,n,o and atom 
 END_DOC
 double precision, intent(in) :: r_inucl,r_jnucl,r_ij
 integer         , intent(in) :: m,n,o,i_charge
 double precision :: delta_mn 
 if(m==n)then
  delta_mn = 0.5d0
 else
  delta_mn = 1.D0
 endif
 jastrow_sm_90_atomic = delta_mn * (r_inucl**m * r_jnucl**n + r_jnucl**m * r_inucl**n)*r_ij**o
end
