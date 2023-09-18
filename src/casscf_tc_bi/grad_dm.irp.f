 BEGIN_PROVIDER [real*8, gradvec_tc_r, (0:3,nMonoEx)]
&BEGIN_PROVIDER [real*8, gradvec_tc_l,  (0:3,nMonoEx)]
 BEGIN_DOC
! gradvec_tc_r(0:3,i) = <Chi_0| H E_q^p(i) |Phi_0> 
! 
! gradvec_tc_l(0:3,i) = <Chi_0| E_p^q(i) H |Phi_0> 
!
! where the indices "i" corresponds to E_q^p(i) 
! 
! i = mat_idx_c_a(q,p) 
!
! and gradvec_tc_r/l(0) = full matrix element
!
!     gradvec_tc_r/l(1) = one-body part 

!     gradvec_tc_r/l(2) = two-body part 

!     gradvec_tc_r/l(3) = three-body part 
 END_DOC
 implicit none
 integer                        :: ii,tt,aa,indx
 integer                        :: i,t,a,fff
 double precision :: res_l(0:3), res_r(0:3)
 gradvec_tc_l = 0.d0
 gradvec_tc_r = 0.d0
 ! computing the core/inactive --> virtual orbitals gradients
  do i=1,n_core_inact_orb
    ii=list_core_inact(i)
    do t=1,n_act_orb
      tt=list_act(t)
      indx = mat_idx_c_a(i,t) 
      call gradvec_tc_it(ii,tt,res_l,res_r)
      do fff = 0,3
       gradvec_tc_l(fff,indx)=res_l(fff)
       gradvec_tc_r(fff,indx)=res_r(fff)
      enddo
    end do
  end do
  
  do i=1,n_core_inact_orb
    ii=list_core_inact(i)
    do a=1,n_virt_orb
      indx = mat_idx_c_v(i,a) 
      aa=list_virt(a)
      call gradvec_tc_ia(ii,aa,res_l,res_r)
      do fff = 0,3
       gradvec_tc_l(fff,indx)=res_l(fff)
       gradvec_tc_r(fff,indx)=res_r(fff)
      enddo
    end do
  end do
  
! print*,'DM  grad'
  do t=1,n_act_orb
    tt=list_act(t)
    do a=1,n_virt_orb
      aa=list_virt(a)
      indx = mat_idx_a_v(t,a) 
!      print*,indx,t,a
      call gradvec_tc_ta(tt,aa,res_l, res_r)
      do fff = 0,3
       gradvec_tc_l(fff,indx)=res_l(fff)
       gradvec_tc_r(fff,indx)=res_r(fff)
      enddo
    end do
  end do
END_PROVIDER 

subroutine gradvec_tc_ia(i,a,res_l, res_r)
 implicit none
 BEGIN_DOC
! doubly occupied --> virtual TC gradient 
!
! Corresponds to res_r = <X0|H E_i^a|Phi_0>, 
!
!                res_l = <X0|E_a^i H|Phi_0>
 END_DOC
 integer, intent(in) :: i,a
 double precision, intent(out) :: res_l(0:3), res_r(0:3)
 res_l = 0.d0
 res_r = 0.d0
 res_l(1) = -2 * mo_bi_ortho_tc_one_e(a,i)
 res_r(1) = -2 * mo_bi_ortho_tc_one_e(i,a)
 integer :: j,t,r,jj,tt,rr
 do jj = 1, n_core_inact_orb
  j = list_core_inact(jj)
  res_r(2) += -2.d0 * ( 2.d0 * mo_bi_ortho_tc_two_e(j,i,j,a) - mo_bi_ortho_tc_two_e(i,j,j,a))
  res_l(2) -= -2.d0 * ( 2.d0 * mo_bi_ortho_tc_two_e(j,a,j,i) - mo_bi_ortho_tc_two_e(j,a,i,j))
 enddo
 do tt = 1, n_act_orb
  t = list_act(tt)
  do rr = 1, n_act_orb
   r = list_act(rr)
   res_r(2) += -0.5d0 * (                                                                                    &
   tc_transition_matrix_mo(r,t,1,1) *(2.d0 * mo_bi_ortho_tc_two_e(r,i,t,a) - mo_bi_ortho_tc_two_e(i,r,t,a))  & 
  +tc_transition_matrix_mo(t,r,1,1) *(2.d0 * mo_bi_ortho_tc_two_e(t,i,r,a) - mo_bi_ortho_tc_two_e(i,t,r,a))  &
   )
   res_l(2) -= -0.5d0 * (                                                                                    &
   tc_transition_matrix_mo(t,r,1,1) *(2.d0 * mo_bi_ortho_tc_two_e(t,a,r,i) - mo_bi_ortho_tc_two_e(t,a,i,r))  &
  +tc_transition_matrix_mo(r,t,1,1) *(2.d0 * mo_bi_ortho_tc_two_e(r,a,t,i) - mo_bi_ortho_tc_two_e(r,a,i,t))  & 
   )
  enddo
 enddo
 res_r(0) = res_r(1) + res_r(2) + res_r(3)
 res_l(0) = res_l(1) + res_l(2) + res_l(3)
end

subroutine gradvec_tc_it(i,t,res_l, res_r)
 implicit none
 BEGIN_DOC
! doubly occupied --> active TC gradient 
!
! Corresponds to res_r = <X0|H E_i^t|Phi_0>
!
!                res_l = <X0|E_t^i H |Phi_0>
 END_DOC
 integer, intent(in) :: i,t
 double precision, intent(out) :: res_l(0:3),res_r(0:3)
 integer :: rr,r,j,jj,u,uu,v,vv
 res_r = 0.d0
 res_l = 0.d0
 res_r(1)  += -2.d0 * mo_bi_ortho_tc_one_e(i,t) 
 res_l(1)  -= -2.D0 * mo_bi_ortho_tc_one_e(t,i) 

 do rr = 1, n_act_orb
  r = list_act(rr)
  res_r(1) +=  mo_bi_ortho_tc_one_e(i,r) * tc_transition_matrix_mo(t,r,1,1)   
  res_l(1) -=  mo_bi_ortho_tc_one_e(r,i) * tc_transition_matrix_mo(r,t,1,1)
 enddo

 do jj = 1, n_core_inact_orb
  j = list_core_inact(jj)
  res_r(2) += -2.d0 * (2d0 * mo_bi_ortho_tc_two_e(i,j,t,j) - mo_bi_ortho_tc_two_e(j,i,t,j)) 
  res_l(2) -= -2.d0 * (2d0 * mo_bi_ortho_tc_two_e(t,j,i,j) - mo_bi_ortho_tc_two_e(t,j,j,i)) 
  do rr = 1, n_act_orb
   r = list_act(rr)
   res_r(2) += tc_transition_matrix_mo(t,r,1,1) * (2.d0 * mo_bi_ortho_tc_two_e(i,j,r,j) - mo_bi_ortho_tc_two_e(i,j,j,r)) 
   res_l(2) -= tc_transition_matrix_mo(r,t,1,1) * (2.d0 * mo_bi_ortho_tc_two_e(r,j,i,j) - mo_bi_ortho_tc_two_e(j,r,j,i)) 
  enddo
 enddo
 do rr = 1, n_act_orb
  r = list_act(rr)
  do uu = 1, n_act_orb
   u = list_act(uu)
   res_r(2) += -0.5d0 * (  & 
    tc_transition_matrix_mo(u,r,1,1) * (2.d0 * mo_bi_ortho_tc_two_e(u,i,r,t) - mo_bi_ortho_tc_two_e(u,i,t,r))  &
  + tc_transition_matrix_mo(r,u,1,1) * (2.d0 * mo_bi_ortho_tc_two_e(i,r,t,u) - mo_bi_ortho_tc_two_e(i,r,u,t))  &
    )
   res_l(2) -= -0.5d0 * (  & 
    tc_transition_matrix_mo(r,u,1,1) * (2.d0 * mo_bi_ortho_tc_two_e(r,t,u,i) - mo_bi_ortho_tc_two_e(t,r,u,i))  &
  + tc_transition_matrix_mo(u,r,1,1) * (2.d0 * mo_bi_ortho_tc_two_e(t,u,i,r) - mo_bi_ortho_tc_two_e(u,t,i,r))  &
    )
   do vv = 1, n_act_orb
    v = list_act(vv)
    res_r(2) +=  0.5d0 * (  & 
    mo_bi_ortho_tc_two_e(i,r,v,u) * tc_two_rdm(t,r,v,u) + mo_bi_ortho_tc_two_e(r,i,v,u) * tc_two_rdm(r,t,v,u) )
    res_l(2) -=  0.5d0 * (  & 
    mo_bi_ortho_tc_two_e(v,u,i,r) * tc_two_rdm(v,u,t,r) + mo_bi_ortho_tc_two_e(v,u,r,i) * tc_two_rdm(v,u,r,t) )
   enddo
  enddo
 enddo
 res_r(0) = res_r(1) + res_r(2) + res_r(3)
 res_l(0) = res_l(1) + res_l(2) + res_l(3)
end

subroutine gradvec_tc_ta(t,a,res_l, res_r)
 implicit none
 BEGIN_DOC
! active --> virtual TC gradient 
!
! Corresponds to res_r = <X0|H E_t^a|Phi_0>
!
!                res_l = <X0|E_a^t H |Phi_0>
 END_DOC
 integer, intent(in) :: t,a
 double precision, intent(out) :: res_l(0:3),res_r(0:3)
 integer :: rr,r,j,jj,u,uu,v,vv
 double precision :: res_r_inact_test, res_r_act_test
 double precision :: res_l_inact_test, res_l_act_test
 res_r = 0.d0
 res_l = 0.d0
 do rr = 1, n_act_orb
  r = list_act(rr) 
  res_l(1) +=  mo_bi_ortho_tc_one_e(a,r) * tc_transition_matrix_mo(t,r,1,1)
  res_r(1) -=  mo_bi_ortho_tc_one_e(r,a) * tc_transition_matrix_mo(r,t,1,1)
 enddo

 res_r_inact_test = 0.d0
 res_l_inact_test = 0.d0
 do jj = 1, n_core_inact_orb
  j = list_core_inact(jj)
  do rr = 1, n_act_orb
   r = list_act(rr)
   res_r_inact_test += -tc_transition_matrix_mo(r,t,1,1) * & 
   (2.d0 * mo_bi_ortho_tc_two_e(r,j,a,j) - mo_bi_ortho_tc_two_e(r,j,j,a))
   res_l_inact_test -= -tc_transition_matrix_mo(t,r,1,1) * & 
   (2.d0 * mo_bi_ortho_tc_two_e(a,j,r,j) - mo_bi_ortho_tc_two_e(j,a,r,j))
  enddo
 enddo
 res_r_act_test = 0.d0
 res_l_act_test = 0.d0
 do rr = 1, n_act_orb
  r = list_act(rr)
  do vv = 1, n_act_orb
   v = list_act(vv)
   do uu = 1, n_act_orb
    u = list_act(uu)
    res_r_act_test += - (mo_bi_ortho_tc_two_e(v,r,u,a) * tc_two_rdm(r,v,t,u) & 
                        +mo_bi_ortho_tc_two_e(v,r,a,u) * tc_two_rdm(r,v,u,t))
    res_l_act_test -= - (mo_bi_ortho_tc_two_e(u,a,v,r) * tc_two_rdm(t,u,r,v) & 
                        +mo_bi_ortho_tc_two_e(a,u,v,r) * tc_two_rdm(u,t,r,v))
   enddo
  enddo
 enddo
 res_r_act_test *= 0.5d0
 res_l_act_test *= 0.5d0
 res_r(2) = res_r_inact_test + res_r_act_test
 res_l(2) = res_l_inact_test + res_l_act_test

 integer :: m,x,y
 double precision :: res_r_inact, res_r_act
 if(.False.)then
   ! test quantities
   res_r_inact = 0.d0
   res_r_act   = 0.d0
   do m = 1, mo_num
    do x = 1, mo_num
     do jj = 1, n_core_inact_orb
      j = list_core_inact(jj) 
      res_r_inact += 0.5d0 * mo_bi_ortho_tc_two_e(t,j,m,x) * tc_two_rdm(a,j,m,x) & 
                  -0.5d0 * mo_bi_ortho_tc_two_e(m,j,a,x) * tc_two_rdm(m,j,t,x) & 
                  +0.5d0 * mo_bi_ortho_tc_two_e(j,t,m,x) * tc_two_rdm(j,a,m,x) & 
                  -0.5d0 * mo_bi_ortho_tc_two_e(x,j,m,a) * tc_two_rdm(x,j,m,t)
     enddo
     do rr = 1, n_act_orb
      r = list_act(rr) 
      res_r_act   += 0.5d0 * mo_bi_ortho_tc_two_e(t,r,m,x) * tc_two_rdm(a,r,m,x) & 
                  -0.5d0 * mo_bi_ortho_tc_two_e(m,r,a,x) * tc_two_rdm(m,r,t,x) & 
                  +0.5d0 * mo_bi_ortho_tc_two_e(r,t,m,x) * tc_two_rdm(r,a,m,x) & 
                  -0.5d0 * mo_bi_ortho_tc_two_e(x,r,m,a) * tc_two_rdm(x,r,m,t)
     enddo
    enddo
   enddo
   if(dabs(res_r_inact).gt.1.d-12)then
    if(dabs(res_r_inact_test - res_r_inact).gt.1.d-10)then
     print*,'inact'
     print*,'t,a',t,a
     print*,res_r_inact_test , res_r_inact, dabs(res_r_inact_test - res_r_inact)
    endif
   endif
   if(dabs(res_r_act).gt.1.d-12)then
    if(dabs(res_r_act_test - res_r_act).gt.1.d-10)then
     print*,'act'
     print*,'t,a',t,a
     print*,res_r_act_test , res_r_act, dabs(res_r_act_test - res_r_act)
    endif
   endif
 endif

 res_r(0) = res_r(1) + res_r(2) + res_r(3)
 res_l(0) = res_l(1) + res_l(2) + res_l(3)
 
end
