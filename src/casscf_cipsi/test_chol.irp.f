program test_chol
 implicit none
 read_wf= .True.
 touch read_wf 
! call routine_bielec_PxxQ_no
! call routine_bielecCI_no
! call test_bielec_PxxQ_chol
! call test_bielecCI

end

subroutine routine_bielec_PQxx_no
 implicit none
 integer :: i_chol, i_act, ii_act, j_act, jj_act, i_core_inact, j_core_inact, ii_core_inact, jj_core_inact
 integer :: i_virt, ii_virt, j_virt, jj_virt, i_mo, j_mo
 double precision :: exact, new, error, accu, bielec_no_basis_chol
 double precision :: bielec_PQxx_no
 
 accu = 0.d0
 do i_core_inact = 1, n_core_inact_act_orb
  ii_core_inact = list_core_inact_act(i_core_inact)
  do j_core_inact = 1, n_core_inact_act_orb
   jj_core_inact = list_core_inact_act(j_core_inact)
   do i_mo = 1, mo_num
    do j_mo = 1, mo_num
     exact = bielec_PQxx_no_array(j_mo,i_mo, j_core_inact, i_core_inact) 
     new   = bielec_PQxx_no(j_mo,i_mo, j_core_inact, i_core_inact) 
     error = dabs(exact-new)
     if(dabs(exact).gt.1.d-10)then
      print*,exact,new,error
     endif
     accu += error
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/(dble(mo_num*mo_num*n_core_inact_act_orb**2))
end

subroutine routine_bielec_PxxQ_no_array
 implicit none
 integer :: i_chol, i_act, ii_act, j_act, jj_act, i_core_inact, j_core_inact, ii_core_inact, jj_core_inact
 integer :: i_virt, ii_virt, j_virt, jj_virt, i_mo, j_mo
 double precision :: exact, new, error, accu, bielec_no_basis_chol
 double precision :: bielec_PxxQ_no
 
 accu = 0.d0
 do i_mo = 1, mo_num
  do i_core_inact = 1, n_core_inact_act_orb
  ii_core_inact = list_core_inact_act(i_core_inact)
   do j_core_inact = 1, n_core_inact_act_orb
   jj_core_inact = list_core_inact_act(j_core_inact)
    do j_mo = 1, mo_num
     exact = bielec_PxxQ_no_array(j_mo, j_core_inact,  i_core_inact,i_mo) 
!     new   = bielec_no_basis_chol(j_mo,i_mo, jj_core_inact, ii_core_inact) 
     new   = bielec_PxxQ_no(j_mo, j_core_inact,  i_core_inact,i_mo) 
     error = dabs(exact-new)
     accu += error
     if(dabs(exact).gt.1.d-10)then
      print*,exact,new,error
     endif
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/(dble(mo_num*mo_num*n_core_inact_act_orb**2))
end

subroutine test_bielec_PQxx(i_mo, j_mo, i_ca, j_ca)
 implicit none
 integer :: i_mo, j_mo, i_ca, j_ca 
 double precision :: exact, new, error, accu
 double precision :: bielec_PQxx
 
 accu = 0.d0
 do j_ca = 1, n_core_inact_act_orb
  do i_ca = 1, n_core_inact_act_orb
   do j_mo = 1, mo_num
    do i_mo = 1, mo_num
     exact = bielec_PQxx_array(i_mo, j_mo, i_ca, j_ca)
     new   = bielec_PQxx(i_mo, j_mo, i_ca, j_ca)
     error = dabs(exact-new)
     accu += error
     if(dabs(exact).gt.1.d-10)then
      print*,exact,new,error
     endif
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/(dble(mo_num*mo_num*n_core_inact_act_orb**2))
end

subroutine test_bielec_PxxQ_chol(i_mo, i_ca, j_ca, j_mo)
 implicit none
 integer :: i_mo, i_ca, j_ca, j_mo
 double precision :: exact, new, error, accu
 double precision :: bielec_PxxQ
 accu = 0.d0
 do j_mo = 1, mo_num
  do j_ca = 1, n_core_inact_act_orb
   do i_ca =1, n_core_inact_act_orb
    do i_mo = 1, mo_num
     exact = bielec_PxxQ_array(i_mo, i_ca, j_ca, j_mo)
     new   = bielec_PxxQ(i_mo, i_ca, j_ca, j_mo)
     error = dabs(exact-new)
     accu += error
     if(dabs(exact).gt.1.d-10)then
      print*,exact,new,error
     endif
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/(dble(mo_num*mo_num*n_core_inact_act_orb**2))
end
