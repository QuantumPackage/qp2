! ---

program bi_ort_ints

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

! call test_3e
! call test_5idx
! call test_5idx2
!  call test_4idx()
  !call test_4idx_n4()
  !call test_4idx2()
  !call test_5idx2
  !call test_5idx

  call test_mos_in_r()
  call test_int2_grad1_u12_bimo_t()

end

subroutine test_5idx2
  PROVIDE three_e_5_idx_cycle_2_bi_ort
end

subroutine test_4idx2()
  !PROVIDE three_e_4_idx_direct_bi_ort 
  PROVIDE three_e_4_idx_exch23_bi_ort
end

subroutine test_3e
 implicit none
 integer :: i,k,j,l,m,n,ipoint
 double precision :: accu, contrib,new,ref
 i = 1
 k = 1
 n = 0
 accu = 0.d0
 do i = 1, mo_num
  do k = 1, mo_num
   do j = 1, mo_num
    do l = 1, mo_num
     do m = 1, mo_num
      do n = 1, mo_num
        call give_integrals_3_body_bi_ort(n, l, k, m, j, i, new)
        call give_integrals_3_body_bi_ort_old(n, l, k, m, j, i, ref)
        contrib = dabs(new - ref)
        accu += contrib
        if(contrib .gt. 1.d-10)then
         print*,'pb !!'
         print*,i,k,j,l,m,n
         print*,ref,new,contrib
         stop
        endif
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(mo_num)**6


end

subroutine test_5idx
 implicit none
 integer :: i,k,j,l,m,n,ipoint
 double precision :: accu, contrib,new,ref
 double precision, external :: three_e_5_idx_exch12_bi_ort
 i = 1
 k = 1
 n = 0
 accu = 0.d0
 PROVIDE three_e_5_idx_direct_bi_ort_old

 do i = 1, mo_num
  do k = 1, mo_num
   do j = 1, mo_num
    do l = 1, mo_num
     do m = 1, mo_num
!      if (dabs(three_e_5_idx_direct_bi_ort(m,l,j,k,i) - three_e_5_idx_exch12_bi_ort(m,l,i,k,j)) > 1.d-10) then 
!         stop
!      endif
      new = three_e_5_idx_direct_bi_ort(m,l,j,k,i)
      ref = three_e_5_idx_direct_bi_ort_old(m,l,j,k,i)
      contrib = dabs(new - ref)
      accu += contrib
      if(contrib .gt. 1.d-10)then
       print*,'direct'
       print*,i,k,j,l,m
       print*,ref,new,contrib
       stop
      endif
!
!      new = three_e_5_idx_exch12_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_exch12_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'exch12'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!
!      new = three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_cycle_1_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'cycle1'

!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!      new = three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_cycle_2_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'cycle2'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!      new = three_e_5_idx_exch23_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_exch23_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'exch23'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!      new = three_e_5_idx_exch13_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_exch13_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'exch13'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!      new = three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_cycle_1_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'cycle1'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!      new = three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_cycle_2_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'cycle2'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!      new = three_e_5_idx_exch23_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_exch23_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'exch23'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
!      new = three_e_5_idx_exch13_bi_ort(m,l,j,k,i)
!      ref = three_e_5_idx_exch13_bi_ort_old(m,l,j,k,i)
!      contrib = dabs(new - ref)
!      accu += contrib
!      if(contrib .gt. 1.d-10)then
!       print*,'exch13'
!       print*,i,k,j,l,m
!       print*,ref,new,contrib
!       stop
!      endif
!
     enddo
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(mo_num)**5


end

! ---

subroutine test_4idx_n4()

  implicit none
  integer          :: i, j, k, l
  double precision :: accu, contrib, new, ref, thr

  thr = 1d-10

  PROVIDE three_e_4_idx_direct_bi_ort_old
  PROVIDE three_e_4_idx_direct_bi_ort_n4

  accu = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = three_e_4_idx_direct_bi_ort_n4 (l,k,j,i)
          ref = three_e_4_idx_direct_bi_ort_old(l,k,j,i)
          contrib = dabs(new - ref)
          accu += contrib
          if(contrib .gt. thr) then
            print*, ' problem in three_e_4_idx_direct_bi_ort_n4'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

        enddo
      enddo
    enddo
  enddo
  print*, ' accu on three_e_4_idx_direct_bi_ort_n4 = ', accu / dble(mo_num)**4

  ! ---

  PROVIDE three_e_4_idx_exch13_bi_ort_old
  PROVIDE three_e_4_idx_exch13_bi_ort_n4

  accu = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = three_e_4_idx_exch13_bi_ort_n4 (l,k,j,i)
          ref = three_e_4_idx_exch13_bi_ort_old(l,k,j,i)
          contrib = dabs(new - ref)
          accu += contrib
          if(contrib .gt. thr) then
            print*, ' problem in three_e_4_idx_exch13_bi_ort_n4'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

        enddo
      enddo
    enddo
  enddo
  print*, ' accu on three_e_4_idx_exch13_bi_ort_n4 = ', accu / dble(mo_num)**4

  ! ---

  PROVIDE three_e_4_idx_cycle_1_bi_ort_old
  PROVIDE three_e_4_idx_cycle_1_bi_ort_n4

  accu = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = three_e_4_idx_cycle_1_bi_ort_n4 (l,k,j,i)
          ref = three_e_4_idx_cycle_1_bi_ort_old(l,k,j,i)
          contrib = dabs(new - ref)
          accu += contrib
          if(contrib .gt. thr) then
            print*, ' problem in three_e_4_idx_cycle_1_bi_ort_n4'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

        enddo
      enddo
    enddo
  enddo
  print*, ' accu on three_e_4_idx_cycle_1_bi_ort_n4 = ', accu / dble(mo_num)**4

  ! ---

  PROVIDE three_e_4_idx_exch23_bi_ort_old
  PROVIDE three_e_4_idx_exch23_bi_ort_n4

  accu = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = three_e_4_idx_exch23_bi_ort_n4  (l,k,j,i)
          ref = three_e_4_idx_exch23_bi_ort_old(l,k,j,i)
          contrib = dabs(new - ref)
          accu += contrib
          if(contrib .gt. thr) then
            print*, ' problem in three_e_4_idx_exch23_bi_ort_n4'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

        enddo
      enddo
    enddo
  enddo
  print*, ' accu on three_e_4_idx_exch23_bi_ort_n4 = ', accu / dble(mo_num)**4

  ! ---

  return
end

! ---

subroutine test_4idx()

  implicit none
  integer          :: i, j, k, l
  double precision :: accu, contrib, new, ref, thr, norm

  thr = 1d-10

  PROVIDE three_e_4_idx_direct_bi_ort_old
  PROVIDE three_e_4_idx_direct_bi_ort 

  accu = 0.d0
  norm = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = three_e_4_idx_direct_bi_ort    (l,k,j,i)
          ref = three_e_4_idx_direct_bi_ort_old(l,k,j,i)
          contrib = dabs(new - ref)
          if(contrib .gt. thr) then
            print*, ' problem in three_e_4_idx_direct_bi_ort'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

          accu += contrib
          norm += dabs(ref)
        enddo
      enddo
    enddo
  enddo

  print*, ' accu on three_e_4_idx_direct_bi_ort (%) = ', 100.d0 * accu / norm

  ! ---

  PROVIDE three_e_4_idx_exch13_bi_ort_old
  PROVIDE three_e_4_idx_exch13_bi_ort 

  accu = 0.d0
  norm = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = three_e_4_idx_exch13_bi_ort   (l,k,j,i)
          ref = three_e_4_idx_exch13_bi_ort_old(l,k,j,i)
          contrib = dabs(new - ref)
          if(contrib .gt. thr) then
            print*, ' problem in three_e_4_idx_exch13_bi_ort'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

          accu += contrib
          norm += dabs(ref)
        enddo
      enddo
    enddo
  enddo

  print*, ' accu on three_e_4_idx_exch13_bi_ort (%) = ', 100.d0 * accu / norm

  ! ---

  PROVIDE three_e_4_idx_cycle_1_bi_ort_old
  PROVIDE three_e_4_idx_cycle_1_bi_ort

  accu = 0.d0
  norm = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = three_e_4_idx_cycle_1_bi_ort    (l,k,j,i)
          ref = three_e_4_idx_cycle_1_bi_ort_old(l,k,j,i)
          contrib = dabs(new - ref)
          if(contrib .gt. thr) then
            print*, ' problem in three_e_4_idx_cycle_1_bi_ort'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

          accu += contrib
          norm += dabs(ref)
        enddo
      enddo
    enddo
  enddo

  print*, ' accu on three_e_4_idx_cycle_1_bi_ort (%) = ', 100.d0 * accu / norm

  ! ---

  PROVIDE three_e_4_idx_exch23_bi_ort_old
  PROVIDE three_e_4_idx_exch23_bi_ort

  accu = 0.d0
  norm = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = three_e_4_idx_exch23_bi_ort    (l,k,j,i)
          ref = three_e_4_idx_exch23_bi_ort_old(l,k,j,i)
          contrib = dabs(new - ref)
          if(contrib .gt. thr) then
            print*, ' problem in three_e_4_idx_exch23_bi_ort'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

          accu += contrib
          norm += dabs(ref)
        enddo
      enddo
    enddo
  enddo

  print*, ' accu on three_e_4_idx_exch23_bi_ort (%) = ', 100.d0 * accu / norm

  ! ---

  return
end

! ---

subroutine test_mos_in_r()

  implicit none

  integer          :: i, j
  double precision :: err_tot, nrm_tot, err_loc, acc_thr

  PROVIDE mos_l_in_r_array_transp_old mos_r_in_r_array_transp_old
  PROVIDE mos_l_in_r_array_transp mos_r_in_r_array_transp

  acc_thr = 1d-13

  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do j = 1, n_points_final_grid
      err_loc = dabs(mos_l_in_r_array_transp_old(j,i) - mos_l_in_r_array_transp(j,i))
      if(err_loc > acc_thr) then
        print*, " error on", j, i
        print*, " old res", mos_l_in_r_array_transp_old(j,i)
        print*, " new res", mos_l_in_r_array_transp    (j,i)
        stop
      endif
      err_tot = err_tot + err_loc
      nrm_tot = nrm_tot + dabs(mos_l_in_r_array_transp_old(j,i))
    enddo
  enddo
  print *, ' absolute accuracy on mos_l_in_r_array_transp (%) =', 100.d0 * err_tot / nrm_tot

  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do j = 1, n_points_final_grid
      err_loc = dabs(mos_r_in_r_array_transp_old(j,i) - mos_r_in_r_array_transp(j,i))
      if(err_loc > acc_thr) then
        print*, " error on", j, i
        print*, " old res", mos_r_in_r_array_transp_old(j,i)
        print*, " new res", mos_r_in_r_array_transp    (j,i)
        stop
      endif
      err_tot = err_tot + err_loc
      nrm_tot = nrm_tot + dabs(mos_r_in_r_array_transp_old(j,i))
    enddo
  enddo
  print *, ' absolute accuracy on mos_r_in_r_array_transp (%) =', 100.d0 * err_tot / nrm_tot

  return
end

! ---

subroutine test_int2_grad1_u12_bimo_t()

  implicit none
  integer          :: i, j, ipoint, m
  double precision :: err_tot, nrm_tot, err_loc, acc_thr

  PROVIDE int2_grad1_u12_bimo_t_old
  PROVIDE int2_grad1_u12_bimo_t

  acc_thr = 1d-13

  err_tot = 0.d0
  nrm_tot = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do m = 1, 3
        do ipoint = 1, n_points_final_grid
          err_loc = dabs(int2_grad1_u12_bimo_t_old(ipoint,m,j,i) - int2_grad1_u12_bimo_t(ipoint,m,j,i))
          if(err_loc > acc_thr) then
            print*, " error on", ipoint, m, j, i
            print*, " old res", int2_grad1_u12_bimo_t_old(ipoint,m,j,i)
            print*, " new res", int2_grad1_u12_bimo_t    (ipoint,m,j,i)
            stop
          endif
          err_tot = err_tot + err_loc
          nrm_tot = nrm_tot + dabs(int2_grad1_u12_bimo_t_old(ipoint,m,j,i))
        enddo
      enddo
    enddo
  enddo
  print *, ' absolute accuracy on int2_grad1_u12_bimo_t (%) =', 100.d0 * err_tot / nrm_tot

  return
end

! ---


