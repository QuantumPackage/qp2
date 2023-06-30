program tc_bi_ortho
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

! call test_h_u0
! call test_slater_tc_opt
! call timing_tot
! call timing_diag
! call timing_single
! call timing_double

  call test_no()
  !call test_no_aba()
  !call test_no_aab()
  !call test_no_aaa()
end

subroutine test_h_u0
 implicit none
 double precision, allocatable :: v_0_ref(:),v_0_new(:),u_0(:), v_0_ref_dagger(:)
 double precision :: accu 
 logical :: do_right
 integer :: i
 allocate(v_0_new(N_det),v_0_ref(N_det),u_0(N_det),v_0_ref_dagger(N_det))
 do_right = .True.
 do i = 1, N_det
  u_0(i) = psi_r_coef_bi_ortho(i,1)
 enddo
 call H_tc_u_0_nstates_openmp(v_0_new,u_0,N_states,N_det, do_right)
 call htc_bi_ortho_calc_tdav_slow (v_0_ref,u_0,N_states,N_det)
 print*,'difference right '
 accu = 0.d0
 do i = 1, N_det
  print*,dabs(v_0_new(i) - v_0_ref(i)),v_0_new(i) , v_0_ref(i)
  accu += dabs(v_0_new(i) - v_0_ref(i))
 enddo
 print*,'accu = ',accu
 do_right = .False.
 v_0_new = 0.d0
 call H_tc_u_0_nstates_openmp(v_0_new,u_0,N_states,N_det, do_right)
 call htcdag_bi_ortho_calc_tdav_slow(v_0_ref_dagger,u_0,N_states,N_det, do_right)
 print*,'difference left'
 accu = 0.d0
 do i = 1, N_det
  print*,dabs(v_0_new(i) - v_0_ref_dagger(i)),v_0_new(i) , v_0_ref_dagger(i)
  accu += dabs(v_0_new(i) - v_0_ref_dagger(i))
 enddo
 print*,'accu = ',accu
end

subroutine test_slater_tc_opt
 implicit none
 integer :: i,j,degree
 double precision :: hmono, htwoe, htot, hthree 
 double precision :: hnewmono, hnewtwoe, hnewthree, hnewtot
 double precision :: accu_d ,i_count, accu
 accu = 0.d0
 accu_d = 0.d0
 i_count = 0.d0
 do i = 1, N_det
  do j = 1,N_det
   call htilde_mu_mat_bi_ortho_slow(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hnewmono, hnewtwoe, hnewthree, hnewtot)
   if(dabs(htot).gt.1.d-15)then
     i_count += 1.D0
     accu += dabs(htot-hnewtot) 
     if(dabs(htot-hnewtot).gt.1.d-8.or.dabs(htot-hnewtot).gt.dabs(htot))then
      call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
      print*,j,i,degree
      call debug_det(psi_det(1,1,i),N_int)
      call debug_det(psi_det(1,1,j),N_int)
      print*,htot,hnewtot,dabs(htot-hnewtot)
      print*,hthree,hnewthree,dabs(hthree-hnewthree)
      stop
     endif
   endif
  enddo
 enddo
 print*,'accu   = ',accu/i_count

end

subroutine timing_tot
 implicit none
 integer :: i,j
 double precision :: wall0, wall1
 double precision, allocatable :: mat_old(:,:),mat_new(:,:)
 double precision :: hmono, htwoe, hthree, htot, i_count
 integer :: degree 
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,2), N_int, hmono, htwoe, hthree, htot)
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,2), N_int, hmono, htwoe, hthree, htot)
 call wall_time(wall0)
 i_count = 0.d0
 do i = 1, N_det
  do j = 1, N_det
!   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   i_count += 1.d0
   call htilde_mu_mat_bi_ortho_slow(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for old hij for total   = ',wall1 - wall0

 call wall_time(wall0)
 i_count = 0.d0
 do i = 1, N_det
  do j = 1, N_det
!   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   i_count += 1.d0
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for new hij for total   = ',wall1 - wall0
 call i_H_j(psi_det(1,1,1), psi_det(1,1,2),N_int,htot)
 call wall_time(wall0)
 i_count = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   call i_H_j(psi_det(1,1,j), psi_det(1,1,i),N_int,htot)
   i_count += 1.d0
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for new hij STANDARD    = ',wall1 - wall0

end

subroutine timing_diag
 implicit none
 integer :: i,j
 double precision :: wall0, wall1
 double precision, allocatable :: mat_old(:,:),mat_new(:,:)
 double precision :: hmono, htwoe, hthree, htot, i_count
 integer :: degree 
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,1), N_int, hmono, htwoe, hthree, htot)
 call wall_time(wall0)
 i_count = 0.d0
 do i = 1, N_det
  do j = i,i 
   i_count += 1.d0
   call htilde_mu_mat_bi_ortho_slow(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for old hij for diagonal= ',wall1 - wall0

 call wall_time(wall0)
 i_count = 0.d0
 do i = 1, N_det
  do j = i,i
   i_count += 1.d0
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for new hij for diagonal= ',wall1 - wall0

end

subroutine timing_single
 implicit none
 integer :: i,j
 double precision :: wall0, wall1,accu
 double precision, allocatable :: mat_old(:,:),mat_new(:,:)
 double precision :: hmono, htwoe, hthree, htot, i_count
 integer :: degree 
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,1), N_int, hmono, htwoe, hthree, htot)
 i_count = 0.d0
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   if(degree.ne.1)cycle
   i_count += 1.d0
   call wall_time(wall0)
   call htilde_mu_mat_bi_ortho_slow(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call wall_time(wall1)
   accu += wall1 - wall0
  enddo
 enddo
 print*,'i_count = ',i_count
 print*,'time for old hij for singles = ',accu

 i_count = 0.d0
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   if(degree.ne.1)cycle
   i_count += 1.d0
   call wall_time(wall0)
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call wall_time(wall1)
   accu += wall1 - wall0
  enddo
 enddo
 print*,'i_count = ',i_count
 print*,'time for new hij for singles = ',accu

end

subroutine timing_double
 implicit none
 integer :: i,j
 double precision :: wall0, wall1,accu
 double precision, allocatable :: mat_old(:,:),mat_new(:,:)
 double precision :: hmono, htwoe, hthree, htot, i_count
 integer :: degree 
 call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,1), psi_det(1,1,1), N_int, hmono, htwoe, hthree, htot)
 i_count = 0.d0
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   if(degree.ne.2)cycle
   i_count += 1.d0
   call wall_time(wall0)
   call htilde_mu_mat_bi_ortho_slow(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call wall_time(wall1)
   accu += wall1 - wall0
  enddo
 enddo
 print*,'i_count = ',i_count
 print*,'time for old hij for doubles = ',accu

 i_count = 0.d0
 accu = 0.d0
 do i = 1, N_det
  do j = 1, N_det
   call get_excitation_degree(psi_det(1,1,j), psi_det(1,1,i),degree,N_int)
   if(degree.ne.2)cycle
   i_count += 1.d0
   call wall_time(wall0)
   call htilde_mu_mat_opt_bi_ortho(psi_det(1,1,j), psi_det(1,1,i), N_int, hmono, htwoe, hthree, htot)
   call wall_time(wall1)
   accu += wall1 - wall0
  enddo
 enddo
 call wall_time(wall1)
 print*,'i_count = ',i_count
 print*,'time for new hij for doubles = ',accu

end

! ---

subroutine test_no()

  implicit none
  integer          :: i, j, k, l
  double precision :: accu, contrib, new, ref, thr

  print*, ' testing normal_two_body_bi_orth ...'

  thr = 1d-8

  PROVIDE normal_two_body_bi_orth_old
  PROVIDE normal_two_body_bi_orth

  accu = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = normal_two_body_bi_orth    (l,k,j,i)
          ref = normal_two_body_bi_orth_old(l,k,j,i)
          contrib = dabs(new - ref)
          accu += contrib
          if(contrib .gt. thr) then
            print*, ' problem on normal_two_body_bi_orth'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

        enddo
      enddo
    enddo
  enddo
  print*, ' accu on normal_two_body_bi_orth = ', accu / dble(mo_num)**4

 return
end

! ---

subroutine test_no_aba()

  implicit none
  integer          :: i, j, k, l
  double precision :: accu, contrib, new, ref, thr

  print*, ' testing no_aba_contraction ...'

  thr = 1d-8

  PROVIDE no_aba_contraction_v0
  PROVIDE no_aba_contraction

  accu = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = no_aba_contraction   (l,k,j,i)
          ref = no_aba_contraction_v0(l,k,j,i)
          contrib = dabs(new - ref)
          accu += contrib
          if(contrib .gt. thr) then
            print*, ' problem on no_aba_contraction'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

        enddo
      enddo
    enddo
  enddo
  print*, ' accu on no_aba_contraction = ', accu / dble(mo_num)**4

 return
end

! ---


subroutine test_no_aab()

  implicit none
  integer          :: i, j, k, l
  double precision :: accu, contrib, new, ref, thr

  print*, ' testing no_aab_contraction ...'

  thr = 1d-8

  PROVIDE no_aab_contraction_v0
  PROVIDE no_aab_contraction

  accu = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = no_aab_contraction   (l,k,j,i)
          ref = no_aab_contraction_v0(l,k,j,i)
          contrib = dabs(new - ref)
          accu += contrib
          if(contrib .gt. thr) then
            print*, ' problem on no_aab_contraction'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

        enddo
      enddo
    enddo
  enddo
  print*, ' accu on no_aab_contraction = ', accu / dble(mo_num)**4

 return
end

! ---

subroutine test_no_aaa()

  implicit none
  integer          :: i, j, k, l
  double precision :: accu, contrib, new, ref, thr

  print*, ' testing no_aaa_contraction ...'

  thr = 1d-8

  PROVIDE no_aaa_contraction_v0
  PROVIDE no_aaa_contraction

  accu = 0.d0
  do i = 1, mo_num
    do j = 1, mo_num
      do k = 1, mo_num
        do l = 1, mo_num

          new = no_aaa_contraction   (l,k,j,i)
          ref = no_aaa_contraction_v0(l,k,j,i)
          contrib = dabs(new - ref)
          accu += contrib
          if(contrib .gt. thr) then
            print*, ' problem on no_aaa_contraction'
            print*, l, k, j, i
            print*, ref, new, contrib
            stop
          endif

        enddo
      enddo
    enddo
  enddo
  print*, ' accu on no_aaa_contraction = ', accu / dble(mo_num)**4

 return
end

! ---
