subroutine H_tc_s2_u_0_with_pure_three(v_0, s_0, u_0, N_st, sze)
  BEGIN_DOC
  ! Computes $v_0 = H^TC | u_0\rangle$ WITH PURE TRIPLE EXCITATION TERMS 
  !
  ! Assumes that the determinants are in psi_det
  !
  ! istart, iend, ishift, istep are used in ZMQ parallelization.
  END_DOC

  use bitmasks
  implicit none

  integer,          intent(in)  :: N_st,sze
  double precision, intent(in)  :: u_0(sze,N_st)
  double precision, intent(out) :: v_0(sze,N_st), s_0(sze,N_st)
  call H_tc_s2_u_0_opt(v_0, s_0, u_0, N_st, sze)
  integer :: i,j,degree,ist
  double precision :: hmono, htwoe, hthree, htot
  do i = 1, N_det
   do j = 1, N_det
    call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
    if(degree .ne. 3)cycle
    call triple_htilde_mu_mat_fock_bi_ortho(N_int, psi_det(1,1,i), psi_det(1,1,j), hmono, htwoe, hthree, htot)
    do ist = 1, N_st
     v_0(i,ist) += htot * u_0(j,ist)
    enddo
   enddo
  enddo
end

subroutine H_tc_s2_u_0_with_pure_three_omp(v_0, s_0, u_0, N_st, sze)
  BEGIN_DOC
  ! Computes $v_0 = H^TC | u_0\rangle$ WITH PURE TRIPLE EXCITATION TERMS 
  !
  ! Assumes that the determinants are in psi_det
  !
  ! istart, iend, ishift, istep are used in ZMQ parallelization.
  END_DOC

  use bitmasks
  implicit none

  integer,          intent(in)  :: N_st,sze
  double precision, intent(in)  :: u_0(sze,N_st)
  double precision, intent(out) :: v_0(sze,N_st), s_0(sze,N_st)
  call H_tc_s2_u_0_opt(v_0, s_0, u_0, N_st, sze)
  integer :: i,j,degree,ist
  double precision :: hmono, htwoe, hthree, htot
  !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
  !$OMP SHARED(N_st, N_det, N_int, psi_det, u_0, v_0)       &
  !$OMP PRIVATE(ist, i, j, degree, hmono, htwoe, hthree,htot)
  do i = 1, N_det
   do j = 1, N_det
    call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
    if(degree .ne. 3)cycle
    call triple_htilde_mu_mat_fock_bi_ortho(N_int, psi_det(1,1,i), psi_det(1,1,j), hmono, htwoe, hthree, htot)
    do ist = 1, N_st
     v_0(i,ist) += htot * u_0(j,ist)
    enddo
   enddo
  enddo
 !$OMP END PARALLEL DO
end

! ---

subroutine H_tc_s2_dagger_u_0_with_pure_three(v_0, s_0, u_0, N_st, sze)
  BEGIN_DOC
  ! Computes $v_0 = (H^TC)^dagger | u_0\rangle$ WITH PURE TRIPLE EXCITATION TERMS 
  !
  ! Assumes that the determinants are in psi_det
  !
  ! istart, iend, ishift, istep are used in ZMQ parallelization.
  END_DOC

  use bitmasks
  implicit none

  integer,          intent(in)  :: N_st,sze
  double precision, intent(in)  :: u_0(sze,N_st)
  double precision, intent(out) :: v_0(sze,N_st), s_0(sze,N_st)
  call H_tc_s2_dagger_u_0_opt(v_0, s_0, u_0, N_st, sze)
  integer :: i,j,degree,ist
  double precision :: hmono, htwoe, hthree, htot
  do i = 1, N_det
   do j = 1, N_det
    call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
    if(degree .ne. 3)cycle
    call triple_htilde_mu_mat_fock_bi_ortho(N_int, psi_det(1,1,j), psi_det(1,1,i), hmono, htwoe, hthree, htot)
    do ist = 1, N_st
     v_0(i,ist) += htot * u_0(j,ist)
    enddo
   enddo
  enddo
end

subroutine H_tc_s2_dagger_u_0_with_pure_three_omp(v_0, s_0, u_0, N_st, sze)
  BEGIN_DOC
  ! Computes $v_0 = (H^TC)^dagger | u_0\rangle$ WITH PURE TRIPLE EXCITATION TERMS 
  !
  ! Assumes that the determinants are in psi_det
  !
  ! istart, iend, ishift, istep are used in ZMQ parallelization.
  END_DOC

  use bitmasks
  implicit none

  integer,          intent(in)  :: N_st,sze
  double precision, intent(in)  :: u_0(sze,N_st)
  double precision, intent(out) :: v_0(sze,N_st), s_0(sze,N_st)
  call H_tc_s2_dagger_u_0_opt(v_0, s_0, u_0, N_st, sze)
  integer :: i,j,degree,ist
  double precision :: hmono, htwoe, hthree, htot
  !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8) &
  !$OMP SHARED(N_st, N_det, N_int, psi_det, u_0, v_0)       &
  !$OMP PRIVATE(ist, i, j, degree, hmono, htwoe, hthree,htot)
  do i = 1, N_det
   do j = 1, N_det
    call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
    if(degree .ne. 3)cycle
    call triple_htilde_mu_mat_fock_bi_ortho(N_int, psi_det(1,1,j), psi_det(1,1,i), hmono, htwoe, hthree, htot)
    do ist = 1, N_st
     v_0(i,ist) += htot * u_0(j,ist)
    enddo
   enddo
  enddo
 !$OMP END PARALLEL DO
end

! ---
subroutine triple_htilde_mu_mat_fock_bi_ortho(Nint, key_j, key_i, hmono, htwoe, hthree, htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> for triple excitation  
!!
!! WARNING !!
! 
! Genuine triple excitations of the same spin are not yet implemented
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  integer, intent(in)            :: Nint
  double precision, intent(out)  :: hmono, htwoe, hthree, htot
  integer                        :: degree
  integer                        :: h1, p1, h2, p2, s1, s2, h3, p3, s3
  integer                        :: holes_array(100,2),particles_array(100,2),degree_array(2)
  double precision               :: phase,sym_3_e_int_from_6_idx_tensor

  hmono  = 0.d0
  htwoe  = 0.d0
  hthree = 0.d0
  htot   = 0.d0
  call get_excitation_general(key_j, key_i, Nint,degree_array,holes_array, particles_array,phase)
  degree = degree_array(1) + degree_array(2)
  if(degree .ne. 3)return
  if(degree_array(1)==3.or.degree_array(2)==3)then
   if(degree_array(1) == 3)then
    h1 = holes_array(1,1) 
    h2 = holes_array(2,1) 
    h3 = holes_array(3,1) 
    p1 = particles_array(1,1) 
    p2 = particles_array(2,1) 
    p3 = particles_array(3,1) 
   else
    h1 = holes_array(1,2) 
    h2 = holes_array(2,2) 
    h3 = holes_array(3,2) 
    p1 = particles_array(1,2) 
    p2 = particles_array(2,2) 
    p3 = particles_array(3,2) 
   endif
   hthree = sym_3_e_int_from_6_idx_tensor(p3, p2, p1, h3, h2, h1)
  else 
   if(degree_array(1) == 2.and.degree_array(2) == 1)then ! double alpha + single beta
    h1 = holes_array(1,1) 
    h2 = holes_array(2,1) 
    h3 = holes_array(1,2) 
    p1 = particles_array(1,1) 
    p2 = particles_array(2,1) 
    p3 = particles_array(1,2) 
   else if(degree_array(2) == 2 .and. degree_array(1) == 1)then ! double beta + single alpha 
    h1 = holes_array(1,2) 
    h2 = holes_array(2,2) 
    h3 = holes_array(1,1) 
    p1 = particles_array(1,2) 
    p2 = particles_array(2,2) 
    p3 = particles_array(1,1) 
   else 
    print*,'PB !!'
    stop
   endif
   hthree = three_body_ints_bi_ort(p3,p2,p1,h3,h2,h1) - three_body_ints_bi_ort(p3,p2,p1,h3,h1,h2)
  endif
  hthree  *= phase
  htot = hthree
 end

