subroutine get_excitation_general(key_i,key_j, Nint,degree_array,holes_array, particles_array,phase)
 use bitmasks
 BEGIN_DOC
! returns the array, for each spin, of holes/particles between key_i and key_j 
!
! with the following convention: a^+_{particle} a_{hole}|key_i> = |key_j>
 END_DOC
  include 'utils/constants.include.F'
 implicit none
 integer, intent(in)            :: Nint
 integer(bit_kind), intent(in)  :: key_j(Nint,2),key_i(Nint,2)
 integer, intent(out)           :: holes_array(100,2),particles_array(100,2),degree_array(2)
 double precision, intent(out)  :: phase
 integer :: ispin,k,i,pos
 integer(bit_kind) :: key_hole, key_particle
 integer(bit_kind) :: xorvec(N_int_max,2)
 holes_array = -1
 particles_array = -1
 degree_array = 0
  do i = 1, N_int
   xorvec(i,1) = xor( key_i(i,1), key_j(i,1))
   xorvec(i,2) = xor( key_i(i,2), key_j(i,2))
   degree_array(1) += popcnt(xorvec(i,1))
   degree_array(2) += popcnt(xorvec(i,2))
  enddo
  degree_array(1) = shiftr(degree_array(1),1)
  degree_array(2) = shiftr(degree_array(2),1)
  
 do ispin = 1, 2
  k = 1
  !!! GETTING THE HOLES 
  do i = 1, N_int
   key_hole = iand(xorvec(i,ispin),key_i(i,ispin))
   do while(key_hole .ne.0_bit_kind)
    pos = trailz(key_hole)
    holes_array(k,ispin) = 1+ bit_kind_size * (i-1) + pos
    key_hole = ibclr(key_hole,pos)
    k += 1
    if(k .gt.100)then
     print*,'WARNING in get_excitation_general'
     print*,'More than a 100-th excitation for spin ',ispin
     print*,'stoping ...'
     stop
    endif
   enddo 
  enddo
 enddo
 do ispin = 1, 2
  k = 1
  !!! GETTING THE PARTICLES
  do i = 1, N_int
   key_particle = iand(xor(key_i(i,ispin),key_j(i,ispin)),key_j(i,ispin))
   do while(key_particle .ne.0_bit_kind)
    pos = trailz(key_particle)
    particles_array(k,ispin) = 1+ bit_kind_size * (i-1) + pos
    key_particle = ibclr(key_particle,pos)
    k += 1
    if(k .gt.100)then
     print*,'WARNING in get_excitation_general '
     print*,'More than a 100-th excitation for spin ',ispin
     print*,'stoping ...'
     stop
    endif
   enddo
  enddo 
 enddo
 integer :: h,p, i_ok
 integer(bit_kind), allocatable :: det_i(:,:),det_ip(:,:)
 integer                        :: exc(0:2,2,2)
 double precision :: phase_tmp
 allocate(det_i(Nint,2),det_ip(N_int,2))
 det_i = key_i
 phase = 1.d0
 do ispin = 1, 2
  do i = 1, degree_array(ispin)
   h = holes_array(i,ispin)
   p = particles_array(i,ispin)
   det_ip = det_i
   call do_single_excitation(det_ip,h,p,ispin,i_ok)
   if(i_ok == -1)then
     print*,'excitation was not possible '
     stop
   endif
   call get_single_excitation(det_i,det_ip,exc,phase_tmp,Nint)
   phase *= phase_tmp
   det_i = det_ip
  enddo
 enddo

end

subroutine get_holes_general(key_i, key_j,Nint, holes_array)
 use bitmasks
 BEGIN_DOC
! returns the array, per spin, of holes between key_i and key_j 
!
! with the following convention: a_{hole}|key_i> --> |key_j>
 END_DOC
 implicit none
 integer, intent(in)            :: Nint
 integer(bit_kind), intent(in)  :: key_j(Nint,2),key_i(Nint,2)
 integer, intent(out)           :: holes_array(100,2)
 integer(bit_kind) :: key_hole
 integer :: ispin,k,i,pos
 holes_array = -1
 do ispin = 1, 2
  k = 1
  do i = 1, N_int
   key_hole = iand(xor(key_i(i,ispin),key_j(i,ispin)),key_i(i,ispin))
   do while(key_hole .ne.0_bit_kind)
    pos = trailz(key_hole)
    holes_array(k,ispin) = 1+ bit_kind_size * (i-1) + pos
    key_hole = ibclr(key_hole,pos)
    k += 1
    if(k .gt.100)then
     print*,'WARNING in get_holes_general'
     print*,'More than a 100-th excitation for spin ',ispin
     print*,'stoping ...'
     stop
    endif
   enddo 
  enddo
 enddo
end

subroutine get_particles_general(key_i, key_j,Nint,particles_array)
 use bitmasks
 BEGIN_DOC
! returns the array, per spin, of particles between key_i and key_j 
!
! with the following convention: a^dagger_{particle}|key_i> --> |key_j>
 END_DOC
 implicit none
 integer, intent(in)            :: Nint
 integer(bit_kind), intent(in)  :: key_j(Nint,2),key_i(Nint,2)
 integer, intent(out)           :: particles_array(100,2)
 integer(bit_kind) :: key_particle
 integer :: ispin,k,i,pos
 particles_array = -1
 do ispin = 1, 2
  k = 1
  do i = 1, N_int
   key_particle = iand(xor(key_i(i,ispin),key_j(i,ispin)),key_j(i,ispin))
   do while(key_particle .ne.0_bit_kind)
    pos = trailz(key_particle)
    particles_array(k,ispin) = 1+ bit_kind_size * (i-1) + pos
    key_particle = ibclr(key_particle,pos)
    k += 1
    if(k .gt.100)then
     print*,'WARNING in get_holes_general'
     print*,'More than a 100-th excitation for spin ',ispin
     print*,'Those are the two determinants'
     call debug_det(key_i, N_int)
     call debug_det(key_j, N_int)
     print*,'stoping ...'
     stop
    endif
   enddo 
  enddo
 enddo
end

subroutine get_phase_general(key_i,Nint,degree, holes_array, particles_array,phase)
 implicit none
 integer, intent(in)            :: degree(2), Nint
 integer(bit_kind), intent(in)  :: key_i(Nint,2)
 integer, intent(in)            :: holes_array(100,2),particles_array(100,2)
 double precision, intent(out)  :: phase
 integer :: i,ispin,h,p, i_ok
 integer(bit_kind), allocatable :: det_i(:,:),det_ip(:,:)
 integer                        :: exc(0:2,2,2)
 double precision :: phase_tmp
 allocate(det_i(Nint,2),det_ip(N_int,2))
 det_i = key_i
 phase = 1.d0
 do ispin = 1, 2
  do i = 1, degree(ispin)
   h = holes_array(i,ispin)
   p = particles_array(i,ispin)
   det_ip = det_i
   call do_single_excitation(det_ip,h,p,ispin,i_ok)
   if(i_ok == -1)then
     print*,'excitation was not possible '
     stop
   endif
   call get_single_excitation(det_i,det_ip,exc,phase_tmp,Nint)
   phase *= phase_tmp
   det_i = det_ip
  enddo
 enddo

end

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

! ---
subroutine triple_htilde_mu_mat_fock_bi_ortho(Nint, key_j, key_i, hmono, htwoe, hthree, htot)
  use bitmasks
  BEGIN_DOC
! <key_j | H_tilde | key_i> for triple excitation  
!!
!! WARNING !!
! 
! Genuine triple excitations of the same spin are not yet implemented
  END_DOC
  implicit none
  integer(bit_kind), intent(in)  :: key_j(N_int,2),key_i(N_int,2)
  integer, intent(in)            :: Nint
  double precision, intent(out)  :: hmono, htwoe, hthree, htot
  integer                        :: occ(N_int*bit_kind_size,2)
  integer                        :: Ne(2),i,j,ii,jj,ispin,jspin,k,kk
  integer                        :: degree,exc_double(0:2,2,2),exc_single(0:2,2,2)
  integer                        :: degree_alpha,degree_beta
  integer                        :: h1, p1, h2, p2, s1, s2, h3, p3, s3, h4, p4, s4
  double precision :: phase_double, phase_single
  integer(bit_kind)  :: key_j_alpha(N_int,2),key_i_alpha(N_int,2)
  integer(bit_kind)  :: key_j_beta(N_int,2),key_i_beta(N_int,2)
  integer :: other_spin(2)

  hmono  = 0.d0
  htwoe  = 0.d0
  hthree = 0.d0
  htot   = 0.d0
  call bitstring_to_list_ab(key_i,occ,Ne,N_int)
  call get_excitation_degree(key_i,key_j,degree,N_int)
  if(degree.ne.3)then
   return
  endif
  other_spin(1) = 2
  other_spin(2) = 1
  do i = 1, N_int
   key_j_alpha(i,1) = key_j(i,1)
   key_j_alpha(i,2) = 0_bit_kind
   key_i_alpha(i,1) = key_i(i,1)
   key_i_alpha(i,2) = 0_bit_kind

   key_j_beta(i,2) = key_j(i,2)
   key_j_beta(i,1) = 0_bit_kind
   key_i_beta(i,2) = key_i(i,2)
   key_i_beta(i,1) = 0_bit_kind
  enddo
  ! check whether it is a triple excitation of the same spin
  
  call get_excitation_degree(key_i_alpha,key_j_alpha,degree_alpha,N_int)
  call get_excitation_degree(key_i_beta,key_j_beta,degree_beta,N_int)
  if(degree_alpha==3.or.degree_beta==3)then
   return
  else 
   if(degree_alpha == 2.and.degree_beta == 1)then ! double alpha + single beta
    call get_double_excitation(key_i_alpha,key_j_alpha,exc_double,phase_double,N_int)
    call decode_exc(exc_double,2,h1,p1,h2,p2,s1,s2)
    call get_single_excitation(key_i_beta,key_j_beta,exc_single,phase_single,N_int)
    call decode_exc(exc_single,1,h3,p3,h4,p4,s3,s4)
   else if(degree_beta == 2 .and. degree_alpha == 1)then ! double beta + single alpha 
    call get_double_excitation(key_i_beta,key_j_beta,exc_double,phase_double,N_int)
    call decode_exc(exc_double,2,h1,p1,h2,p2,s1,s2)
    call get_single_excitation(key_i_alpha,key_j_alpha,exc_single,phase_single,N_int)
    call decode_exc(exc_single,1,h3,p3,h4,p4,s3,s4)
   else 
    print*,'PB !!'
    print*,'degree_beta, degree_alpha',degree_beta, degree_alpha
    print*,'degree',degree
    stop
   endif
   hthree = three_body_ints_bi_ort(p3,p2,p1,h3,h2,h1) - three_body_ints_bi_ort(p3,p2,p1,h3,h1,h2)
   hthree  *= phase_single * phase_double
  endif
  htot = hthree
 end

