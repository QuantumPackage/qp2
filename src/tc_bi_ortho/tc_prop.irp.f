
BEGIN_PROVIDER [ double precision, tc_transition_matrix, (mo_num, mo_num,N_states,N_states) ]
 implicit none
 BEGIN_DOC
 ! tc_transition_matrix(p,h,istate,jstate) = <Chi_istate| a^\dagger_p a_h |Phi_jstate>
 !
 ! where <Chi_istate| and |Phi_jstate> are the left/right eigenvectors on a bi-ortho basis
 END_DOC
 integer :: i,j,istate,jstate,m,n,p,h
 double precision :: phase
 integer, allocatable           :: occ(:,:)
 integer                        :: n_occ_ab(2),degree,exc(0:2,2,2)
 allocate(occ(N_int*bit_kind_size,2))
 tc_transition_matrix = 0.d0
 do istate = 1, N_states
  do jstate = 1, N_states
   do i = 1, N_det
    do j = 1, N_det
     call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
     if(degree.gt.1)then
      cycle
     else if (degree == 0)then
      call bitstring_to_list_ab(psi_det(1,1,i), occ, n_occ_ab, N_int)
      do p = 1, n_occ_ab(1) ! browsing the alpha electrons
       m = occ(p,1)
       tc_transition_matrix(m,m,istate,jstate)+= psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,jstate)
      enddo
      do p = 1, n_occ_ab(2) ! browsing the beta electrons
       m = occ(p,1)
       tc_transition_matrix(m,m,istate,jstate)+= psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,jstate)
      enddo
     else
      call get_single_excitation(psi_det(1,1,j),psi_det(1,1,i),exc,phase,N_int)
      if (exc(0,1,1) == 1) then
        ! Single alpha
        h = exc(1,1,1) ! hole in psi_det(1,1,j) 
        p = exc(1,2,1) ! particle in psi_det(1,1,j) 
      else
        ! Single beta
        h = exc(1,1,2) ! hole in psi_det(1,1,j) 
        p = exc(1,2,2) ! particle in psi_det(1,1,j) 
      endif
      tc_transition_matrix(p,h,istate,jstate)+= phase * psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,jstate)
     endif
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER


  BEGIN_PROVIDER [ double precision, natorb_tc_reigvec_mo, (mo_num, mo_num)]
 &BEGIN_PROVIDER [ double precision, natorb_tc_leigvec_mo, (mo_num, mo_num)]
 &BEGIN_PROVIDER [ double precision, natorb_tc_eigval, (mo_num)]
  implicit none
 BEGIN_DOC
 ! natorb_tc_reigvec_mo : RIGHT eigenvectors of the ground state transition matrix (equivalent of natural orbitals)
 ! natorb_tc_leigvec_mo : LEFT  eigenvectors of the ground state transition matrix (equivalent of natural orbitals)
 ! natorb_tc_eigval     : eigenvalues of the ground state transition matrix (equivalent of the occupation numbers). WARNINING :: can be negative !!
 END_DOC
  double precision, allocatable :: dm_tmp(:,:)
  integer :: i,j,k,n_real
  allocate( dm_tmp(mo_num,mo_num))
  dm_tmp(:,:) = -tc_transition_matrix(:,:,1,1)
  print*,'dm_tmp'
  do i = 1, mo_num
   write(*,'(100(F16.10,X))')-dm_tmp(:,i)
  enddo
!    call non_hrmt_diag_split_degen( mo_num, dm_tmp&
    call non_hrmt_fock_mat( mo_num, dm_tmp&
!   call non_hrmt_bieig( mo_num, dm_tmp&
                      , natorb_tc_leigvec_mo, natorb_tc_reigvec_mo& 
                      , n_real, natorb_tc_eigval )
 double precision :: accu
  accu = 0.d0
  do i = 1, n_real
   print*,'natorb_tc_eigval(i) = ',-natorb_tc_eigval(i)
   accu += -natorb_tc_eigval(i)
  enddo
  print*,'accu = ',accu
  dm_tmp = 0.d0
  do i = 1, n_real
   accu = 0.d0
   do k = 1, mo_num
    accu += natorb_tc_reigvec_mo(k,i) * natorb_tc_leigvec_mo(k,i)
   enddo
   accu = 1.d0/dsqrt(dabs(accu))
   natorb_tc_reigvec_mo(:,i) *= accu
   natorb_tc_leigvec_mo(:,i) *= accu
   do j = 1, n_real
    do k = 1, mo_num
     dm_tmp(j,i) += natorb_tc_reigvec_mo(k,i) * natorb_tc_leigvec_mo(k,j)
    enddo
   enddo
  enddo
  double precision :: accu_d, accu_nd
  accu_d = 0.d0
  accu_nd = 0.d0
  do i = 1, mo_num
   accu_d += dm_tmp(i,i)
 !  write(*,'(100(F16.10,X))')dm_tmp(:,i)
   do j = 1, mo_num
    if(i==j)cycle
    accu_nd += dabs(dm_tmp(j,i))
   enddo
  enddo
  print*,'Trace of the overlap between TC natural orbitals     ',accu_d
  print*,'L1 norm of extra diagonal elements of overlap matrix ',accu_nd
 
 
 END_PROVIDER 
 
  BEGIN_PROVIDER [ double precision, fock_diag_sorted_r_natorb, (mo_num, mo_num)]
 &BEGIN_PROVIDER [ double precision, fock_diag_sorted_l_natorb, (mo_num, mo_num)]
 &BEGIN_PROVIDER [ double precision, fock_diag_sorted_v_natorb, (mo_num)]
  implicit none
  integer ::i,j,k
  print*,'Diagonal elements of the Fock matrix before '
  do i = 1, mo_num
   write(*,*)i,Fock_matrix_tc_mo_tot(i,i)
  enddo
  double precision, allocatable :: fock_diag(:)
  allocate(fock_diag(mo_num))
  fock_diag = 0.d0
  do i = 1, mo_num
   fock_diag(i) = 0.d0
   do j = 1, mo_num
    do k = 1, mo_num
     fock_diag(i) += natorb_tc_leigvec_mo(k,i) * Fock_matrix_tc_mo_tot(k,j) * natorb_tc_reigvec_mo(j,i) 
    enddo
   enddo
  enddo
  integer, allocatable :: iorder(:)
  allocate(iorder(mo_num))
  do i = 1, mo_num
   iorder(i) = i
  enddo 
  call dsort(fock_diag,iorder,mo_num)
  print*,'Diagonal elements of the Fock matrix after '
  do i = 1, mo_num
   write(*,*)i,fock_diag(i)
  enddo
  do i = 1, mo_num 
   fock_diag_sorted_v_natorb(i) = natorb_tc_eigval(iorder(i))
   do j = 1, mo_num
    fock_diag_sorted_r_natorb(j,i) = natorb_tc_reigvec_mo(j,iorder(i))
    fock_diag_sorted_l_natorb(j,i) = natorb_tc_leigvec_mo(j,iorder(i))
   enddo
  enddo
 
 END_PROVIDER 
 
 
 
  BEGIN_PROVIDER [ double precision, natorb_tc_reigvec_ao, (ao_num, mo_num)]
 &BEGIN_PROVIDER [ double precision, natorb_tc_leigvec_ao, (ao_num, mo_num)]
 &BEGIN_PROVIDER [ double precision, overlap_natorb_tc_eigvec_ao, (mo_num, mo_num) ]
 
   BEGIN_DOC
   ! EIGENVECTORS OF FOCK MATRIX ON THE AO BASIS and their OVERLAP
   !
   ! THE OVERLAP SHOULD BE THE SAME AS overlap_natorb_tc_eigvec_mo
   END_DOC
 
   implicit none
   integer                       :: i, j, k, q, p
   double precision              :: accu, accu_d
   double precision, allocatable :: tmp(:,:)
 
 
 !  ! MO_R x R
   call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0          &
             , mo_r_coef, size(mo_r_coef, 1)                   &
             , fock_diag_sorted_r_natorb, size(fock_diag_sorted_r_natorb, 1) &
             , 0.d0, natorb_tc_reigvec_ao, size(natorb_tc_reigvec_ao, 1) )
 !
   ! MO_L x L
   call dgemm( 'N', 'N', ao_num, mo_num, mo_num, 1.d0          &
             , mo_l_coef, size(mo_l_coef, 1)                   &
             , fock_diag_sorted_l_natorb, size(fock_diag_sorted_l_natorb, 1) &
             , 0.d0, natorb_tc_leigvec_ao, size(natorb_tc_leigvec_ao, 1) )
 
 
   allocate( tmp(mo_num,ao_num) )
 
   ! tmp <-- L.T x S_ao
   call dgemm( "T", "N", mo_num, ao_num, ao_num, 1.d0                                           &
             , natorb_tc_leigvec_ao, size(natorb_tc_leigvec_ao, 1), ao_overlap, size(ao_overlap, 1) &
             , 0.d0, tmp, size(tmp, 1) )
 
   ! S <-- tmp x R
   call dgemm( "N", "N", mo_num, mo_num, ao_num, 1.d0                             &
             , tmp, size(tmp, 1), natorb_tc_reigvec_ao, size(natorb_tc_reigvec_ao, 1) &
             , 0.d0, overlap_natorb_tc_eigvec_ao, size(overlap_natorb_tc_eigvec_ao, 1) )
 
   deallocate( tmp )
 
   ! ---
   double precision :: norm
   do i = 1, mo_num
    norm = 1.d0/dsqrt(dabs(overlap_natorb_tc_eigvec_ao(i,i)))
    do j = 1, mo_num
     natorb_tc_reigvec_ao(j,i) *= norm
     natorb_tc_leigvec_ao(j,i) *= norm
    enddo
   enddo
 
   allocate( tmp(mo_num,ao_num) )
 
   ! tmp <-- L.T x S_ao
   call dgemm( "T", "N", mo_num, ao_num, ao_num, 1.d0                                           &
             , natorb_tc_leigvec_ao, size(natorb_tc_leigvec_ao, 1), ao_overlap, size(ao_overlap, 1) &
             , 0.d0, tmp, size(tmp, 1) )
 
   ! S <-- tmp x R
   call dgemm( "N", "N", mo_num, mo_num, ao_num, 1.d0                             &
             , tmp, size(tmp, 1), natorb_tc_reigvec_ao, size(natorb_tc_reigvec_ao, 1) &
             , 0.d0, overlap_natorb_tc_eigvec_ao, size(overlap_natorb_tc_eigvec_ao, 1) )
 
 
 
   deallocate( tmp )
 
   accu_d = 0.d0
   accu = 0.d0
   do i = 1, mo_num
     accu_d += overlap_natorb_tc_eigvec_ao(i,i)
     do j = 1, mo_num
       if(i==j)cycle
       accu += dabs(overlap_natorb_tc_eigvec_ao(j,i))
     enddo
   enddo
   print*,'Trace of the overlap_natorb_tc_eigvec_ao           = ',accu_d
   print*,'mo_num                                             = ',mo_num
   print*,'L1 norm of extra diagonal elements of overlap matrix ',accu
   accu = accu / dble(mo_num**2)
 
 END_PROVIDER

 BEGIN_PROVIDER [double precision, tc_bi_ortho_dipole, (3,N_states)]
 implicit none
 integer :: i,j,istate,m
 double precision :: nuclei_part(3)
 tc_bi_ortho_dipole = 0.d0
 do istate = 1, N_states
  do i = 1, mo_num
   do j = 1, mo_num
    tc_bi_ortho_dipole(1,istate) += -(tc_transition_matrix(j,i,istate,istate)) *  mo_bi_orth_bipole_x(j,i)
    tc_bi_ortho_dipole(2,istate) += -(tc_transition_matrix(j,i,istate,istate)) *  mo_bi_orth_bipole_y(j,i)
    tc_bi_ortho_dipole(3,istate) += -(tc_transition_matrix(j,i,istate,istate)) *  mo_bi_orth_bipole_z(j,i)
   enddo
  enddo
 enddo

 nuclei_part = 0.d0
 do m = 1, 3
  do i = 1,nucl_num
   nuclei_part(m) += nucl_charge(i) * nucl_coord(i,m)
  enddo
 enddo
!
 do istate = 1, N_states
  do m = 1, 3
    tc_bi_ortho_dipole(m,istate) += nuclei_part(m)
  enddo
 enddo
 END_PROVIDER

