
! ---

 BEGIN_PROVIDER [ double precision, natorb_tc_reigvec_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, natorb_tc_leigvec_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, natorb_tc_eigval, (mo_num)]

  BEGIN_DOC
  !
  ! natorb_tc_reigvec_mo : RIGHT eigenvectors of the ground state transition matrix (equivalent of natural orbitals)
  ! natorb_tc_leigvec_mo : LEFT  eigenvectors of the ground state transition matrix (equivalent of natural orbitals)
  ! natorb_tc_eigval     : eigenvalues of the ground state transition matrix (equivalent of the occupation numbers). WARNINING :: can be negative !!
  !
  END_DOC

  implicit none
  integer                       :: i, j, k, n_real
  double precision              :: thr_d, thr_nd, thr_deg, accu
  double precision              :: accu_d, accu_nd
  double precision, allocatable :: dm_tmp(:,:), fock_diag(:)

  allocate(dm_tmp(mo_num,mo_num), fock_diag(mo_num))

  dm_tmp(:,:) = -tc_transition_matrix(:,:,1,1)

  print *, ' dm_tmp'
  do i = 1, mo_num
    fock_diag(i) = fock_matrix_tc_mo_tot(i,i)
    write(*, '(100(F16.10,X))') -dm_tmp(:,i)
  enddo

  thr_d   = 1.d-6
  thr_nd  = 1.d-6
  thr_deg = 1.d-3
  call diag_mat_per_fock_degen( fock_diag, dm_tmp, mo_num, thr_d, thr_nd, thr_deg & 
                              , natorb_tc_leigvec_mo, natorb_tc_reigvec_mo, natorb_tc_eigval)
!   call non_hrmt_bieig( mo_num, dm_tmp&
!                      , natorb_tc_leigvec_mo, natorb_tc_reigvec_mo& 
!                      , n_real, natorb_tc_eigval )

  accu = 0.d0
  do i = 1, n_real
    print*,'natorb_tc_eigval(i) = ',-natorb_tc_eigval(i)
    accu += -natorb_tc_eigval(i)
  enddo
  print *, ' accu = ', accu

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

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, mo_num
    accu_d += dm_tmp(i,i)
    !write(*,'(100(F16.10,X))')dm_tmp(:,i)
    do j = 1, mo_num
      if(i==j)cycle
      accu_nd += dabs(dm_tmp(j,i))
    enddo
  enddo
  print *, ' Trace of the overlap between TC natural orbitals     ', accu_d
  print *, ' L1 norm of extra diagonal elements of overlap matrix ', accu_nd

  deallocate(dm_tmp, fock_diag)
 
END_PROVIDER 

! ---
 
 BEGIN_PROVIDER [ double precision, fock_diag_sorted_r_natorb, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, fock_diag_sorted_l_natorb, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, fock_diag_sorted_v_natorb, (mo_num)]

  implicit none
  integer                       :: i,j,k
  integer,          allocatable :: iorder(:)
  double precision, allocatable :: fock_diag(:)

  print *, ' Diagonal elements of the Fock matrix before '

  do i = 1, mo_num
   write(*,*) i, Fock_matrix_tc_mo_tot(i,i)
  enddo

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

  allocate(iorder(mo_num))
  do i = 1, mo_num
   iorder(i) = i
  enddo 
  call dsort(fock_diag, iorder, mo_num)

  print *, ' Diagonal elements of the Fock matrix after '
  do i = 1, mo_num
   write(*,*) i, fock_diag(i)
  enddo
  deallocate(fock_diag)

  do i = 1, mo_num 
    fock_diag_sorted_v_natorb(i) = natorb_tc_eigval(iorder(i))
    do j = 1, mo_num
      fock_diag_sorted_r_natorb(j,i) = natorb_tc_reigvec_mo(j,iorder(i))
      fock_diag_sorted_l_natorb(j,i) = natorb_tc_leigvec_mo(j,iorder(i))
    enddo
  enddo
  deallocate(iorder)
 
END_PROVIDER 
 
! --- 
 
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

