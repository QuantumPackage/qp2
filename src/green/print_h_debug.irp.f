program print_h_debug
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end

subroutine routine
 use bitmasks
 implicit none
 integer :: i,j
 integer, allocatable :: H_matrix_degree(:,:)
 double precision, allocatable :: H_matrix_phase(:,:)
 integer :: degree
 integer(bit_kind), allocatable :: keys_tmp(:,:,:)
 allocate(keys_tmp(N_int,2,N_det))
 do i = 1, N_det
  print*,''
  call debug_det(psi_det(1,1,i),N_int)
  do j = 1, N_int
   keys_tmp(j,1,i) = psi_det(j,1,i)
   keys_tmp(j,2,i) = psi_det(j,2,i)
  enddo
 enddo
 if(N_det.gt.10000)then
  print*,'Warning !!!'
  print*,'Number of determinants is ',N_det
  print*,'It means that the H matrix will be enormous !'
  print*,'stoppping ..'
  stop
 endif
 print*,''
 print*,'Determinants '
 do i = 1, N_det
 enddo
 allocate(H_matrix_degree(N_det,N_det),H_matrix_phase(N_det,N_det))
 integer         :: exc(0:2,2,2)
 double precision  :: phase
 do i = 1, N_det
  do j = i, N_det 
   call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
   H_matrix_degree(i,j) = degree
   H_matrix_degree(j,i) = degree
   phase = 0.d0
   if(degree==1.or.degree==2)then
    call get_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,degree,phase,N_int)
   endif
   H_matrix_phase(i,j) = phase
   H_matrix_phase(j,i) = phase
  enddo
 enddo
 print*,'H matrix '
 double precision :: s2
 complex*16 :: ref_h_matrix
 ref_h_matrix = h_matrix_all_dets_complex(1,1)
 print*,'HF like determinant energy = ',ref_bitmask_energy+nuclear_repulsion
 print*,'Ref element of H_matrix    = ',ref_h_matrix+nuclear_repulsion
 print*,'Printing the H matrix ...'
 print*,''
 print*,''
!do i = 1, N_det
! H_matrix_all_dets(i,i) -= ref_h_matrix
!enddo

 do i = 1, N_det
  H_matrix_all_dets_complex(i,i) += nuclear_repulsion
 enddo

!do i = 5,N_det
! H_matrix_all_dets(i,3) = 0.d0
! H_matrix_all_dets(3,i) = 0.d0
! H_matrix_all_dets(i,4) = 0.d0
! H_matrix_all_dets(4,i) = 0.d0
!enddo




! TODO: change for complex
 do i = 1, N_det
  write(*,'(I3,X,A3,2000(E24.15))')i,' | ',H_matrix_all_dets_complex(i,:)
 enddo

! print*,''
! print*,''
! print*,''
! print*,'Printing the degree of excitations within the H matrix'
! print*,''
! print*,''
! do i = 1, N_det
!  write(*,'(I3,X,A3,X,1000(I1,X))')i,' | ',H_matrix_degree(i,:)
! enddo
!
!
! print*,''
! print*,''
! print*,'Printing the phase of the Hamiltonian matrix elements '
! print*,''
! print*,''
! do i = 1, N_det
!  write(*,'(I3,X,A3,X,1000(F3.0,X))')i,' | ',H_matrix_phase(i,:)
! enddo
! print*,''


! double precision, allocatable  :: eigenvalues(:)
! complex*16, allocatable  :: eigenvectors(:,:)
! double precision, allocatable  :: s2_eigvalues(:)
! allocate (eigenvectors(size(H_matrix_all_dets,1),N_det))
! allocate (eigenvalues(N_det),s2_eigvalues(N_det))
! call lapack_diag_complex(eigenvalues,eigenvectors,                       &
!     H_matrix_all_dets,size(H_matrix_all_dets,1),N_det)
! print*,'Two first eigenvectors '
! call u_0_S2_u_0(s2_eigvalues,eigenvectors,n_det,keys_tmp,N_int,N_det,size(eigenvectors,1))
! do j =1, N_states
!   print*,'s2 = ',s2_eigvalues(j)
!   print*,'e  = ',eigenvalues(j)
!   print*,'coefs : '
!   do i = 1, N_det
!    print*,'i = ',i,eigenvectors(i,j)
!   enddo
!   if(j>1)then
!    print*,'Delta E(H)  = ',eigenvalues(1) - eigenvalues(j)
!    print*,'Delta E(eV) = ',(eigenvalues(1) - eigenvalues(j))*27.2114d0
!   endif
! enddo
! complex*16               :: get_mo_bielec_integral,k_a_iv,k_b_iv
! integer :: h1,p1,h2,p2
! h1 = 10
! p1 = 16 
! h2 = 14 
! p2 = 14
!!h1 = 1
!!p1 = 4
!!h2 = 2
!!p2 = 2
! k_a_iv = get_mo_bielec_integral(h1,h2,p2,p1,mo_integrals_map)
! h2 = 15 
! p2 = 15
! k_b_iv = get_mo_bielec_integral(h1,h2,p2,p1,mo_integrals_map)
! print*,'k_a_iv = ',k_a_iv
! print*,'k_b_iv = ',k_b_iv
! complex*16 :: k_av,k_bv,k_ai,k_bi
! h1 = 16
! p1 = 14 
! h2 = 14 
! p2 = 16
! k_av = get_mo_bielec_integral(h1,h2,p1,p2,mo_integrals_map)
! h1 = 16
! p1 = 15 
! h2 = 15 
! p2 = 16
! k_bv = get_mo_bielec_integral(h1,h2,p1,p2,mo_integrals_map)
! 
! h1 = 10
! p1 = 14 
! h2 = 14 
! p2 = 10
! k_ai = get_mo_bielec_integral(h1,h2,p1,p2,mo_integrals_map)
! 
! h1 = 10
! p1 = 15 
! h2 = 15 
! p2 = 10
! k_bi = get_mo_bielec_integral(h1,h2,p1,p2,mo_integrals_map)
! 
! print*,'k_av, k_bv = ',k_av,k_bv
! print*,'k_ai, k_bi = ',k_ai,k_bi
! complex*16 :: k_iv
!
! h1 = 10
! p1 = 16 
! h2 = 16 
! p2 = 10
! k_iv = get_mo_bielec_integral(h1,h2,p1,p2,mo_integrals_map)
! print*,'k_iv       = ',k_iv
end
