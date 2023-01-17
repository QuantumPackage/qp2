subroutine non_hrmt_diag_split_degen_bi_orthog(n, A, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  !
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  !
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)
  double precision, allocatable :: reigvec_tmp(:,:), leigvec_tmp(:,:)

  integer                       :: i, j, n_degen,k , iteration
  double precision              :: shift_current
  double precision              :: r,thr,accu_d, accu_nd
  integer,          allocatable :: iorder_origin(:),iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:),S(:,:)
  double precision, allocatable :: Aw(:,:),diag_elem(:),A_save(:,:)
  double precision, allocatable :: im_part(:),re_part(:)
  double precision :: accu,thr_cut, thr_norm=1d0


  thr_cut = 1.d-15
  print*,'Computing the left/right eigenvectors ...'
  print*,'Using the degeneracy splitting algorithm'
 ! initialization
  shift_current = 1.d-15
  iteration = 0 
  print*,'***** iteration = ',iteration


  ! pre-processing the matrix :: sorting by diagonal elements
  allocate(reigvec_tmp(n,n), leigvec_tmp(n,n))
  allocate(diag_elem(n),iorder_origin(n),A_save(n,n))
!  print*,'Aw'
  do i = 1, n
   iorder_origin(i) = i
   diag_elem(i) = A(i,i)
!   write(*,'(100(F16.10,X))')A(:,i)
  enddo
  call dsort(diag_elem, iorder_origin, n)
  do i = 1, n
   do j = 1, n
    A_save(j,i) = A(iorder_origin(j),iorder_origin(i))
   enddo
  enddo

  allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
  allocate(im_part(n),iorder(n))
  allocate( S(n,n) )


  Aw = A_save
  call cancel_small_elmts(aw,n,thr_cut)
  call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
  do i = 1, n
   im_part(i) = -dabs(WI(i))
   iorder(i) = i
  enddo
  call dsort(im_part, iorder, n)
  n_real_eigv = 0
  do i = 1, n
    if(dabs(WI(i)).lt.1.d-20)then
      n_real_eigv += 1
    else
!      print*,'Found an imaginary component to eigenvalue'
!      print*,'Re(i) + Im(i)',WR(i),WI(i)
    endif
  enddo
  if(n_real_eigv.ne.n)then
   shift_current = max(10.d0 * dabs(im_part(1)),shift_current*10.d0)
   print*,'Largest imaginary part found in eigenvalues = ',im_part(1)
   print*,'Splitting the degeneracies by ',shift_current
  else
   print*,'All eigenvalues are real !'
  endif


  do while(n_real_eigv.ne.n)
   iteration += 1
   print*,'***** iteration = ',iteration
   if(shift_current.gt.1.d-3)then
    print*,'shift_current > 1.d-3 !!'
    print*,'Your matrix intrinsically contains complex eigenvalues'
    stop
   endif
   Aw = A_save
   call cancel_small_elmts(Aw,n,thr_cut)
   call split_matrix_degen(Aw,n,shift_current)
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
   n_real_eigv = 0
   do i = 1, n
     if(dabs(WI(i)).lt.1.d-20)then
       n_real_eigv+= 1
     else
!       print*,'Found an imaginary component to eigenvalue'
!       print*,'Re(i) + Im(i)',WR(i),WI(i)
     endif
   enddo
   if(n_real_eigv.ne.n)then
    do i = 1, n
     im_part(i) = -dabs(WI(i))
     iorder(i) = i
    enddo
    call dsort(im_part, iorder, n)
    shift_current = max(10.d0 * dabs(im_part(1)),shift_current*10.d0)
    print*,'Largest imaginary part found in eigenvalues = ',im_part(1)
    print*,'Splitting the degeneracies by ',shift_current
   else
    print*,'All eigenvalues are real !'
   endif
  enddo
  !!!!!!!!!!!!!!!! SORTING THE EIGENVALUES 
  do i = 1, n
   eigval(i) = WR(i)
   iorder(i) = i
  enddo
  call dsort(eigval,iorder,n)
  do i = 1, n
!   print*,'eigval(i) = ',eigval(i)
   reigvec_tmp(:,i) = VR(:,iorder(i))
   leigvec_tmp(:,i) = Vl(:,iorder(i))
  enddo

!!! ONCE ALL EIGENVALUES ARE REAL ::: CHECK BI-ORTHONORMALITY
  !                               check bi-orthogonality
  call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
  print *, ' accu_nd bi-orthog = ', accu_nd
  if(accu_nd .lt. thresh_biorthog_nondiag) then
    print *, ' bi-orthogonality: ok'
  else
    print *, ' '
    print *, ' bi-orthogonality: not imposed yet'
    print *, ' '
    print *, ' '
    print *, ' orthog between degen eigenvect' 
    print *, ' '
    double precision, allocatable :: S_nh_inv_half(:,:)
    allocate(S_nh_inv_half(n,n))
    logical :: complex_root
    deallocate(S_nh_inv_half)
    call impose_orthog_degen_eigvec(n, eigval, reigvec_tmp)
    call impose_orthog_degen_eigvec(n, eigval, leigvec_tmp)
    call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
    if(accu_nd .lt. thresh_biorthog_nondiag) then
      print *, ' bi-orthogonality: ok'
    else 
     print*,'New vectors not bi-orthonormals at ',accu_nd
     call impose_biorthog_qr(n, n, leigvec_tmp, reigvec_tmp, S)
     call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
     if(accu_nd .lt. thresh_biorthog_nondiag) then
       print *, ' bi-orthogonality: ok'
     else 
      print*,'New vectors not bi-orthonormals at ',accu_nd
      print*,'Must be a deep problem ...'
      stop
     endif
    endif
  endif
 
  !! EIGENVECTORS SORTED AND BI-ORTHONORMAL
  do i = 1, n
   do j = 1, n
    VR(iorder_origin(j),i) = reigvec_tmp(j,i)
    VL(iorder_origin(j),i) = leigvec_tmp(j,i)
   enddo
  enddo

  !! RECOMPUTING THE EIGENVALUES 
  eigval = 0.d0
  do i = 1, n
   iorder(i) = i
   accu = 0.d0
   do j = 1, n
    accu += VL(j,i) * VR(j,i) 
    do k = 1, n
     eigval(i) +=  VL(j,i) * A(j,k) * VR(k,i) 
    enddo
   enddo
   eigval(i) *= 1.d0/accu
!   print*,'eigval(i) = ',eigval(i)
  enddo
  !! RESORT JUST TO BE SURE
  call dsort(eigval, iorder, n)
  do i = 1, n
   do j = 1, n
    reigvec(j,i) = VR(j,iorder(i))
    leigvec(j,i) = VL(j,iorder(i))
   enddo
  enddo
  print*,'Checking for final reigvec/leigvec'
  shift_current = max(1.d-10,shift_current)
  print*,'Thr for eigenvectors = ',shift_current
  call check_EIGVEC(n, n, A, eigval, leigvec, reigvec, shift_current, thr_norm, .false.)
  call check_biorthog(n, n, leigvec, reigvec, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
  print *, ' accu_nd bi-orthog = ', accu_nd
  
  if(accu_nd .lt. thresh_biorthog_nondiag) then
    print *, ' bi-orthogonality: ok'
  else 
   print*,'Something went wrong in non_hrmt_diag_split_degen_bi_orthog'
   print*,'Eigenvectors are not bi orthonormal ..'
   print*,'accu_nd = ',accu_nd
   stop
  endif

end 



subroutine non_hrmt_diag_split_degen_s_inv_half(n, A, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  !
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  !
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)
  double precision, allocatable :: reigvec_tmp(:,:), leigvec_tmp(:,:)

  integer                       :: i, j, n_degen,k , iteration
  double precision              :: shift_current
  double precision              :: r,thr,accu_d, accu_nd
  integer,          allocatable :: iorder_origin(:),iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:),S(:,:)
  double precision, allocatable :: Aw(:,:),diag_elem(:),A_save(:,:)
  double precision, allocatable :: im_part(:),re_part(:)
  double precision :: accu,thr_cut, thr_norm=1.d0
  double precision, allocatable :: S_nh_inv_half(:,:)
  logical :: complex_root


  thr_cut = 1.d-15
  print*,'Computing the left/right eigenvectors ...'
  print*,'Using the degeneracy splitting algorithm'
 ! initialization
  shift_current = 1.d-15
  iteration = 0 
  print*,'***** iteration = ',iteration


  ! pre-processing the matrix :: sorting by diagonal elements
  allocate(reigvec_tmp(n,n), leigvec_tmp(n,n))
  allocate(diag_elem(n),iorder_origin(n),A_save(n,n))
!  print*,'Aw'
  do i = 1, n
   iorder_origin(i) = i
   diag_elem(i) = A(i,i)
!   write(*,'(100(F16.10,X))')A(:,i)
  enddo
  call dsort(diag_elem, iorder_origin, n)
  do i = 1, n
   do j = 1, n
    A_save(j,i) = A(iorder_origin(j),iorder_origin(i))
   enddo
  enddo

  allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
  allocate(im_part(n),iorder(n))
  allocate( S(n,n) )
  allocate(S_nh_inv_half(n,n))


  Aw = A_save
  call cancel_small_elmts(aw,n,thr_cut)
  call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
  do i = 1, n
   im_part(i) = -dabs(WI(i))
   iorder(i) = i
  enddo
  call dsort(im_part, iorder, n)
  n_real_eigv = 0
  do i = 1, n
    if(dabs(WI(i)).lt.1.d-20)then
      n_real_eigv += 1
    else
!      print*,'Found an imaginary component to eigenvalue'
!      print*,'Re(i) + Im(i)',WR(i),WI(i)
    endif
  enddo
  if(n_real_eigv.ne.n)then
   shift_current = max(10.d0 * dabs(im_part(1)),shift_current*10.d0)
   print*,'Largest imaginary part found in eigenvalues = ',im_part(1)
   print*,'Splitting the degeneracies by ',shift_current
  else
   print*,'All eigenvalues are real !'
  endif


  do while(n_real_eigv.ne.n)
   iteration += 1
   print*,'***** iteration = ',iteration
   if(shift_current.gt.1.d-3)then
    print*,'shift_current > 1.d-3 !!'
    print*,'Your matrix intrinsically contains complex eigenvalues'
    stop
   endif
   Aw = A_save
!   thr_cut = shift_current
   call cancel_small_elmts(Aw,n,thr_cut)
   call split_matrix_degen(Aw,n,shift_current)
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
   n_real_eigv = 0
   do i = 1, n
     if(dabs(WI(i)).lt.1.d-20)then
       n_real_eigv+= 1
     else
!       print*,'Found an imaginary component to eigenvalue'
!       print*,'Re(i) + Im(i)',WR(i),WI(i)
     endif
   enddo
   if(n_real_eigv.ne.n)then
    do i = 1, n
     im_part(i) = -dabs(WI(i))
     iorder(i) = i
    enddo
    call dsort(im_part, iorder, n)
    shift_current = max(10.d0 * dabs(im_part(1)),shift_current*10.d0)
    print*,'Largest imaginary part found in eigenvalues = ',im_part(1)
    print*,'Splitting the degeneracies by ',shift_current
   else
    print*,'All eigenvalues are real !'
   endif
  enddo
  !!!!!!!!!!!!!!!! SORTING THE EIGENVALUES 
  do i = 1, n
   eigval(i) = WR(i)
   iorder(i) = i
  enddo
  call dsort(eigval,iorder,n)
  do i = 1, n
!   print*,'eigval(i) = ',eigval(i)
   reigvec_tmp(:,i) = VR(:,iorder(i))
   leigvec_tmp(:,i) = Vl(:,iorder(i))
  enddo

!!! ONCE ALL EIGENVALUES ARE REAL ::: CHECK BI-ORTHONORMALITY
  !                               check bi-orthogonality
  call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
  print *, ' accu_nd bi-orthog = ', accu_nd
  if(accu_nd .lt. thresh_biorthog_nondiag) then
    print *, ' bi-orthogonality: ok'
  else
    print *, ' '
    print *, ' bi-orthogonality: not imposed yet'
    if(complex_root) then 
     print *, ' '
     print *, ' '
     print *, ' orthog between degen eigenvect' 
     print *, ' '
     ! bi-orthonormalization using orthogonalization of left, right and then QR between left and right
     call impose_orthog_degen_eigvec(n, eigval, reigvec_tmp) ! orthogonalization of reigvec
     call impose_orthog_degen_eigvec(n, eigval, leigvec_tmp) ! orthogonalization of leigvec
     call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S,  thresh_biorthog_diag, thresh_biorthog_nondiag, .false.) 

     if(accu_nd .lt. thresh_biorthog_nondiag) then
       print *, ' bi-orthogonality: ok'
     else 
      print*,'New vectors not bi-orthonormals at ', accu_nd
      call get_inv_half_nonsymmat_diago(S, n, S_nh_inv_half, complex_root)
      if(complex_root)then 
       call impose_biorthog_qr(n, n, leigvec_tmp, reigvec_tmp) ! bi-orthonormalization using QR
      else
       print*,'S^{-1/2} exists !!'
       call bi_ortho_s_inv_half(n,leigvec_tmp,reigvec_tmp,S_nh_inv_half) ! use of S^{-1/2} bi-orthonormalization 
      endif
     endif
    else ! the matrix S^{-1/2} exists
     print*,'S^{-1/2} exists !!'
     call bi_ortho_s_inv_half(n,leigvec_tmp,reigvec_tmp,S_nh_inv_half) ! use of S^{-1/2} bi-orthonormalization 
    endif
    call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
    if(accu_nd .lt. thresh_biorthog_nondiag) then
      print *, ' bi-orthogonality: ok'
    else 
     print*,'New vectors not bi-orthonormals at ',accu_nd
     print*,'Must be a deep problem ...'
     stop
    endif
  endif
 
  !! EIGENVECTORS SORTED AND BI-ORTHONORMAL
  do i = 1, n
   do j = 1, n
    VR(iorder_origin(j),i) = reigvec_tmp(j,i)
    VL(iorder_origin(j),i) = leigvec_tmp(j,i)
   enddo
  enddo

  !! RECOMPUTING THE EIGENVALUES 
  eigval = 0.d0
  do i = 1, n
   iorder(i) = i
   accu = 0.d0
   do j = 1, n
    accu += VL(j,i) * VR(j,i) 
    do k = 1, n
     eigval(i) +=  VL(j,i) * A(j,k) * VR(k,i) 
    enddo
   enddo
   eigval(i) *= 1.d0/accu
!   print*,'eigval(i) = ',eigval(i)
  enddo
  !! RESORT JUST TO BE SURE
  call dsort(eigval, iorder, n)
  do i = 1, n
   do j = 1, n
    reigvec(j,i) = VR(j,iorder(i))
    leigvec(j,i) = VL(j,iorder(i))
   enddo
  enddo
  print*,'Checking for final reigvec/leigvec'
  shift_current = max(1.d-10,shift_current)
  print*,'Thr for eigenvectors = ',shift_current
  call check_EIGVEC(n, n, A, eigval, leigvec, reigvec, shift_current, thr_norm, .false.)
  call check_biorthog(n, n, leigvec, reigvec, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
  print *, ' accu_nd bi-orthog = ', accu_nd
  
  if(accu_nd .lt. thresh_biorthog_nondiag) then
    print *, ' bi-orthogonality: ok'
  else 
   print*,'Something went wrong in non_hrmt_diag_split_degen_bi_orthog'
   print*,'Eigenvectors are not bi orthonormal ..'
   print*,'accu_nd = ',accu_nd
   stop
  endif

end 


subroutine non_hrmt_fock_mat(n, A, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  !
  ! routine returning the eigenvalues and left/right eigenvectors of the TC fock matrix 
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)
  double precision, allocatable :: reigvec_tmp(:,:), leigvec_tmp(:,:)

  integer                       :: i, j, n_degen,k , iteration
  double precision              :: shift_current
  double precision              :: r,thr,accu_d, accu_nd
  integer,          allocatable :: iorder_origin(:),iorder(:)
  double precision, allocatable :: WR(:), WI(:), Vl(:,:), VR(:,:),S(:,:)
  double precision, allocatable :: Aw(:,:),diag_elem(:),A_save(:,:)
  double precision, allocatable :: im_part(:),re_part(:)
  double precision :: accu,thr_cut
  double precision, allocatable :: S_nh_inv_half(:,:)
  logical :: complex_root
  double precision :: thr_norm=1d0


  thr_cut = 1.d-15
  print*,'Computing the left/right eigenvectors ...'
  print*,'Using the degeneracy splitting algorithm'
 ! initialization
  shift_current = 1.d-15
  iteration = 0 
  print*,'***** iteration = ',iteration


  ! pre-processing the matrix :: sorting by diagonal elements
  allocate(reigvec_tmp(n,n), leigvec_tmp(n,n))
  allocate(diag_elem(n),iorder_origin(n),A_save(n,n))
!  print*,'Aw'
  do i = 1, n
   iorder_origin(i) = i
   diag_elem(i) = A(i,i)
!   write(*,'(100(F16.10,X))')A(:,i)
  enddo
  call dsort(diag_elem, iorder_origin, n)
  do i = 1, n
   do j = 1, n
    A_save(j,i) = A(iorder_origin(j),iorder_origin(i))
   enddo
  enddo

  allocate(WR(n), WI(n), VL(n,n), VR(n,n), Aw(n,n))
  allocate(im_part(n),iorder(n))
  allocate( S(n,n) )
  allocate(S_nh_inv_half(n,n))


  Aw = A_save
  call cancel_small_elmts(aw,n,thr_cut)
  call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
  do i = 1, n
   im_part(i) = -dabs(WI(i))
   iorder(i) = i
  enddo
  call dsort(im_part, iorder, n)
  n_real_eigv = 0
  do i = 1, n
    if(dabs(WI(i)).lt.1.d-20)then
      n_real_eigv += 1
    else
!      print*,'Found an imaginary component to eigenvalue'
!      print*,'Re(i) + Im(i)',WR(i),WI(i)
    endif
  enddo
  if(n_real_eigv.ne.n)then
   shift_current = max(10.d0 * dabs(im_part(1)),shift_current*10.d0)
   print*,'Largest imaginary part found in eigenvalues = ',im_part(1)
   print*,'Splitting the degeneracies by ',shift_current
  else
   print*,'All eigenvalues are real !'
  endif


  do while(n_real_eigv.ne.n)
   iteration += 1
   print*,'***** iteration = ',iteration
   if(shift_current.gt.1.d-3)then
    print*,'shift_current > 1.d-3 !!'
    print*,'Your matrix intrinsically contains complex eigenvalues'
    stop
   endif
   Aw = A_save
!   thr_cut = shift_current
   call cancel_small_elmts(Aw,n,thr_cut)
   call split_matrix_degen(Aw,n,shift_current)
   call lapack_diag_non_sym(n,Aw,WR,WI,VL,VR)
   n_real_eigv = 0
   do i = 1, n
     if(dabs(WI(i)).lt.1.d-20)then
       n_real_eigv+= 1
     else
!       print*,'Found an imaginary component to eigenvalue'
!       print*,'Re(i) + Im(i)',WR(i),WI(i)
     endif
   enddo
   if(n_real_eigv.ne.n)then
    do i = 1, n
     im_part(i) = -dabs(WI(i))
     iorder(i) = i
    enddo
    call dsort(im_part, iorder, n)
    shift_current = max(10.d0 * dabs(im_part(1)),shift_current*10.d0)
    print*,'Largest imaginary part found in eigenvalues = ',im_part(1)
    print*,'Splitting the degeneracies by ',shift_current
   else
    print*,'All eigenvalues are real !'
   endif
  enddo
  !!!!!!!!!!!!!!!! SORTING THE EIGENVALUES 
  do i = 1, n
   eigval(i) = WR(i)
   iorder(i) = i
  enddo
  call dsort(eigval,iorder,n)
  do i = 1, n
!   print*,'eigval(i) = ',eigval(i)
   reigvec_tmp(:,i) = VR(:,iorder(i))
   leigvec_tmp(:,i) = Vl(:,iorder(i))
  enddo

!!! ONCE ALL EIGENVALUES ARE REAL ::: CHECK BI-ORTHONORMALITY
  !                               check bi-orthogonality
  call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
  print *, ' accu_nd bi-orthog = ', accu_nd
  if(accu_nd .lt. thresh_biorthog_nondiag) then
    print *, ' bi-orthogonality: ok'
  else
    print *, ' '
    print *, ' bi-orthogonality: not imposed yet'
    print *, ' '
    print *, ' '
    print *, ' Using impose_unique_biorthog_degen_eigvec' 
    print *, ' '
    ! bi-orthonormalization using orthogonalization of left, right and then QR between left and right
    call impose_unique_biorthog_degen_eigvec(n, eigval, mo_coef, leigvec_tmp, reigvec_tmp)
    call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
    print*,'accu_nd = ',accu_nd
    if(accu_nd .lt. thresh_biorthog_nondiag) then
      print *, ' bi-orthogonality: ok'
    else 
     print*,'New vectors not bi-orthonormals at ',accu_nd
     call get_inv_half_nonsymmat_diago(S, n, S_nh_inv_half,complex_root)
     if(complex_root)then 
      print*,'S^{-1/2} does not exits, using QR bi-orthogonalization'
      call impose_biorthog_qr(n, n, leigvec_tmp, reigvec_tmp, S) ! bi-orthonormalization using QR
     else
      print*,'S^{-1/2} exists !!'
      call bi_ortho_s_inv_half(n,leigvec_tmp,reigvec_tmp,S_nh_inv_half) ! use of S^{-1/2} bi-orthonormalization 
     endif
    endif
    call check_biorthog(n, n, leigvec_tmp, reigvec_tmp, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
    if(accu_nd .lt. thresh_biorthog_nondiag) then
      print *, ' bi-orthogonality: ok'
    else 
     print*,'New vectors not bi-orthonormals at ',accu_nd
     print*,'Must be a deep problem ...'
     stop
    endif
  endif
 
  !! EIGENVECTORS SORTED AND BI-ORTHONORMAL
  do i = 1, n
   do j = 1, n
    VR(iorder_origin(j),i) = reigvec_tmp(j,i)
    VL(iorder_origin(j),i) = leigvec_tmp(j,i)
   enddo
  enddo

  !! RECOMPUTING THE EIGENVALUES 
  eigval = 0.d0
  do i = 1, n
   iorder(i) = i
   accu = 0.d0
   do j = 1, n
    accu += VL(j,i) * VR(j,i) 
    do k = 1, n
     eigval(i) +=  VL(j,i) * A(j,k) * VR(k,i) 
    enddo
   enddo
   eigval(i) *= 1.d0/accu
!   print*,'eigval(i) = ',eigval(i)
  enddo
  !! RESORT JUST TO BE SURE
  call dsort(eigval, iorder, n)
  do i = 1, n
   do j = 1, n
    reigvec(j,i) = VR(j,iorder(i))
    leigvec(j,i) = VL(j,iorder(i))
   enddo
  enddo
  print*,'Checking for final reigvec/leigvec'
  shift_current = max(1.d-10,shift_current)
  print*,'Thr for eigenvectors = ',shift_current
  call check_EIGVEC(n, n, A, eigval, leigvec, reigvec, shift_current, thr_norm, .false.)
  call check_biorthog(n, n, leigvec, reigvec, accu_d, accu_nd, S, thresh_biorthog_diag, thresh_biorthog_nondiag, .false.)
  print *, ' accu_nd bi-orthog = ', accu_nd
  
  if(accu_nd .lt. thresh_biorthog_nondiag) then
    print *, ' bi-orthogonality: ok'
  else 
   print*,'Something went wrong in non_hrmt_diag_split_degen_bi_orthog'
   print*,'Eigenvectors are not bi orthonormal ..'
   print*,'accu_nd = ',accu_nd
   stop
  endif

end 


