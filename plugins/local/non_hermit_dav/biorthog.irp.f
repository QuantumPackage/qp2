subroutine non_hrmt_bieig_degen(n, A, s_right_in, thr_d, thr_nd, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the EIGENVALUES and corresponding LEFT/RIGHT eigenvetors 
  ! of a non hermitian matrix A(n,n). 
  ! 
  ! IT ALSO TAKES AS INPUT THE OVERLAP MATRIX BETWEEN THE RIGHT BASIS: s_right_in(i,j) = <phi_i|phi_j> where phi_i is a right basis vector
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  ! It also orthogonalizes right degenerate eigenvectors 
  ! 
  ! the left vectors are taken as the inverse of the right eigenvectors in order to have bi orthonormality
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n), s_right_in(n,n)
  double precision, intent(in)  :: thr_d, thr_nd
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j,k, l, n_im_even, ii, m, n_dim
  integer                       :: n_real, n_im
  double precision              :: thr, thr_cut, thr_diag, thr_norm
  double precision              :: accu_d, accu_nd

  integer,          allocatable :: list_real(:), iorder_real(:), deg_num(:), iorder(:)
  integer,          allocatable :: list_im(:), iorder_im(:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:), eigv_real(:), eigv_im(:), eigv_im_im(:)
  double precision, allocatable :: S(:,:), eigv_srtd_real(:), eigv_srtd_im(:)
  double precision, allocatable :: reigvec_srtd_real(:,:), reigvec_srtd_im(:,:)
  double precision, allocatable :: leigvec_srtd_real(:,:), leigvec_srtd_im(:,:)
  double precision, allocatable  :: phi_1_tilde(:),phi_2_tilde(:),chi_1_tilde(:),chi_2_tilde(:)
  logical, allocatable :: is_real(:)
  double precision, allocatable :: vec_tmp(:)
  double precision :: accu 
  double precision, allocatable :: eigval_tmp(:)
  double precision, allocatable :: seigval(:), seigvec(:,:),sphichi(:,:)

  allocate(phi_1_tilde(n),phi_2_tilde(n),chi_1_tilde(n),chi_2_tilde(n))

  allocate(WR(n), WI(n), VL(n,n), VR(n,n), is_real(n),iorder(n) ) 

  call lapack_diag_non_sym(n, A, WR, WI, VL, VR)
!  

  thr_diag = 1d-06
  thr_norm = 1d+10

  ! ---

  ! track & sort the real eigenvalues 

  allocate(list_real(n), list_im(n), eigv_real(n), eigv_im(n),eigv_im_im(n))
  thr    = 0.d0
  n_real = 0
  n_im   = 0
  do i = 1, n
    if(dabs(WI(i)) .le. thr) then
      is_real(i) = .True.
      n_real +=1 
      list_real(n_real) = i
      eigv_real(n_real) = WR(i)
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
      is_real(i) = .False.
      n_im +=1 
      eigv_im(n_im) = WR(i)
      eigv_im_im(n_im) = WI(i)
      list_im(n_im) = i
    endif
  enddo
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REAL PART ONLY

 ! sort the real eigenvalues and copy the eigenvalues and eigenvectors 
 ! in eigv_srtd_real and (leigvec_srtd_real , reigvec_srtd_real)
  allocate(leigvec_srtd_real(n,n_real),reigvec_srtd_real(n,n_real),eigv_srtd_real(n_real)) 
  do i = 1, n_real
   iorder(i) = list_real(i)
  enddo
  call dsort(eigv_real, iorder, n_real)
  do i = 1, n_real
   eigv_srtd_real(i) = eigv_real(i)
   do j = 1, n
    reigvec_srtd_real(j,i) = VR(j,iorder(i))
    leigvec_srtd_real(j,i) = VL(j,iorder(i))
   enddo
  enddo

!!!!! COMPLEX-VALUED EIGENVALUES 
  if(n_im>0)then
   print*,'Number of imaginary eigenvalues found = ',n_im
    n_im_even = n_im/2
    allocate(eigval_tmp(n_im_even))
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMAGINARY COMPONENTS 
   ! sort the complex eigenvalues by their real part and copy the eigenvalues and eigenvectors 
   ! in eigv_srtd_im and (leigvec_srtd_im , reigvec_srtd_im)
    if (iand(n_im,1)==1)then
     call qp_bug(irp_here,n_im,'n_im is odd ')
    endif
    allocate(leigvec_srtd_im(n,n_im),reigvec_srtd_im(n,n_im),eigv_srtd_im(n_im)) 
    do i = 1, n_im, 2
     iorder(i/2+1) = list_im(i)
     eigval_tmp(i/2+1) = eigv_im(i)
    enddo
    call dsort(eigval_tmp, iorder, n_im_even)
    ! first couple 
    do i = 1, n_im
     eigv_srtd_im(i) = eigval_tmp((i-1)/2+1)
    enddo
    do ii = 1,n_im_even
     reigvec_srtd_im(1:n,2*ii-1) =  VR(1:n,iorder(ii))
     reigvec_srtd_im(1:n,2*ii)   =  VR(1:n,iorder(ii)+1)
   enddo
   deallocate(eigval_tmp)
  endif
  
  !! Sorting the union of the real and imaginary eigenvalues 
  double precision, allocatable :: real_eigv_total(:), reigvec_total(:,:), leigvec_total(:,:)
  allocate(real_eigv_total(n), reigvec_total(n,n), leigvec_total(n,n))
  do i = 1, n_real
   real_eigv_total(i) = eigv_real(i)
   reigvec_total(1:n,i) = reigvec_srtd_real(1:n,i)
  enddo
  do ii = 1, n_im
   i=n_real + ii
   real_eigv_total(i) = eigv_srtd_im(ii) 
   reigvec_total(1:n,i) = reigvec_srtd_im(1:n,ii)
  enddo
  do i = 1, n
   iorder(i) = i
  enddo
  call dsort(real_eigv_total, iorder, n)
  do i = 1, n
   eigval(i) = real_eigv_total(i)
   reigvec(1:n,i) = reigvec_total(1:n,iorder(i))
  enddo
 deallocate(iorder)


  ! orthogonalize the degenerate right eigenvectors 
  thr = thr_degen_tc
  allocate(deg_num(n),vec_tmp(n))
  call get_degen(eigval, thr, n, deg_num)
  double precision, allocatable :: reigvec_tmp(:,:)
  integer :: count_mo
  count_mo = 0
  do i = 1, n
    n_dim = deg_num(i)
    count_mo += n_dim
    if(n_dim .le.1.)cycle
    ! buiding the overlap matrix between the degenerate right eigenvectors 
    allocate(sphichi(n_dim, n_dim),seigval(n_dim), seigvec(n_dim, n_dim), reigvec_tmp(n,n_dim), iorder(n_dim))
    sphichi = 0.D0
    k = 0
    do j = i, i+n_dim-1
      k+=1
      iorder(k) = j
    enddo
    do j = 1, n_dim
     do k = 1, n_dim
      do l = 1, n
       do m = 1, n
        sphichi(k,j) += reigvec(l,iorder(k)) * s_right_in(m,l) * reigvec(m,iorder(j))
       enddo
      enddo
     enddo
    enddo
    ! orthogonalize the right eigenvectors by diagonalizing the S matrix
!   call lapack_diagd(seigval,seigvec,sphichi,n_dim,n_dim)
    call pivoted_cholesky( sphichi, n_dim, -1.d0, n_dim, seigvec)
    ! get the orthogonal eigenvectors in the big basis 
    reigvec_tmp = 0.d0
    do j = 1, n_dim
     do l = 1, n_dim 
      do k = 1, n
       reigvec_tmp(k,j) += seigvec(l,j) * reigvec(k,iorder(l))
      enddo
     enddo
    enddo
    ! copy them in the reigvec array 
    do j = 1, n_dim
     reigvec(1:n,iorder(j)) = reigvec_tmp(1:n,j)
    enddo
    deallocate(sphichi,seigval, seigvec, reigvec_tmp, iorder)
    
  enddo
  if(count_mo.ne.n)then
    call qp_bug(irp_here,count_mo,'count_mo is not equal to n !')
  endif
 
  double precision :: accu_tot
  accu_tot = 0.d0
  do i = 1, n ! loop over eigenvectors 
   vec_tmp = 0.d0
   do j = 1, n
    do k = 1, n
     vec_tmp(j) += A(j,k) * reigvec(k,i)
    enddo
   enddo
   do j = 1,n
    vec_tmp(j) -= reigvec(j,i) * eigval(i)
   enddo
   accu = 0.d0
   do j = 1,n
    accu += dabs(vec_tmp(j))
   enddo
   accu_tot += accu
  enddo
  print*,'Average residual of the right eigenvectors ',accu_tot/dble(n)


  double precision, allocatable :: leigvec_tmp(:,:)
  allocate(leigvec_tmp(n,n))
!  call get_pseudo_inverse(reigvec,size(reigvec,1),n,n,leigvec_tmp,size(leigvec_tmp,1),0.d0)
  call get_inverse(reigvec,size(reigvec,1),n,leigvec_tmp,size(leigvec_tmp,1))
  do i = 1, n
   do j = 1, n
    leigvec(i,j) = leigvec_tmp(j,i)
   enddo
  enddo
 

!  ! check bi-orthogonality
!
  thr_diag = 10.d0
  thr_norm = 1d+10
!
  n_real_eigv = n
  allocate( S(n_real_eigv,n_real_eigv) )
  call check_biorthog(n, n, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .True.)
!
  return

end
! ---

subroutine non_hrmt_bieig(n, A, thr_d, thr_nd, leigvec, reigvec, n_real_eigv, eigval)

  BEGIN_DOC
  ! 
  ! routine which returns the sorted REAL EIGENVALUES ONLY and corresponding LEFT/RIGHT eigenvetors 
  ! of a non hermitian matrix A(n,n)
  !
  ! n_real_eigv is the number of real eigenvalues, which might be smaller than the dimension "n" 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(in)  :: thr_d, thr_nd
  integer,          intent(out) :: n_real_eigv
  double precision, intent(out) :: reigvec(n,n), leigvec(n,n), eigval(n)

  integer                       :: i, j,k 
  integer                       :: n_good
  double precision              :: thr, thr_cut, thr_diag, thr_norm
  double precision              :: accu_d, accu_nd

  integer,          allocatable :: list_good(:), iorder(:), deg_num(:)
  double precision, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  double precision, allocatable :: S(:,:)
  double precision, allocatable  :: phi_1_tilde(:),phi_2_tilde(:),chi_1_tilde(:),chi_2_tilde(:)

  allocate(phi_1_tilde(n),phi_2_tilde(n),chi_1_tilde(n),chi_2_tilde(n))

  allocate(WR(n), WI(n), VL(n,n), VR(n,n)) 

  call lapack_diag_non_sym(n, A, WR, WI, VL, VR)

  thr_diag = 1d-06
  thr_norm = 1d+10

  ! ---

  ! track & sort the real eigenvalues 

  n_good = 0
  thr    = Im_thresh_tc
  do i = 1, n
    if(dabs(WI(i)) .lt. thr) then
      n_good += 1
    else
      print*, 'Found an imaginary component to eigenvalue on i = ', i
      print*, 'Re(i) + Im(i)', WR(i), WI(i)
    endif
  enddo

  if(n_good.ne.n) then
    print*,'there are some imaginary eigenvalues '
    thr_diag = 1d-03
    n_good = n
  endif

  allocate(list_good(n_good), iorder(n_good))

  n_good = 0
  do i = 1, n
    n_good += 1
    list_good(n_good) = i
    eigval(n_good) = WR(i)
  enddo

  deallocate( WR, WI )

  n_real_eigv = n_good 
  do i = 1, n_good
    iorder(i) = i
  enddo
  call dsort(eigval, iorder, n_good)
      
  reigvec(:,:) = 0.d0 
  leigvec(:,:) = 0.d0 
  do i = 1, n_real_eigv
    do j = 1, n
      reigvec(j,i) = VR(j,list_good(iorder(i)))
      leigvec(j,i) = VL(j,list_good(iorder(i)))
    enddo
  enddo

  deallocate( list_good, iorder )
  deallocate( VL, VR )

  ASSERT(n==n_real_eigv)

  ! ---

  ! check bi-orthogonality

  thr_diag = 10.d0
  thr_norm = 1d+10

  allocate( S(n_real_eigv,n_real_eigv) )
  call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .false.)

  if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .lt. thr_d) ) then

    print *, ' lapack vectors are normalized and bi-orthogonalized'
    deallocate(S)
    return

  ! accu_nd is modified after adding the normalization
  elseif( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv))/dble(n_real_eigv) .gt. thr_d) ) then

    print *, ' lapack vectors are not normalized but bi-orthogonalized'
    call check_biorthog_binormalize(n, n_real_eigv, leigvec, reigvec, thr_d, thr_nd, .true.)

    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .true.)
    call check_EIGVEC(n, n, A, eigval, leigvec, reigvec, thr_diag, thr_norm, .true.)

    deallocate(S)
    return

  else

    print *, ' lapack vectors are not normalized neither bi-orthogonalized'

    allocate(deg_num(n))
    call reorder_degen_eigvec(n, deg_num, eigval, leigvec, reigvec)
    call impose_biorthog_degen_eigvec(n, deg_num, eigval, leigvec, reigvec)
    deallocate(deg_num)

    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .false.)
    if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(n_real_eigv)) .gt. thr_d) ) then
      call check_biorthog_binormalize(n, n_real_eigv, leigvec, reigvec, thr_d, thr_nd, .true.)
    endif
    call check_biorthog(n, n_real_eigv, leigvec, reigvec, accu_d, accu_nd, S, thr_d, thr_nd, .true.)

    deallocate(S)

  endif

  return

end

! ---

subroutine check_bi_ortho(reigvec, leigvec, n, S, accu_nd)

  BEGIN_DOC
  ! retunrs the overlap matrix S = Leigvec^T Reigvec 
  !
  ! and the square root of the sum of the squared off-diagonal elements of S
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: reigvec(n,n), leigvec(n,n)
  double precision, intent(out) :: S(n,n), accu_nd

  integer :: i,j

  ! S = VL x VR
  call dgemm( 'T', 'N', n, n, n, 1.d0 &
            , leigvec, size(leigvec, 1), reigvec, size(reigvec, 1)  &
            , 0.d0, S, size(S, 1) )

  accu_nd = 0.d0
  do i = 1, n
    do j = 1, n
      if(i.ne.j) then
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

end



subroutine get_degen(energies, thr, n, deg_num)
 implicit none
 BEGIN_DOC 
 ! you enter with energies which are sorted and you get out with 
 ! 
 ! deg_num(i) = degeneraci of the corresponding eigenvalue 
 ! 
 END_DOC
 integer, intent(in) ::  n
 double precision, intent(in) :: energies(n), thr
 integer, intent(out):: deg_num(n)
 
 integer :: i,j
 double precision :: de_thr,ei,ej,de
 do i = 1, n
   deg_num(i) = 1
 enddo

 de_thr = thr

 do i = 1, n-1
   ei = energies(i)

   ! already considered in degen vectors
   if(deg_num(i).eq.0) cycle

   do j = i+1, n
     ej = energies(j)
     de = dabs(ei - ej)

     if(de .lt. de_thr) then
       deg_num(i) = deg_num(i) + 1 
       deg_num(j) = 0
     endif

   enddo
 enddo

end
