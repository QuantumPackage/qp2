
subroutine diag_mat_per_fock_degen(fock_diag, mat_ref, n, thr_d, thr_nd, thr_deg, leigvec, reigvec, eigval)


  BEGIN_DOC
  !
  ! subroutine that diagonalizes a matrix mat_ref BY BLOCK
  !
  ! the blocks are defined by the elements having the SAME DEGENERACIES in the entries "fock_diag"
  !
  ! examples : all elements having degeneracy 1 in fock_diag (i.e. not being degenerated) will be treated together
  !
  !          : all elements having degeneracy 2 in fock_diag (i.e. two elements are equal) will be treated together
  !
  !          : all elements having degeneracy 3 in fock_diag (i.e. two elements are equal) will be treated together
  !
  ! etc... the advantage is to guarentee no spurious mixing because of numerical problems. 
  !
  END_DOC

  implicit none
  integer,          intent(in)  :: n
  double precision, intent(in)  :: fock_diag(n), mat_ref(n,n), thr_d, thr_nd, thr_deg
  double precision, intent(out) :: leigvec(n,n), reigvec(n,n), eigval(n)

  integer                       :: n_degen_list, n_degen,size_mat, i, j, k, icount, m, index_degen
  integer                       :: ii, jj, i_good, j_good, n_real
  integer                       :: icount_eigval
  logical,          allocatable :: is_ok(:)
  integer,          allocatable :: list_degen(:,:), list_same_degen(:)
  integer,          allocatable :: iorder(:), list_degen_sorted(:)
  double precision, allocatable :: leigvec_unsrtd(:,:), reigvec_unsrtd(:,:), eigval_unsrtd(:)
  double precision, allocatable :: mat_tmp(:,:), eigval_tmp(:), leigvec_tmp(:,:), reigvec_tmp(:,:)

  allocate(leigvec_unsrtd(n,n), reigvec_unsrtd(n,n), eigval_unsrtd(n))
  leigvec_unsrtd = 0.d0
  reigvec_unsrtd = 0.d0
  eigval_unsrtd  = 0.d0

  ! obtain degeneracies 
  allocate(list_degen(n,0:n))
  call give_degen_full_list(fock_diag, n, thr_deg, list_degen, n_degen_list)

  allocate(iorder(n_degen_list), list_degen_sorted(n_degen_list))
  do i = 1, n_degen_list
    n_degen = list_degen(i,0) 
    list_degen_sorted(i) = n_degen
    iorder(i) = i
  enddo

  ! sort by number of degeneracies 
  call isort(list_degen_sorted, iorder, n_degen_list)

  allocate(is_ok(n_degen_list))
  is_ok = .True.
  icount_eigval = 0

  ! loop over degeneracies 
  do i = 1, n_degen_list
    if(.not.is_ok(i)) cycle

    is_ok(i) = .False.
    n_degen  = list_degen_sorted(i)

    print *, ' diagonalizing for n_degen = ', n_degen

    k = 1

   ! group all the entries having the same degeneracies 
!!  do while (list_degen_sorted(i+k)==n_degen)
    do m = i+1, n_degen_list
      if(list_degen_sorted(m)==n_degen) then
        is_ok(i+k) = .False.
        k += 1
      endif
    enddo

    print *, ' number of identical degeneracies = ', k
    size_mat = k*n_degen  
    print *, ' size_mat = ', size_mat
    allocate(mat_tmp(size_mat,size_mat), list_same_degen(size_mat))
    allocate(eigval_tmp(size_mat), leigvec_tmp(size_mat,size_mat), reigvec_tmp(size_mat,size_mat))
    ! group all the elements sharing the same degeneracy
    icount = 0
    do j = 1, k ! jth set of degeneracy
      index_degen = iorder(i+j-1)
      do m = 1, n_degen
        icount += 1
        list_same_degen(icount) = list_degen(index_degen,m)
      enddo
    enddo

    print *, ' list of elements '
    do icount = 1, size_mat
      print *, icount, list_same_degen(icount)
    enddo

    ! you copy subset of matrix elements having all the same degeneracy in mat_tmp
    do ii = 1, size_mat
      i_good = list_same_degen(ii)
      do jj = 1, size_mat
        j_good = list_same_degen(jj)
        mat_tmp(jj,ii) = mat_ref(j_good,i_good)
      enddo
    enddo

    call non_hrmt_bieig( size_mat, mat_tmp, thr_d, thr_nd &
                       , leigvec_tmp, reigvec_tmp         & 
                       , n_real, eigval_tmp )

    do ii = 1, size_mat
      icount_eigval += 1
      eigval_unsrtd(icount_eigval) = eigval_tmp(ii) ! copy eigenvalues 
      do jj = 1, size_mat ! copy the eigenvectors 
        j_good = list_same_degen(jj)
        leigvec_unsrtd(j_good,icount_eigval) = leigvec_tmp(jj,ii)
        reigvec_unsrtd(j_good,icount_eigval) = reigvec_tmp(jj,ii)
      enddo
    enddo

    deallocate(mat_tmp, list_same_degen)
    deallocate(eigval_tmp, leigvec_tmp, reigvec_tmp)
  enddo

  if(icount_eigval .ne. n) then
    print *, ' pb !! (icount_eigval.ne.n)'
    print *, ' icount_eigval,n', icount_eigval, n
    stop
  endif
 
  deallocate(iorder)
  allocate(iorder(n))
  do i = 1, n
    iorder(i) = i
  enddo
  call dsort(eigval_unsrtd, iorder, n)

  do i = 1, n
    print*,'sorted eigenvalues '
    i_good = iorder(i)
    eigval(i) = eigval_unsrtd(i)
    print*,'i,eigval(i) = ',i,eigval(i)
    do j = 1, n
      leigvec(j,i) = leigvec_unsrtd(j,i_good)
      reigvec(j,i) = reigvec_unsrtd(j,i_good)
    enddo
  enddo

  deallocate(leigvec_unsrtd, reigvec_unsrtd, eigval_unsrtd)
  deallocate(list_degen)
  deallocate(iorder, list_degen_sorted)
  deallocate(is_ok)

end

! ---

subroutine give_degen_full_list(A, n, thr, list_degen, n_degen_list)

  BEGIN_DOC
  ! you enter with an array A(n) and spits out all the elements degenerated up to thr
  !
  ! the elements of A(n) DON'T HAVE TO BE SORTED IN THE ENTRANCE: TOTALLY GENERAL 
  !
  ! list_degen(i,0) = number of degenerate entries 
  !
  ! list_degen(i,1) = index of the first degenerate entry
  !
  ! list_degen(i,2:list_degen(i,0)) = list of all other dengenerate entries 
  !
  ! if list_degen(i,0) == 1 it means that there is no degeneracy for that element
  END_DOC

  implicit none

  double precision, intent(in)  :: A(n)
  double precision, intent(in)  :: thr
  integer,          intent(in)  :: n
  integer,          intent(out) :: list_degen(n,0:n), n_degen_list
  integer                       :: i, j, icount, icheck
  logical, allocatable          :: is_ok(:)


  allocate(is_ok(n))
  n_degen_list = 0
  is_ok = .True.
  do i = 1, n
    if(.not.is_ok(i)) cycle
    n_degen_list +=1
    is_ok(i) = .False.
    list_degen(n_degen_list,1) = i
    icount = 1
    do j = i+1, n
      if(dabs(A(i)-A(j)).lt.thr.and.is_ok(j)) then
        is_ok(j) = .False.
        icount += 1
        list_degen(n_degen_list,icount) = j
      endif
    enddo

    list_degen(n_degen_list,0) = icount
  enddo

  icheck = 0
  do i = 1, n_degen_list
    icheck += list_degen(i,0)
  enddo

  if(icheck.ne.n)then
    print *, ' pb ! :: icheck.ne.n'
    print *, icheck, n
    stop
  endif

end

! ---

