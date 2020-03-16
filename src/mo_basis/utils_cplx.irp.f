subroutine mo_as_eigvectors_of_mo_matrix_complex(matrix,n,m,label,sign,output)
  !TODO: test this
  implicit none
  integer,intent(in)             :: n,m, sign
  character*(64), intent(in)     :: label
  complex*16, intent(in)   :: matrix(n,m)
  logical, intent(in)            :: output

  integer :: i,j
  double precision, allocatable  :: eigvalues(:)
  complex*16, allocatable  :: mo_coef_new(:,:), R(:,:), A(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: mo_coef_new, R

  call write_time(6)
  if (m /= mo_num) then
    print *, irp_here, ': Error : m/= mo_num'
    stop 1
  endif
  allocate(A(n,m),R(n,m),mo_coef_new(ao_num,m),eigvalues(m))
  if (sign == -1) then
    do j=1,m
      do i=1,n
        A(i,j) = -matrix(i,j)
      enddo
    enddo
  else
    do j=1,m
      do i=1,n
        A(i,j) = matrix(i,j)
      enddo
    enddo
  endif
  mo_coef_new = mo_coef_complex

  call lapack_diag_complex(eigvalues,R,A,n,m)
  if (output) then
    write (6,'(A)')  'MOs are now **'//trim(label)//'**'
    write (6,'(A)') ''
    write (6,'(A)')  'Eigenvalues'
    write (6,'(A)') '-----------'
    write (6,'(A)')  ''
    write (6,'(A)') '======== ================'
  endif
  if (sign == -1) then
    do i=1,m
      eigvalues(i) = -eigvalues(i)
    enddo
  endif
  if (output) then
    do i=1,m
      write (6,'(I8,1X,F16.10)')  i,eigvalues(i)
    enddo
    write (6,'(A)') '======== ================'
    write (6,'(A)')  ''
    write (6,'(A)')  'Fock Matrix'
    write (6,'(A)') '-----------'
    do i=1,n
      write(*,'(200(E24.15))') A(i,:)
    enddo
  endif

  call zgemm('N','N',ao_num,m,m,(1.d0,0.d0),mo_coef_new,size(mo_coef_new,1),R,size(R,1),(0.d0,0.d0),mo_coef_complex,size(mo_coef_complex,1))
  deallocate(A,mo_coef_new,R,eigvalues)
  call write_time(6)

  mo_label = label
end

subroutine mo_as_svd_vectors_of_mo_matrix_complex(matrix,lda,m,n,label)
  !TODO: test this
  implicit none
  integer,intent(in)             :: lda,m,n
  character*(64), intent(in)     :: label
  complex*16, intent(in)   :: matrix(lda,n)

  integer :: i,j
  double precision :: accu
  double precision, allocatable  :: D(:)
  complex*16, allocatable  :: mo_coef_new(:,:), U(:,:), A(:,:), Vt(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: mo_coef_new, U, Vt, A

  call write_time(6)
  if (m /= mo_num) then
    print *, irp_here, ': Error : m/= mo_num'
    stop 1
  endif

  allocate(A(lda,n),U(lda,n),mo_coef_new(ao_num,m),D(m),Vt(lda,n))

  do j=1,n
    do i=1,m
      A(i,j) = matrix(i,j)
    enddo
  enddo
  mo_coef_new = mo_coef_complex

  call svd_complex(A,lda,U,lda,D,Vt,lda,m,n)

  write (6,'(A)') 'MOs are now **'//trim(label)//'**'
  write (6,'(A)')  ''
  write (6,'(A)') 'Eigenvalues'
  write (6,'(A)')  '-----------'
  write (6,'(A)') ''
  write (6,'(A)')  '======== ================ ================'
  write (6,'(A)')  '   MO       Eigenvalue       Cumulative   '
  write (6,'(A)')  '======== ================ ================'

  accu = 0.d0
  do i=1,m
    accu = accu + D(i)
    write (6,'(I8,1X,F16.10,1X,F16.10)')  i,D(i), accu
  enddo
  write (6,'(A)')  '======== ================ ================'
  write (6,'(A)')  ''

  call zgemm('N','N',ao_num,m,m,(1.d0,0.d0),mo_coef_new,size(mo_coef_new,1),U,size(U,1),(0.d0,0.d0),mo_coef_complex,size(mo_coef_complex,1))
  deallocate(A,mo_coef_new,U,Vt,D)
  call write_time(6)

  mo_label = label
end


subroutine mo_as_svd_vectors_of_mo_matrix_eig_complex(matrix,lda,m,n,eig,label)
  !TODO: test this
  implicit none
  integer,intent(in)             :: lda,m,n
  character*(64), intent(in)     :: label
  complex*16, intent(in)   :: matrix(lda,n)
  double precision, intent(out)  :: eig(m)

  integer :: i,j
  double precision :: accu
  double precision, allocatable  :: D(:)
  complex*16, allocatable  :: mo_coef_new(:,:), U(:,:), A(:,:), Vt(:,:), work(:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: mo_coef_new, U, Vt, A

  call write_time(6)
  if (m /= mo_num) then
    print *, irp_here, ': Error : m/= mo_num'
    stop 1
  endif

  allocate(A(lda,n),U(lda,n),mo_coef_new(ao_num,m),D(m),Vt(lda,n))

  do j=1,n
    do i=1,m
      A(i,j) = matrix(i,j)
    enddo
  enddo
  mo_coef_new = mo_coef_complex

  call svd_complex(A,lda,U,lda,D,Vt,lda,m,n)

  write (6,'(A)') 'MOs are now **'//trim(label)//'**'
  write (6,'(A)')  ''
  write (6,'(A)') 'Eigenvalues'
  write (6,'(A)')  '-----------'
  write (6,'(A)') ''
  write (6,'(A)')  '======== ================ ================'
  write (6,'(A)')  '   MO       Eigenvalue       Cumulative   '
  write (6,'(A)')  '======== ================ ================'

  accu = 0.d0
  do i=1,m
    accu = accu + D(i)
    write (6,'(I8,1X,F16.10,1X,F16.10)')  i,D(i), accu
  enddo
  write (6,'(A)')  '======== ================ ================'
  write (6,'(A)')  ''

  call zgemm('N','N',ao_num,m,m,(1.d0,0.d0),mo_coef_new,size(mo_coef_new,1),U,size(U,1),(0.d0,0.d0),mo_coef_complex,size(mo_coef_complex,1))

  do i=1,m
    eig(i) = D(i)
  enddo

  deallocate(A,mo_coef_new,U,Vt,D)
  call write_time(6)

  mo_label = label

end


subroutine mo_coef_new_as_svd_vectors_of_mo_matrix_eig_complex(matrix,lda,m,n,mo_coef_before,eig,mo_coef_new)
  implicit none
  BEGIN_DOC
! You enter with matrix in the MO basis defined with the mo_coef_before. 
!
! You SVD the matrix and set the eigenvectors as mo_coef_new ordered by increasing singular values 
  END_DOC
  integer,intent(in)             :: lda,m,n
  complex*16, intent(in)   :: matrix(lda,n),mo_coef_before(ao_num,m)
  double precision, intent(out)  :: eig(m)
  complex*16, intent(out)  :: mo_coef_new(ao_num,m)

  integer :: i,j
  double precision :: accu
  double precision, allocatable  ::  D(:)
  complex*16, allocatable  ::  mo_coef_tmp(:,:), U(:,:), A(:,:), Vt(:,:), work(:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, Vt, A

  call write_time(6)
  if (m /= mo_num) then
    print *, irp_here, ': Error : m/= mo_num'
    stop 1
  endif

  allocate(A(lda,n),U(lda,n),D(m),Vt(lda,n),mo_coef_tmp(ao_num,mo_num))

  do j=1,n
    do i=1,m
      A(i,j) = matrix(i,j)
    enddo
  enddo
  mo_coef_tmp = mo_coef_before

  call svd_complex(A,lda,U,lda,D,Vt,lda,m,n)

  write (6,'(A)')  ''
  write (6,'(A)') 'Eigenvalues'
  write (6,'(A)')  '-----------'
  write (6,'(A)') ''
  write (6,'(A)')  '======== ================ ================'
  write (6,'(A)')  '   MO       Eigenvalue       Cumulative   '
  write (6,'(A)')  '======== ================ ================'

  accu = 0.d0
  do i=1,m
    accu = accu + D(i)
    write (6,'(I8,1X,F16.10,1X,F16.10)')  i,D(i), accu
  enddo
  write (6,'(A)')  '======== ================ ================'
  write (6,'(A)')  ''

  call zgemm('N','N',ao_num,m,m,(1.d0,0.d0),mo_coef_tmp,size(mo_coef_new,1),U,size(U,1),(0.d0,0.d0),mo_coef_new,size(mo_coef_new,1))

  do i=1,m
    eig(i) = D(i)
  enddo

  deallocate(A,U,Vt,D,mo_coef_tmp)
  call write_time(6)

end


