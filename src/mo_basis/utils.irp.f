subroutine save_mos
  implicit none
  double precision, allocatable  :: buffer(:,:)
  integer                        :: i,j

  call ezfio_set_mo_basis_mo_num(mo_num)
  call ezfio_set_mo_basis_mo_label(mo_label)
  call ezfio_set_mo_basis_ao_md5(ao_md5)
  allocate ( buffer(ao_num,mo_num) )
  buffer = 0.d0
  do j = 1, mo_num
    do i = 1, ao_num
      buffer(i,j) = mo_coef(i,j)
    enddo
  enddo
  call ezfio_set_mo_basis_mo_coef(buffer)
  call ezfio_set_mo_basis_mo_occ(mo_occ)
  call ezfio_set_mo_basis_mo_class(mo_class)
  deallocate (buffer)

end


subroutine save_mos_no_occ
  implicit none
  double precision, allocatable  :: buffer(:,:)
  integer                        :: i,j

!  call system('$QP_ROOT/scripts/save_current_mos.sh '//trim(ezfio_filename))
 !call ezfio_set_mo_basis_mo_num(mo_num)
 !call ezfio_set_mo_basis_mo_label(mo_label)
 !call ezfio_set_mo_basis_ao_md5(ao_md5)
  allocate ( buffer(ao_num,mo_num) )
  buffer = 0.d0
  do j = 1, mo_num
    do i = 1, ao_num
      buffer(i,j) = mo_coef(i,j)
    enddo
  enddo
  call ezfio_set_mo_basis_mo_coef(buffer)
  deallocate (buffer)

end

subroutine save_mos_truncated(n)
  implicit none
  double precision, allocatable  :: buffer(:,:)
  integer                        :: i,j,n

!  call system('$QP_ROOT/scripts/save_current_mos.sh '//trim(ezfio_filename))

  call ezfio_set_mo_basis_mo_num(n)
  call ezfio_set_mo_basis_mo_label(mo_label)
  call ezfio_set_mo_basis_ao_md5(ao_md5)
  allocate ( buffer(ao_num,n) )
  buffer = 0.d0
  do j = 1, n
    do i = 1, ao_num
      buffer(i,j) = mo_coef(i,j)
    enddo
  enddo
  call ezfio_set_mo_basis_mo_coef(buffer)
  call ezfio_set_mo_basis_mo_occ(mo_occ)
  call ezfio_set_mo_basis_mo_class(mo_class)
  deallocate (buffer)

end

subroutine mo_as_eigvectors_of_mo_matrix(matrix,n,m,label,sign,output)
  implicit none
  integer,intent(in)             :: n,m, sign
  character*(64), intent(in)     :: label
  double precision, intent(in)   :: matrix(n,m)
  logical, intent(in)            :: output

  integer :: i,j
  double precision, allocatable  :: mo_coef_new(:,:), R(:,:),eigvalues(:), A(:,:)
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
  mo_coef_new = mo_coef

  call lapack_diag(eigvalues,R,A,n,m)
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
  endif

  call dgemm('N','N',ao_num,m,m,1.d0,mo_coef_new,size(mo_coef_new,1),R,size(R,1),0.d0,mo_coef,size(mo_coef,1))
  deallocate(A,mo_coef_new,R,eigvalues)
  call write_time(6)

  mo_label = label
end

subroutine mo_as_svd_vectors_of_mo_matrix(matrix,lda,m,n,label)
  implicit none
  integer,intent(in)             :: lda,m,n
  character*(64), intent(in)     :: label
  double precision, intent(in)   :: matrix(lda,n)

  integer :: i,j
  double precision :: accu
  double precision, allocatable  :: mo_coef_new(:,:), U(:,:),D(:), A(:,:), Vt(:,:), work(:)
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
  mo_coef_new = mo_coef

  call svd(A,lda,U,lda,D,Vt,lda,m,n)

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

  call dgemm('N','N',ao_num,m,m,1.d0,mo_coef_new,size(mo_coef_new,1),U,size(U,1),0.d0,mo_coef,size(mo_coef,1))
  deallocate(A,mo_coef_new,U,Vt,D)
  call write_time(6)

  mo_label = label
end

subroutine mo_as_svd_vectors_of_mo_matrix_eig(matrix,lda,m,n,eig,label)
  implicit none
  integer,intent(in)             :: lda,m,n
  character*(64), intent(in)     :: label
  double precision, intent(in)   :: matrix(lda,n)
  double precision, intent(out)  :: eig(m)

  integer :: i,j
  double precision :: accu
  double precision, allocatable  :: mo_coef_new(:,:), U(:,:),D(:), A(:,:), Vt(:,:), work(:)
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
  mo_coef_new = mo_coef

  call svd(A,lda,U,lda,D,Vt,lda,m,n)

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

  call dgemm('N','N',ao_num,m,m,1.d0,mo_coef_new,size(mo_coef_new,1),U,size(U,1),0.d0,mo_coef,size(mo_coef,1))

  do i=1,m
    eig(i) = D(i)
  enddo

  deallocate(A,mo_coef_new,U,Vt,D)
  call write_time(6)

  mo_label = label

end


subroutine mo_coef_new_as_svd_vectors_of_mo_matrix_eig(matrix,lda,m,n,mo_coef_before,eig,mo_coef_new)
  implicit none
  BEGIN_DOC
! You enter with matrix in the MO basis defined with the mo_coef_before. 
!
! You SVD the matrix and set the eigenvectors as mo_coef_new ordered by increasing singular values 
  END_DOC
  integer,intent(in)             :: lda,m,n
  double precision, intent(in)   :: matrix(lda,n),mo_coef_before(ao_num,m)
  double precision, intent(out)  :: eig(m),mo_coef_new(ao_num,m)

  integer :: i,j
  double precision :: accu
  double precision, allocatable  ::  mo_coef_tmp(:,:), U(:,:),D(:), A(:,:), Vt(:,:), work(:)
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

  call svd(A,lda,U,lda,D,Vt,lda,m,n)

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

  call dgemm('N','N',ao_num,m,m,1.d0,mo_coef_tmp,size(mo_coef_new,1),U,size(U,1),0.d0,mo_coef_new,size(mo_coef_new,1))

  do i=1,m
    eig(i) = D(i)
  enddo

  deallocate(A,U,Vt,D,mo_coef_tmp)
  call write_time(6)

end


