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
    !write (6,'(A)')  'Fock Matrix'
    !write (6,'(A)') '-----------'
    !do i=1,n
    !  write(*,'(200(E24.15))') A(i,:)
    !enddo
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

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!

subroutine mo_as_eigvectors_of_mo_matrix_kpts(matrix,n,m,nk,label,sign,output)
  !TODO: test this
  implicit none
  integer,intent(in)             :: n,m,nk, sign
  character*(64), intent(in)     :: label
  complex*16, intent(in)   :: matrix(n,m,nk)
  logical, intent(in)            :: output

  integer :: i,j,k
  double precision, allocatable  :: eigvalues(:)
  complex*16, allocatable  :: mo_coef_new(:,:), R(:,:), A(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: mo_coef_new, R

  call write_time(6)
  if (m /= mo_num_per_kpt) then
    print *, irp_here, ': Error : m/= mo_num_per_kpt'
    stop 1
  endif
  if (nk /= kpt_num) then
    print *, irp_here, ': Error : nk/= kpt_num'
    stop 1
  endif
  allocate(A(n,m),R(n,m),mo_coef_new(ao_num_per_kpt,m),eigvalues(m))
  do k=1,nk
    if (sign == -1) then
      do j=1,m
        do i=1,n
          A(i,j) = -matrix(i,j,k)
        enddo
      enddo
    else
      do j=1,m
        do i=1,n
          A(i,j) = matrix(i,j,k)
        enddo
      enddo
    endif
    mo_coef_new = mo_coef_kpts(:,:,k)

    call lapack_diag_complex(eigvalues,R,A,n,m)
    if (sign == -1) then
      do i=1,m
        eigvalues(i) = -eigvalues(i)
      enddo
    endif
    if (output) then
      do i=1,m
        write (6,'(2(I8),1X,F16.10)') k,i,eigvalues(i)
      enddo
      write (6,'(A)') '======== ================'
      write (6,'(A)')  ''
      !write (6,'(A)')  'Fock Matrix'
      !write (6,'(A)') '-----------'
      !do i=1,n
      !  write(*,'(200(E24.15))') A(i,:)
      !enddo
    endif

    call zgemm('N','N',ao_num_per_kpt,m,m,(1.d0,0.d0), &
        mo_coef_new,size(mo_coef_new,1),R,size(R,1),(0.d0,0.d0), &
        mo_coef_kpts(:,:,k),size(mo_coef_kpts,1))
  enddo
  deallocate(A,mo_coef_new,R,eigvalues)
  call write_time(6)

  mo_label = label
    if (output) then
      write (6,'(A)')  'MOs are now **'//trim(label)//'**'
      write (6,'(A)') ''
      write (6,'(A)')  'Eigenvalues'
      write (6,'(A)') '-----------'
      write (6,'(A)')  ''
      write (6,'(A)') '======== ================'
    endif
end

subroutine mo_as_eigvectors_of_mo_matrix_kpts_real(matrix,n,m,nk,label,sign,output)
  !TODO: test this
  implicit none
  integer,intent(in)             :: n,m,nk, sign
  character*(64), intent(in)     :: label
  double precision, intent(in)   :: matrix(n,m,nk)
  logical, intent(in)            :: output

  integer :: i,j,k
  double precision, allocatable  :: eigvalues(:)
  !complex*16, allocatable  :: mo_coef_new(:,:)
  double precision, allocatable  :: mo_coef_new(:,:),mo_coef_tmp(:,:),R(:,:), A(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: mo_coef_new, R

  call write_time(6)
  if (m /= mo_num_per_kpt) then
    print *, irp_here, ': Error : m/= mo_num_per_kpt'
    stop 1
  endif
  if (nk /= kpt_num) then
    print *, irp_here, ': Error : nk/= kpt_num'
    stop 1
  endif
  allocate(A(n,m),R(n,m),mo_coef_tmp(ao_num_per_kpt,m),mo_coef_new(ao_num_per_kpt,m),eigvalues(m))
  do k=1,nk
    if (sign == -1) then
      do j=1,m
        do i=1,n
          A(i,j) = -matrix(i,j,k)
        enddo
      enddo
    else
      do j=1,m
        do i=1,n
          A(i,j) = matrix(i,j,k)
        enddo
      enddo
    endif
    mo_coef_new = dble(mo_coef_kpts(:,:,k))

    call lapack_diag(eigvalues,R,A,n,m)
    if (sign == -1) then
      do i=1,m
        eigvalues(i) = -eigvalues(i)
      enddo
    endif
    if (output) then
      do i=1,m
        write (6,'(2(I8),1X,F16.10)') k,i,eigvalues(i)
      enddo
      write (6,'(A)') '======== ================'
      write (6,'(A)')  ''
      !write (6,'(A)')  'Fock Matrix'
      !write (6,'(A)') '-----------'
      !do i=1,n
      !  write(*,'(200(E24.15))') A(i,:)
      !enddo
    endif

    call dgemm('N','N',ao_num_per_kpt,m,m,1.d0, &
        mo_coef_new,size(mo_coef_new,1),R,size(R,1),0.d0, &
        mo_coef_tmp,size(mo_coef_tmp,1))
    call zlacp2('N',ao_num_per_kpt,m,mo_coef_tmp,size(mo_coef_tmp,1), &
                mo_coef_kpts(:,:,k),size(mo_coef_kpts,1))
  enddo
  deallocate(A,mo_coef_new,mo_coef_tmp,R,eigvalues)
  call write_time(6)

  mo_label = label
    if (output) then
      write (6,'(A)')  'MOs are now **'//trim(label)//'**'
      write (6,'(A)') ''
      write (6,'(A)')  'Eigenvalues'
      write (6,'(A)') '-----------'
      write (6,'(A)')  ''
      write (6,'(A)') '======== ================'
    endif
end

subroutine mo_as_svd_vectors_of_mo_matrix_kpts(matrix,lda,m,n,label)
  !TODO: implement
  print *, irp_here, ' not implemented for kpts'
  stop 1
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


subroutine mo_as_svd_vectors_of_mo_matrix_eig_kpts(matrix,lda,m,n,nk,eig,label)
  !TODO: implement
  !print *, irp_here, ' not implemented for kpts'
  !stop 1
  implicit none
  integer,intent(in)             :: lda,m,n,nk
  character*(64), intent(in)     :: label
  complex*16, intent(in)   :: matrix(lda,n,nk)
  double precision, intent(out)  :: eig(m,nk)

  integer :: i,j,k
  double precision :: accu
  double precision, allocatable  :: D(:)
  complex*16, allocatable  :: mo_coef_new(:,:), U(:,:), A(:,:), Vt(:,:), work(:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: mo_coef_new, U, Vt, A

  call write_time(6)
  if (m /= mo_num_per_kpt) then
    print *, irp_here, ': Error : m/= mo_num_per_kpt'
    stop 1
  endif

  
  allocate(A(lda,n),U(lda,n),mo_coef_new(ao_num_per_kpt,m),D(m),Vt(lda,n))

  do k=1,nk
    do j=1,n
      do i=1,m
        A(i,j) = matrix(i,j,k)
      enddo
    enddo
    mo_coef_new(1:ao_num_per_kpt,1:m) = mo_coef_kpts(1:ao_num_per_kpt,1:m,k)

    call svd_complex(A,lda,U,lda,D,Vt,lda,m,n)



    call zgemm('N','N',ao_num_per_kpt,m,m,  &
               (1.d0,0.d0),mo_coef_new,size(mo_coef_new,1),U,size(U,1),&
               (0.d0,0.d0),mo_coef_kpts(1,1,k),size(mo_coef_kpts,1))

    do i=1,m
      eig(i,k) = D(i)
    enddo
    !do j=1,mo_num_per_kpt
    !  do i=1,mo_num_per_kpt
    !    print'(3(I5),2(E25.15))',i,j,k,mo_coef_kpts(i,j,k)
    !  enddo
    !enddo
  enddo

  deallocate(A,mo_coef_new,U,Vt,D)

  write (6,'(A)') 'MOs are now **'//trim(label)//'**'
  write (6,'(A)')  ''
  write (6,'(A)') 'Eigenvalues '
  write (6,'(A)')  '-----------'
  write (6,'(A)') ''
  write (6,'(A)')  '======== ================ ================'
  write (6,'(A)')  '   MO       Eigenvalue       Cumulative   '
  write (6,'(A)')  '======== ================ ================'
 
  do k=1,nk
    accu = 0.d0
    do i=1,m
      accu = accu + eig(i,k)
      write (6,'(I8,1X,F16.10,1X,F16.10)')  i,eig(i,k), accu
    enddo
  write (6,'(A)')  '-------- ---------------- ----------------'
  enddo
  write (6,'(A)')  '======== ================ ================'
  write (6,'(A)')  ''

  call write_time(6)

  mo_label = label

end


subroutine mo_coef_new_as_svd_vectors_of_mo_matrix_eig_kpts(matrix,lda,m,n,mo_coef_before,eig,mo_coef_new)
  !TODO: implement
  print *, irp_here, ' not implemented for kpts'
  stop 1
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

