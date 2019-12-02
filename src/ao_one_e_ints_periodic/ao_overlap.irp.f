BEGIN_PROVIDER [ complex*16, ao_overlap,(ao_num,ao_num) ]
  implicit none
  BEGIN_DOC  
! Overlap between atomic basis functions:
! :math:`\int \chi_i(r) \chi_j(r) dr)`
  END_DOC
  if (read_ao_integrals_overlap) then
    call read_one_e_integrals_complex('ao_overlap', ao_overlap,&
        size(ao_overlap,1), size(ao_overlap,2))
    print *,  'AO overlap integrals read from disk'
  else
    print *, 'complex AO overlap integrals must be provided'
  endif
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_overlap_abs,(ao_num,ao_num) ]
  implicit none
  BEGIN_DOC  
! Overlap between absolute value of atomic basis functions:
! :math:`\int |\chi_i(r)| |\chi_j(r)| dr)`
  END_DOC
  integer :: i,j
  !$OMP PARALLEL DO SCHEDULE(GUIDED) &
  !$OMP DEFAULT(NONE) &
  !$OMP PRIVATE(i,j) &
  !$OMP SHARED(ao_overlap_abs,ao_overlap,ao_num)
  do j=1,ao_num
   do i= 1,ao_num
    ao_overlap_abs(i,j)= cdabs(ao_overlap(i,j))
   enddo
  enddo
  !$OMP END PARALLEL DO
END_PROVIDER

BEGIN_PROVIDER [ complex*16, S_inv,(ao_num,ao_num) ]
 implicit none
 BEGIN_DOC
! Inverse of the overlap matrix
 END_DOC
 call get_pseudo_inverse_complex(ao_overlap,size(ao_overlap,1),ao_num,ao_num,S_inv,size(S_inv,1))
END_PROVIDER

BEGIN_PROVIDER [ complex*16, S_half_inv, (AO_num,AO_num) ]

  BEGIN_DOC
!   :math:`X = S^{-1/2}` obtained by SVD
  END_DOC

  implicit none

  integer                         :: num_linear_dependencies
  integer                         :: LDA, LDC
  double precision, allocatable   :: D(:)
  complex*16, allocatable   :: U(:,:),Vt(:,:)
  integer                         :: info, i, j, k
  double precision, parameter     :: threshold_overlap_AO_eigenvalues = 1.d-6

  LDA = size(AO_overlap,1)
  LDC = size(S_half_inv,1)

  allocate(         &
    U(LDC,AO_num),  &
    Vt(LDA,AO_num), &
    D(AO_num))

  call svd_complex(    &
       AO_overlap,LDA, &
       U,LDC,          &
       D,              &
       Vt,LDA,         &
       AO_num,AO_num)

  num_linear_dependencies = 0
  do i=1,AO_num
    print*,D(i)
    if(abs(D(i)) <= threshold_overlap_AO_eigenvalues) then
      D(i) = 0.d0
      num_linear_dependencies += 1
    else
      ASSERT (D(i) > 0.d0)
      D(i) = 1.d0/sqrt(D(i))
    endif
    do j=1,AO_num
      S_half_inv(j,i) = (0.d0,0.d0)
    enddo
  enddo
  write(*,*) 'linear dependencies',num_linear_dependencies

  do k=1,AO_num
    if(D(k) /= 0.d0) then
      do j=1,AO_num
        do i=1,AO_num
          S_half_inv(i,j) = S_half_inv(i,j) + U(i,k)*D(k)*Vt(k,j)
        enddo
      enddo
    endif
  enddo


END_PROVIDER

BEGIN_PROVIDER [ complex*16, S_half, (ao_num,ao_num)  ]
 implicit none
 BEGIN_DOC
 ! :math:`S^{1/2}`
 END_DOC

  integer :: i,j,k
  complex*16, allocatable  :: U(:,:)
  complex*16, allocatable  :: Vt(:,:)
  double precision, allocatable  :: D(:)

  allocate(U(ao_num,ao_num),Vt(ao_num,ao_num),D(ao_num))

  call svd_complex(ao_overlap,size(ao_overlap,1),U,size(U,1),D,Vt,size(Vt,1),ao_num,ao_num)

  do i=1,ao_num
    D(i) = dsqrt(D(i))
    do j=1,ao_num
      S_half(j,i) = (0.d0,0.d0)
    enddo
  enddo

  do k=1,ao_num
      do j=1,ao_num
        do i=1,ao_num
          S_half(i,j) = S_half(i,j) + U(i,k)*D(k)*Vt(k,j)
        enddo
      enddo
  enddo

  deallocate(U,Vt,D)

END_PROVIDER

