
! ---

 BEGIN_PROVIDER [double precision, ao_overlap  , (ao_num, ao_num)]

  BEGIN_DOC
  ! Overlap between atomic basis functions:
  !
  ! :math:`\int \chi_i(r) \chi_j(r) dr`
  END_DOC

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta, c
  double precision :: A_center(3), B_center(3)

  ao_overlap   = 0.d0

  if(read_ao_integrals_overlap) then
   call ezfio_get_ao_one_e_ints_ao_integrals_overlap(ao_overlap(1:ao_num, 1:ao_num))
   print *,  'AO overlap integrals read from disk'
  else
   call ao_cart_to_ao_basis(ao_cart_overlap, ao_cart_num, ao_overlap, ao_num)
  endif

  if (write_ao_integrals_overlap) then
     call ezfio_set_ao_one_e_ints_ao_integrals_overlap(ao_overlap(1:ao_num, 1:ao_num))
     print *,  'AO overlap integrals written to disk'
  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ao_overlap_imag, (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Imaginary part of the overlap
 END_DOC
 ao_overlap_imag = 0.d0
END_PROVIDER

! ---

BEGIN_PROVIDER [ complex*16, ao_overlap_complex, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! Overlap for complex AOs
  END_DOC
  integer                        :: i,j
  do j=1,ao_num
    do i=1,ao_num
      ao_overlap_complex(i,j) = dcmplx( ao_overlap(i,j), ao_overlap_imag(i,j) )
    enddo
  enddo
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ao_overlap_abs, (ao_num, ao_num) ]

  BEGIN_DOC
  ! Overlap between absolute values of atomic basis functions:
  !
  ! :math:`\int |\chi_i(r)| |\chi_j(r)| dr`
  END_DOC

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta
  double precision :: A_center(3), B_center(3)
  double precision :: lower_exp_val, dx

  if(is_periodic) then

    do j = 1, ao_num
      do i = 1, ao_num
        ao_overlap_abs(i,j) = cdabs(ao_overlap_complex(i,j))
      enddo
    enddo

  else
   print*,'todo ! numerical integration on DFT grid !'
   stop

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, S_inv,(ao_num,ao_num) ]
 implicit none
 BEGIN_DOC
! Inverse of the overlap matrix
 END_DOC
 call get_pseudo_inverse(ao_overlap,size(ao_overlap,1),ao_num,ao_num,S_inv, &
    size(S_inv,1),lin_dep_cutoff)
END_PROVIDER

BEGIN_PROVIDER [ complex*16, S_inv_complex,(ao_num,ao_num) ]
 implicit none
 BEGIN_DOC
! Inverse of the overlap matrix
 END_DOC
 call get_pseudo_inverse_complex(ao_overlap_complex, size(ao_overlap_complex,1),&
    ao_num,ao_num,S_inv_complex,size(S_inv_complex,1),lin_dep_cutoff)
END_PROVIDER

BEGIN_PROVIDER [ double precision, S_half_inv, (AO_num,AO_num) ]

  BEGIN_DOC
!   :math:`X = S^{-1/2}` obtained by SVD
  END_DOC

  implicit none

  integer                         :: num_linear_dependencies
  integer                         :: LDA, LDC
  double precision, allocatable   :: U(:,:),Vt(:,:), D(:)
  integer                         :: info, i, j, k
  double precision, parameter     :: threshold_overlap_AO_eigenvalues = 1.d-6

  LDA = size(AO_overlap,1)
  LDC = size(S_half_inv,1)

  allocate(         &
    U(LDC,AO_num),  &
    Vt(LDA,AO_num), &
    D(AO_num))

  call svd(            &
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
      S_half_inv(j,i) = 0.d0
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


BEGIN_PROVIDER [ double precision, S_half, (ao_num,ao_num)  ]
 implicit none
 BEGIN_DOC
 ! :math:`S^{1/2}`
 END_DOC

  integer :: i,j,k
  double precision, allocatable  :: U(:,:)
  double precision, allocatable  :: Vt(:,:)
  double precision, allocatable  :: D(:)

  allocate(U(ao_num,ao_num),Vt(ao_num,ao_num),D(ao_num))

  call svd(ao_overlap,size(ao_overlap,1),U,size(U,1),D,Vt,size(Vt,1),ao_num,ao_num)

  do i=1,ao_num
    D(i) = dsqrt(D(i))
    do j=1,ao_num
      S_half(j,i) = 0.d0
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

