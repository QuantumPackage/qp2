
! ---

 BEGIN_PROVIDER [ double precision, ao_overlap  , (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_overlap_x, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_overlap_y, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_overlap_z, (ao_num, ao_num) ]

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
  ao_overlap_x = 0.d0
  ao_overlap_y = 0.d0
  ao_overlap_z = 0.d0

  if(read_ao_integrals_overlap) then

    call ezfio_get_ao_one_e_ints_ao_integrals_overlap(ao_overlap(1:ao_num, 1:ao_num))
    print *,  'AO overlap integrals read from disk'

  else

    if(use_cosgtos) then
      !print*, ' use_cosgtos for ao_overlap ?', use_cosgtos

      do j = 1, ao_num
        do i = 1, ao_num
          ao_overlap  (i,j) = ao_overlap_cosgtos  (i,j) 
          ao_overlap_x(i,j) = ao_overlap_cosgtos_x(i,j)
          ao_overlap_y(i,j) = ao_overlap_cosgtos_y(i,j)
          ao_overlap_z(i,j) = ao_overlap_cosgtos_z(i,j)
        enddo
      enddo

    else

      dim1=100
      !$OMP PARALLEL DO SCHEDULE(GUIDED) &
      !$OMP DEFAULT(NONE) &
      !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
      !$OMP  overlap_x,overlap_y, overlap_z, overlap, &
      !$OMP  alpha, beta,i,j,c) &
      !$OMP SHARED(nucl_coord,ao_power,ao_prim_num, &
      !$OMP  ao_overlap_x,ao_overlap_y,ao_overlap_z,ao_overlap,ao_num,ao_coef_normalized_ordered_transp,ao_nucl, &
      !$OMP  ao_expo_ordered_transp,dim1)
      do j=1,ao_num
      A_center(1) = nucl_coord( ao_nucl(j), 1 )
      A_center(2) = nucl_coord( ao_nucl(j), 2 )
      A_center(3) = nucl_coord( ao_nucl(j), 3 )
      power_A(1)  = ao_power( j, 1 )
      power_A(2)  = ao_power( j, 2 )
      power_A(3)  = ao_power( j, 3 )
      do i= 1,ao_num
        B_center(1) = nucl_coord( ao_nucl(i), 1 )
        B_center(2) = nucl_coord( ao_nucl(i), 2 )
        B_center(3) = nucl_coord( ao_nucl(i), 3 )
        power_B(1)  = ao_power( i, 1 )
        power_B(2)  = ao_power( i, 2 )
        power_B(3)  = ao_power( i, 3 )
        do n = 1,ao_prim_num(j)
        alpha = ao_expo_ordered_transp(n,j)
        do l = 1, ao_prim_num(i)
          beta = ao_expo_ordered_transp(l,i)
          call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
          c = ao_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i)
          ao_overlap(i,j) += c * overlap
          if(isnan(ao_overlap(i,j)))then
           print*,'i,j',i,j
           print*,'l,n',l,n
           print*,'c,overlap',c,overlap
           print*,overlap_x,overlap_y,overlap_z
           stop
          endif
          ao_overlap_x(i,j) += c * overlap_x
          ao_overlap_y(i,j) += c * overlap_y
          ao_overlap_z(i,j) += c * overlap_z
        enddo
        enddo
      enddo
      enddo
      !$OMP END PARALLEL DO

    endif

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

    dim1=100
    lower_exp_val = 40.d0
    !$OMP PARALLEL DO SCHEDULE(GUIDED)                               &
        !$OMP DEFAULT(NONE)                                          &
        !$OMP PRIVATE(A_center,B_center,power_A,power_B,             &
        !$OMP  overlap_x,overlap_y, overlap_z,                       &
        !$OMP  alpha, beta,i,j,dx)                                   &
        !$OMP SHARED(nucl_coord,ao_power,ao_prim_num,                &
        !$OMP  ao_overlap_abs,ao_num,ao_coef_normalized_ordered_transp,ao_nucl,&
        !$OMP  ao_expo_ordered_transp,dim1,lower_exp_val)
    do j=1,ao_num
      A_center(1) = nucl_coord( ao_nucl(j), 1 )
      A_center(2) = nucl_coord( ao_nucl(j), 2 )
      A_center(3) = nucl_coord( ao_nucl(j), 3 )
      power_A(1)  = ao_power( j, 1 )
      power_A(2)  = ao_power( j, 2 )
      power_A(3)  = ao_power( j, 3 )
      do i= 1,ao_num
        ao_overlap_abs(i,j)= 0.d0
        B_center(1) = nucl_coord( ao_nucl(i), 1 )
        B_center(2) = nucl_coord( ao_nucl(i), 2 )
        B_center(3) = nucl_coord( ao_nucl(i), 3 )
        power_B(1)  = ao_power( i, 1 )
        power_B(2)  = ao_power( i, 2 )
        power_B(3)  = ao_power( i, 3 )
        do n = 1,ao_prim_num(j)
          alpha = ao_expo_ordered_transp(n,j)
          do l = 1, ao_prim_num(i)
            beta = ao_expo_ordered_transp(l,i)
            call overlap_x_abs(A_center(1),B_center(1),alpha,beta,power_A(1),power_B(1),overlap_x,lower_exp_val,dx,dim1)
            call overlap_x_abs(A_center(2),B_center(2),alpha,beta,power_A(2),power_B(2),overlap_y,lower_exp_val,dx,dim1)
            call overlap_x_abs(A_center(3),B_center(3),alpha,beta,power_A(3),power_B(3),overlap_z,lower_exp_val,dx,dim1)
            ao_overlap_abs(i,j) += abs(ao_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i)) * overlap_x * overlap_y * overlap_z
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

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

