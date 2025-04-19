
! ---

 BEGIN_PROVIDER [double precision, ao_cart_overlap  , (ao_cart_num, ao_cart_num)]
&BEGIN_PROVIDER [double precision, ao_cart_overlap_x, (ao_cart_num, ao_cart_num)]
&BEGIN_PROVIDER [double precision, ao_cart_overlap_y, (ao_cart_num, ao_cart_num)]
&BEGIN_PROVIDER [double precision, ao_cart_overlap_z, (ao_cart_num, ao_cart_num)]

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

  ao_cart_overlap   = 0.d0
  ao_cart_overlap_x = 0.d0
  ao_cart_overlap_y = 0.d0
  ao_cart_overlap_z = 0.d0

  if(read_ao_cart_integrals_overlap) then

    call ezfio_get_ao_cart_one_e_ints_ao_cart_integrals_overlap(ao_cart_overlap(1:ao_cart_num, 1:ao_cart_num))
    print *,  'AO overlap integrals read from disk'

  else

    if(use_cgtos) then

      do j = 1, ao_cart_num
        do i = 1, ao_cart_num
          ao_cart_overlap  (i,j) = ao_cart_overlap_cgtos  (i,j) 
          ao_cart_overlap_x(i,j) = ao_cart_overlap_cgtos_x(i,j)
          ao_cart_overlap_y(i,j) = ao_cart_overlap_cgtos_y(i,j)
          ao_cart_overlap_z(i,j) = ao_cart_overlap_cgtos_z(i,j)
        enddo
      enddo

    else

      dim1=100
      !$OMP PARALLEL DO SCHEDULE(GUIDED) &
      !$OMP DEFAULT(NONE) &
      !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
      !$OMP  overlap_x,overlap_y, overlap_z, overlap, &
      !$OMP  alpha, beta,i,j,n,l,c) &
      !$OMP SHARED(nucl_coord,ao_cart_power,ao_cart_prim_num, &
      !$OMP  ao_cart_overlap_x,ao_cart_overlap_y,ao_cart_overlap_z,ao_cart_overlap,ao_cart_num,ao_cart_coef_normalized_ordered_transp,ao_cart_nucl, &
      !$OMP  ao_cart_expo_ordered_transp,dim1)
      do j=1,ao_cart_num
      A_center(1) = nucl_coord( ao_cart_nucl(j), 1 )
      A_center(2) = nucl_coord( ao_cart_nucl(j), 2 )
      A_center(3) = nucl_coord( ao_cart_nucl(j), 3 )
      power_A(1)  = ao_cart_power( j, 1 )
      power_A(2)  = ao_cart_power( j, 2 )
      power_A(3)  = ao_cart_power( j, 3 )
      do i= 1,ao_cart_num
        B_center(1) = nucl_coord( ao_cart_nucl(i), 1 )
        B_center(2) = nucl_coord( ao_cart_nucl(i), 2 )
        B_center(3) = nucl_coord( ao_cart_nucl(i), 3 )
        power_B(1)  = ao_cart_power( i, 1 )
        power_B(2)  = ao_cart_power( i, 2 )
        power_B(3)  = ao_cart_power( i, 3 )
        do n = 1,ao_cart_prim_num(j)
        alpha = ao_cart_expo_ordered_transp(n,j)
        do l = 1, ao_cart_prim_num(i)
          beta = ao_cart_expo_ordered_transp(l,i)
          call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
          c = ao_cart_coef_normalized_ordered_transp(n,j) * ao_cart_coef_normalized_ordered_transp(l,i)
          ao_cart_overlap(i,j) += c * overlap
          if(isnan(ao_cart_overlap(i,j)))then
           print*,'i,j',i,j
           print*,'l,n',l,n
           print*,'c,overlap',c,overlap
           print*,overlap_x,overlap_y,overlap_z
           stop
          endif
          ao_cart_overlap_x(i,j) += c * overlap_x
          ao_cart_overlap_y(i,j) += c * overlap_y
          ao_cart_overlap_z(i,j) += c * overlap_z
        enddo
        enddo
      enddo
      enddo
      !$OMP END PARALLEL DO

    endif

  endif

  if (write_ao_cart_integrals_overlap) then
     call ezfio_set_ao_cart_one_e_ints_ao_cart_integrals_overlap(ao_cart_overlap(1:ao_cart_num, 1:ao_cart_num))
     print *,  'AO overlap integrals written to disk'
  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ao_cart_overlap_imag, (ao_cart_num, ao_cart_num) ]
 implicit none
 BEGIN_DOC
 ! Imaginary part of the overlap
 END_DOC
 ao_cart_overlap_imag = 0.d0
END_PROVIDER

! ---

BEGIN_PROVIDER [ complex*16, ao_cart_overlap_complex, (ao_cart_num, ao_cart_num) ]
  implicit none
  BEGIN_DOC
  ! Overlap for complex AOs
  END_DOC
  integer                        :: i,j
  do j=1,ao_cart_num
    do i=1,ao_cart_num
      ao_cart_overlap_complex(i,j) = dcmplx( ao_cart_overlap(i,j), ao_cart_overlap_imag(i,j) )
    enddo
  enddo
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ao_cart_overlap_abs, (ao_cart_num, ao_cart_num) ]

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

    do j = 1, ao_cart_num
      do i = 1, ao_cart_num
        ao_cart_overlap_abs(i,j) = cdabs(ao_cart_overlap_complex(i,j))
      enddo
    enddo

  else

    dim1=100
    lower_exp_val = 40.d0
 !$OMP PARALLEL DO SCHEDULE(GUIDED)               &
 !$OMP DEFAULT(NONE)                              &
 !$OMP PRIVATE(A_center,B_center,power_A,power_B, &
 !$OMP  overlap_x,overlap_y, overlap_z,           &
 !$OMP  alpha, beta,i,j,dx)                       &
 !$OMP SHARED(nucl_coord,ao_cart_power,ao_cart_prim_num,    &
 !$OMP  ao_cart_overlap_abs,ao_cart_num,ao_cart_coef_normalized_ordered_transp,ao_cart_nucl,&
 !$OMP  ao_cart_expo_ordered_transp,dim1,lower_exp_val)
    do j=1,ao_cart_num
      A_center(1) = nucl_coord( ao_cart_nucl(j), 1 )
      A_center(2) = nucl_coord( ao_cart_nucl(j), 2 )
      A_center(3) = nucl_coord( ao_cart_nucl(j), 3 )
      power_A(1)  = ao_cart_power( j, 1 )
      power_A(2)  = ao_cart_power( j, 2 )
      power_A(3)  = ao_cart_power( j, 3 )
      do i= 1,ao_cart_num
        ao_cart_overlap_abs(i,j)= 0.d0
        B_center(1) = nucl_coord( ao_cart_nucl(i), 1 )
        B_center(2) = nucl_coord( ao_cart_nucl(i), 2 )
        B_center(3) = nucl_coord( ao_cart_nucl(i), 3 )
        power_B(1)  = ao_cart_power( i, 1 )
        power_B(2)  = ao_cart_power( i, 2 )
        power_B(3)  = ao_cart_power( i, 3 )
        do n = 1,ao_cart_prim_num(j)
          alpha = ao_cart_expo_ordered_transp(n,j)
          do l = 1, ao_cart_prim_num(i)
            beta = ao_cart_expo_ordered_transp(l,i)
            call overlap_x_abs(A_center(1),B_center(1),alpha,beta,power_A(1),power_B(1),overlap_x,lower_exp_val,dx,dim1)
            call overlap_x_abs(A_center(2),B_center(2),alpha,beta,power_A(2),power_B(2),overlap_y,lower_exp_val,dx,dim1)
            call overlap_x_abs(A_center(3),B_center(3),alpha,beta,power_A(3),power_B(3),overlap_z,lower_exp_val,dx,dim1)
            ao_cart_overlap_abs(i,j) += abs(ao_cart_coef_normalized_ordered_transp(n,j) * ao_cart_coef_normalized_ordered_transp(l,i)) * overlap_x * overlap_y * overlap_z
          enddo
        enddo
      enddo
    enddo
 !$OMP END PARALLEL DO

  endif

END_PROVIDER

! ---


BEGIN_PROVIDER [ double precision, S_cart_half, (ao_cart_num,ao_cart_num)  ]
 implicit none
 BEGIN_DOC
 ! :math:`S^{1/2}`
 END_DOC

  integer :: i,j,k
  double precision, allocatable  :: U(:,:)
  double precision, allocatable  :: Vt(:,:)
  double precision, allocatable  :: D(:)

  allocate(U(ao_cart_num,ao_cart_num),Vt(ao_cart_num,ao_cart_num),D(ao_cart_num))

  call svd(ao_cart_overlap,size(ao_cart_overlap,1),U,size(U,1),D,Vt,size(Vt,1),ao_cart_num,ao_cart_num)

  do i=1,ao_cart_num
    D(i) = dsqrt(D(i))
    do j=1,ao_cart_num
      S_cart_half(j,i) = 0.d0
    enddo
  enddo

  do k=1,ao_cart_num
      do j=1,ao_cart_num
        do i=1,ao_cart_num
          S_cart_half(i,j) = S_cart_half(i,j) + U(i,k)*D(k)*Vt(k,j)
        enddo
      enddo
  enddo

  deallocate(U,Vt,D)

END_PROVIDER
