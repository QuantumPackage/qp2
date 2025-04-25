BEGIN_PROVIDER [ integer, basis_prod_num ]
  implicit none
  BEGIN_DOC
  ! Maximum number of basis functions for the projector
  END_DOC
  basis_prod_num = ao_num * (ao_num+1) / 2
END_PROVIDER

BEGIN_PROVIDER [ integer, basis_prod_idx, (2,basis_prod_num) ]
  implicit none
  BEGIN_DOC
  ! Indices of main basis to create projector basis
  END_DOC
  integer :: i,j,a
  a=0
  do j=1,ao_num
     do i=1,j
        a = a+1
        basis_prod_idx(1,a) = i
        basis_prod_idx(2,a) = j
     enddo
  enddo
END_PROVIDER

double precision function general_overlap_integral(dim,            &
     P_new,P_center,fact_p,p,p_inv,iorder_p,                       &
     Q_new,Q_center,fact_q,q,q_inv,iorder_q)
  implicit none
  BEGIN_DOC
  ! Computes the overlap integral <p|q> where p,q are products of primitives
  END_DOC
  integer,intent(in)             :: dim
  include 'utils/constants.include.F'
  double precision, intent(in)   :: P_new(0:max_dim,3),P_center(3),fact_p,p,p_inv
  double precision, intent(in)   :: Q_new(0:max_dim,3),Q_center(3),fact_q,q,q_inv
  integer, intent(in)            :: iorder_p(3)
  integer, intent(in)            :: iorder_q(3)

  double precision :: overlap, tmp
  double precision, external :: overlap_gaussian_x
  integer :: i,j,k

  general_overlap_integral = 1.d0
  do k=1,3
     tmp = 0.d0
     do i=0,iorder_p(k)
        do j=0,iorder_q(k)
           tmp += P_new(i,k)* Q_new(j,k) * overlap_gaussian_x(P_center(k), Q_center(k), p, q, i, j, dim)
        enddo
     enddo
     general_overlap_integral *= tmp
     if (tmp == 0.d0) exit
  enddo
  general_overlap_integral *= fact_p * fact_q
end function general_overlap_integral

double precision function basis_prod_overlap_func(a,b)
  implicit none
  BEGIN_DOC
  ! Overlap matrix element of the basis of products
  END_DOC
  include 'constants.include.F'

  integer, intent(in)            :: a,b
  double precision, external     :: general_overlap_integral
  double precision, external     :: ao_two_e_integral_general

  integer :: i,j,k,l
  i = basis_prod_idx(1,a)
  k = basis_prod_idx(2,a)
  j = basis_prod_idx(1,b)
  l = basis_prod_idx(2,b)
  basis_prod_overlap_func = ao_two_e_integral_general(i,j,k,l,general_overlap_integral)

end function basis_prod_overlap_func

 BEGIN_PROVIDER [ double precision, basis_prod_overlap, (basis_prod_num, basis_prod_num) ]
&BEGIN_PROVIDER [ double precision, basis_prod_overlap_inv, (basis_prod_num, basis_prod_num) ]
  implicit none
  BEGIN_DOC
  ! Overlap matrix of the basis of products
  END_DOC
  integer                        :: a,b
  double precision, external     :: basis_prod_overlap_func

  !$OMP PARALLEL DO PRIVATE(a,b) SCHEDULE(guided)
  do a=1,basis_prod_num
     do b=1,a
        basis_prod_overlap(a,b) = basis_prod_overlap_func(a,b)
        basis_prod_overlap(b,a) = basis_prod_overlap(a,b)
     enddo
     print *, a, '/', basis_prod_num
  enddo
  !$OMP END PARALLEL DO

  ! Eliminate linear dependencies

  double precision, allocatable :: U(:,:), D(:), Vt(:,:)
  allocate(U(basis_prod_num,basis_prod_num),Vt(basis_prod_num,basis_prod_num), D(basis_prod_num))

print *, 'svd'
  call svd( &
       basis_prod_overlap, size(basis_prod_overlap,1), &
       U,  size(U,1), D, &
       Vt, size(Vt,1), &
       basis_prod_num,basis_prod_num)

  double precision :: local_cutoff
  integer          :: i, j, k, mm

  local_cutoff = 1.d-10
  mm=basis_prod_num
  do i=1,basis_prod_num
    if ( D(i) < local_cutoff) then
      mm = mm-1
      D(i) = 0.d0
    endif
    do j=1,basis_prod_num
      basis_prod_overlap(j,i) = 0.d0
      basis_prod_overlap_inv(j,i) = 0.d0
    enddo
  enddo

  if (mm < basis_prod_num) then
     print *,  'Removed Linear dependencies below ', local_cutoff
  endif

print *, 'gemm1'
  double precision :: dk, dk_inv
  double precision, allocatable :: U2(:,:)
  allocate(U2(basis_prod_num,basis_prod_num))
  do k=1,basis_prod_num
     dk = D(k)
     do i=1,basis_prod_num
        U2(i,k) = U(i,k)*dk
     enddo
  enddo

  call dgemm('N','N', basis_prod_num, basis_prod_num, basis_prod_num, 1.d0, &
       U2, basis_prod_num, Vt, basis_prod_num, 0.d0, basis_prod_overlap, basis_prod_num)

print *, 'gemm2'
  do k=1,basis_prod_num
     if (D(k) == 0.d0) then
       dk = 0.d0
     else
       dk = 1.d0/D(k)
    endif
     do i=1,basis_prod_num
        U2(i,k) = U(i,k)*dk
     enddo
  enddo

  call dgemm('N','N', basis_prod_num, basis_prod_num, basis_prod_num, 1.d0, &
       U2, basis_prod_num, Vt, basis_prod_num, 0.d0, basis_prod_overlap_inv, basis_prod_num)


  deallocate(U,U2,D,Vt)
END_PROVIDER

