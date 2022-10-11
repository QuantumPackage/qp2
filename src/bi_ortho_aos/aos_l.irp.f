 BEGIN_PROVIDER [ double precision, ao_coef_l , (ao_num,ao_prim_num_max) ]
 implicit none
 BEGIN_DOC
! Primitive coefficients and exponents for each atomic orbital. Copied from shell info.
 END_DOC

 integer :: i, l
 do i=1,ao_num
   l = ao_shell(i)
   ao_coef_l(i,:) = shell_coef(l,:)
 end do
END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_coef_l_normalized, (ao_num,ao_prim_num_max) ]
&BEGIN_PROVIDER [ double precision, ao_coef_l_normalization_factor, (ao_num) ]
  implicit none
  BEGIN_DOC
  ! Coefficients including the |AO| normalization
  END_DOC

  do i=1,ao_num
    l = ao_shell(i)
    ao_coef_l_normalized(i,:) = shell_coef(l,:) * shell_normalization_factor(l)
  end do

  double precision               :: norm,overlap_x,overlap_y,overlap_z,C_A(3), c
  integer                        :: l, powA(3), nz
  integer                        :: i,j,k
  nz=100
  C_A = 0.d0

  do i=1,ao_num

    powA(1) = ao_power(i,1)
    powA(2) = ao_power(i,2)
    powA(3) = ao_power(i,3)

    ! Normalization of the primitives
    if (primitives_normalized) then
      do j=1,ao_prim_num(i)
        call overlap_gaussian_xyz(C_A,C_A,ao_expo(i,j),ao_expo(i,j), &
           powA,powA,overlap_x,overlap_y,overlap_z,norm,nz)
        ao_coef_l_normalized(i,j) = ao_coef_l_normalized(i,j)/dsqrt(norm)
      enddo
    endif
    ! Normalization of the contracted basis functions
    if (ao_normalized) then
      norm = 0.d0
      do j=1,ao_prim_num(i)
        do k=1,ao_prim_num(i)
          call overlap_gaussian_xyz(C_A,C_A,ao_expo(i,j),ao_expo(i,k),powA,powA,overlap_x,overlap_y,overlap_z,c,nz)
          norm = norm+c*ao_coef_l_normalized(i,j)*ao_coef_l_normalized(i,k)
        enddo
      enddo
      ao_coef_l_normalization_factor(i) = 1.d0/dsqrt(norm)
    else
      ao_coef_l_normalization_factor(i) = 1.d0
    endif
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_coef_l_normalized_ordered, (ao_num,ao_prim_num_max) ]
  implicit none
  BEGIN_DOC
  ! Sorted primitives to accelerate 4 index |MO| transformation
  END_DOC

  integer                        :: iorder(ao_prim_num_max)
  double precision               :: d(ao_prim_num_max,2)
  integer                        :: i,j
  do i=1,ao_num
    do j=1,ao_prim_num(i)
      iorder(j) = j
      d(j,2) = ao_coef_l_normalized(i,j)
    enddo
    call dsort(d(1,1),iorder,ao_prim_num(i))
    call dset_order(d(1,2),iorder,ao_prim_num(i))
    do j=1,ao_prim_num(i)
      ao_coef_l_normalized_ordered(i,j) = d(j,2)
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_coef_l_normalized_ordered_transp, (ao_prim_num_max,ao_num) ]
  implicit none
  BEGIN_DOC
  ! Transposed :c:data:`ao_coef_l_normalized_ordered`
  END_DOC
  integer                        :: i,j
  do j=1, ao_num
    do i=1, ao_prim_num_max
      ao_coef_l_normalized_ordered_transp(i,j) = ao_coef_l_normalized_ordered(j,i)
    enddo
  enddo
END_PROVIDER

