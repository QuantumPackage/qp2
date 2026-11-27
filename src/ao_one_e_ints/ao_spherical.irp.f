use sizet
BEGIN_PROVIDER [ double precision, ao_cart_to_sphe_coef_transp, (ao_sphe_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Transpose of :c:data:`ao_cart_to_sphe_coef` 
  END_DOC
  integer(size_t) :: i,j
  do i = 1, ao_num
    do j = 1, ao_sphe_num
      ao_cart_to_sphe_coef_transp(j,i) = ao_cart_to_sphe_coef(i,j)
    enddo
  enddo
END_PROVIDER


subroutine ao_cart_to_ao_sphe(A_ao_cart,LDA_ao_cart,A_ao_sphe,LDA_ao_sphe)
 use sizet
  implicit none
  BEGIN_DOC
  ! Transform matrix A from the |AO| cartesian basis to the |AO| spherical basis
  !
  ! :math:`(B^{sc})^T.A^c.B^{sc}`
  !
  ! where :math:`B^{sc}` is :c:data:`ao_cart_to_sphe_coef`, 
  ! the matrix of coefficients from the cartesian AO basis to spherical one,
  ! and :math:`B^{sc}` is :c:data:`ao_cart_to_sphe_coef_transp`, its transpose.
  END_DOC
  integer(size_t), intent(in)    :: LDA_ao_cart,LDA_ao_sphe
  double precision, intent(in)   :: A_ao_cart(LDA_ao_cart,ao_num)
  double precision, intent(out)  :: A_ao_sphe(LDA_ao_sphe,ao_sphe_num)
  double precision, allocatable  :: T(:,:)
  !
  allocate (T(ao_num,ao_sphe_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call dgemm('N','N', int(ao_num,4), int(ao_sphe_num,4), int(ao_num,4), &
      1.d0, A_ao_cart, int(LDA_ao_cart,4), &
      ao_cart_to_sphe_coef, int(size(ao_cart_to_sphe_coef,1),4), &
      0.d0, T, int(size(T,1),4))
  ! Notice that for the following dgemm we could have used
  ! ao_cart_to_sphe_coef_transp, but instead we transposed with the 'T' argument
  call dgemm('T','N', int(ao_sphe_num,4), int(ao_sphe_num,4), int(ao_num,4), &
      1.d0, ao_cart_to_sphe_coef, int(size(ao_cart_to_sphe_coef,1),4), &
      T, int(ao_num,4), 0.d0, A_ao_sphe, int(size(A_ao_sphe,1),4))
  !
  ! call restore_symmetry(mo_num,mo_num,A_ao_sphe,size(A_ao_sphe,1),1.d-15)
  deallocate(T)
end


BEGIN_PROVIDER [ double precision, ao_cart_to_sphe_overlap, (ao_sphe_num,ao_sphe_num) ]
  implicit none
  BEGIN_DOC
  ! |AO| overlap matrix in the spherical basis set obtained as
  !
  ! :math:`(B^{sc})^T.S^c.B^{sc}`
  !
  ! where :math:`S^c` is the overlap matrix in the cartesian AO basis
  END_DOC
  !
  call ao_cart_to_ao_sphe(ao_overlap, int(ao_num, size_t), &
      ao_cart_to_sphe_overlap, int(ao_sphe_num, size_t))
  !
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_cart_to_sphe_overlap_inv, (ao_sphe_num,ao_sphe_num) ]
 implicit none
 BEGIN_DOC
 ! Inverse of :c:data:`ao_cart_to_sphe_overlap`
 END_DOC
 !
  call get_pseudo_inverse(ao_cart_to_sphe_overlap, &
                         int(size(ao_cart_to_sphe_overlap,1),size_t), &
                         int(ao_sphe_num,size_t), int(ao_sphe_num,size_t), &
                         ao_cart_to_sphe_overlap_inv, &
                         int(size(ao_cart_to_sphe_overlap_inv,1),size_t), &
                         lin_dep_cutoff)
END_PROVIDER
