
! ---

double precision function ao_value(i, r)

  BEGIN_DOC
  ! Returns the value of the i-th ao at point $\textbf{r}$
  END_DOC

  implicit none
  integer,          intent(in) :: i
  double precision, intent(in) :: r(3)
  double precision, allocatable :: tmp_array_cart(:)
  double precision, external :: ddot
  ! TODO if in the cartesian basis transformation matrix is identity 
  allocate(tmp_array_cart(ao_cart_num))
  call give_all_aos_cart_at_r(r, tmp_array_cart) 
  ao_value = ddot(ao_cart_num,ao_cart_to_ao_basis_mat_transp(1,i),1,tmp_array_cart,1)
end


! ---

subroutine give_all_aos_at_r(r, tmp_array)

  BEGIN_dOC
  !
  ! input  : r == r(1) = x and so on
  !
  ! output : tmp_array(i) = aos(i) evaluated in $\textbf{r}$
  !
  END_DOC

  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: tmp_array(ao_num)
  double precision, allocatable :: tmp_array_cart(:)
  allocate(tmp_array_cart(ao_cart_num))
  call give_all_aos_cart_at_r(r, tmp_array_cart)
  call ao_cart_to_ao_basis_vec(tmp_array_cart, tmp_array)
end

! ---

subroutine give_all_aos_and_grad_at_r(r, aos_array, aos_grad_array)

  BEGIN_DOC
  !
  ! input : r(1) ==> r(1) = x, r(2) = y, r(3) = z
  !
  ! output : 
  !
  ! * aos_array(i) = ao(i) evaluated at ro
  ! * aos_grad_array(1,i) = gradient X of the ao(i) evaluated at $\textbf{r}$
  !
  END_DOC

  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: aos_array(ao_num)
  double precision, intent(out) :: aos_grad_array(3,ao_num)
  double precision, allocatable :: aos_cart_array(:), aos_cart_grad_array(:,:)
  allocate(aos_cart_array(ao_cart_num), aos_cart_grad_array(3,ao_cart_num))
  call give_all_aos_cart_and_grad_at_r(r, aos_cart_array, aos_cart_grad_array)
  call ao_cart_to_ao_basis_vec(aos_cart_array, aos_array)

  call dgemm('N','T',3,ao_num,ao_cart_num,1.d0, &
             aos_cart_grad_array, size(aos_cart_grad_array,1), &
             ao_cart_to_ao_basis_mat,size(ao_cart_to_ao_basis_mat,1), 0.d0,&
             aos_grad_array, size(aos_grad_array,1))



end

! ---

subroutine give_all_aos_and_grad_and_lapl_at_r(r, aos_array, aos_grad_array, aos_lapl_array)

  BEGIN_DOC
  !
  ! input  : r(1) ==> r(1) = x, r(2) = y, r(3) = z
  !
  ! output :
  !
  ! * aos_array(i) = ao(i) evaluated at $\textbf{r}$
  !
  ! * aos_grad_array(1,i) = $\nabla_x$ of the ao(i) evaluated at $\textbf{r}$
  !
  ! * aos_lapl_array(1,i) = $d/dx^2$ of the ao(i) evaluated at $\textbf{r}$
  END_DOC

  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: aos_array(ao_num)
  double precision, intent(out) :: aos_grad_array(3,ao_num)
  double precision, intent(out) :: aos_lapl_array(3,ao_num)
  double precision, allocatable :: aos_cart_array(:), aos_cart_grad_array(:,:), aos_cart_lapl_array(:,:)
  allocate(aos_cart_array(ao_cart_num), aos_cart_grad_array(3,ao_cart_num))
  call give_all_aos_cart_and_grad_and_lapl_at_r(r, aos_cart_array, aos_cart_grad_array, aos_cart_lapl_array)
  call ao_cart_to_ao_basis_vec(aos_cart_array, aos_array)

  call dgemm('N','T',3,ao_num,ao_cart_num,1.d0, &
             aos_cart_grad_array, size(aos_cart_grad_array,1), &
             ao_cart_to_ao_basis_mat,size(ao_cart_to_ao_basis_mat,1), 0.d0,&
             aos_grad_array, size(aos_grad_array,1))

  call dgemm('N','T',3,ao_num,ao_cart_num,1.d0, &
             aos_cart_lapl_array, size(aos_cart_lapl_array,1), &
             ao_cart_to_ao_basis_mat,size(ao_cart_to_ao_basis_mat,1), 0.d0,&
             aos_lapl_array, size(aos_lapl_array,1))

end

! ---

