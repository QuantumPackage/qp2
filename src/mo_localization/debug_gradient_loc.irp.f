program debug_gradient_loc

  !BEGIN_DOC
  ! Check if the gradient is correct
  !END_DOC

  implicit none

  integer :: list_size, n
  integer, allocatable :: list(:)
  double precision, allocatable :: v_grad(:), v_grad2(:)
  double precision :: norm, max_elem, threshold, max_error
  integer :: i, nb_error

  threshold = 1d-12

  list = list_act
  list_size = dim_list_act_orb

  n = list_size*(list_size-1)/2
  
  allocate(v_grad(n),v_grad2(n))

  if (localization_method == 'boys') then
    print*,'Foster-Boys'
    call gradient_FB(n,list_size,list,v_grad,max_elem,norm)
    call gradient_FB_omp(n,list_size,list,v_grad2,max_elem,norm)
  elseif (localization_method == 'pipek') then
    print*,'Pipek-Mezey'
    call gradient_PM(n,list_size,list,v_grad,max_elem,norm)
    call gradient_PM(n,list_size,list,v_grad2,max_elem,norm) 
  else
    print*,'Unknown localization_method, please select boys or pipek'
    call abort
  endif
 
  do i = 1, n
    print*,i,v_grad(i)
  enddo

  v_grad = v_grad - v_grad2

  nb_error = 0
  max_elem = 0d0

  do i = 1, n
    if (dabs(v_grad(i)) > threshold) then
      print*,v_grad(i)
      nb_error = nb_error + 1
      if (dabs(v_grad(i)) > max_elem) then
        max_elem = v_grad(i)
      endif
    endif
  enddo

  print*,'Threshold error', threshold
  print*, 'Nb error', nb_error
  print*,'Max error', max_elem

  deallocate(v_grad,v_grad2)
 
end
