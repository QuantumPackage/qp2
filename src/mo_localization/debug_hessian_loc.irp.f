program debug_hessian_loc

  !BEGIN_DOC
  ! Check if the hessian is correct
  !END_DOC
  
  implicit none

  integer :: list_size, n
  integer, allocatable :: list(:)
  double precision, allocatable :: H(:,:), H2(:,:)
  double precision :: threshold, max_error, max_elem
  integer :: i, nb_error

  threshold = 1d-12

  list = list_act
  list_size = dim_list_act_orb

  n = list_size*(list_size-1)/2
  
  allocate(H(n,n),H2(n,n))

  if (localization_method == 'boys') then
    print*,'Foster-Boys'
    call hessian_FB(n,list_size,list,H)
    call hessian_FB_omp(n,list_size,list,H2)
  elseif(localization_method == 'pipek') then
    print*,'Pipek-Mezey'
    call hessian_PM(n,list_size,list,H)
    call hessian_PM(n,list_size,list,H2)
  else
    print*,'Unknown localization_method, please select boys or pipek'
    call abort
  endif
 
  do i = 1, n
    print*,i,H(i,i)
  enddo

  H = H - H2

  nb_error = 0
  max_elem = 0d0

  do i = 1, n
    if (dabs(H(i,i)) > threshold) then
      print*,H(i,i)
      nb_error = nb_error + 1
      if (dabs(H(i,i)) > max_elem) then
        max_elem = H(i,i)
      endif
    endif
  enddo

  print*,'Threshold error', threshold
  print*, 'Nb error', nb_error
  print*,'Max error', max_elem

  deallocate(H,H2)
  
end
