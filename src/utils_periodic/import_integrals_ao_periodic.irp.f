program print_integrals
  print *, 'Number of AOs?'
  read(*,*) ao_num
  TOUCH ao_num
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  double precision :: integral
  double precision, allocatable :: A(:,:), B(:,:)
  double precision :: tmp_re, tmp_im

  integer             :: n_integrals 
  integer(key_kind), allocatable   :: buffer_i(:) 
  real(integral_kind), allocatable :: buffer_values(:)

  call ezfio_set_ao_basis_ao_num(ao_num)

  allocate (A(ao_num,ao_num), B(ao_num,ao_num) )
  
  A(1,1) = huge(1.d0)
  iunit = getunitandopen('E.qp','r')
  read (iunit,*,end=9) A(1,1)
  9 continue
  close(iunit)
  if (A(1,1) /= huge(1.d0)) then
    call ezfio_set_nuclei_nuclear_repulsion(A(1,1))
    call ezfio_set_nuclei_io_nuclear_repulsion("Read")
  endif

  A = 0.d0
  B = 0.d0
  iunit = getunitandopen('T.qp','r')
  do 
    read (iunit,*,end=10) i,j, tmp_re, tmp_im
    A(i,j) = tmp_re
    B(i,j) = tmp_im
    if (i.ne.j) then
      A(j,i) =  tmp_re
      B(j,i) = -tmp_im
    endif
  enddo
  10 continue
  close(iunit)
  call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(A(1:ao_num, 1:ao_num))
  call ezfio_set_ao_one_e_ints_ao_integrals_kinetic_imag(B(1:ao_num, 1:ao_num))
  call ezfio_set_ao_one_e_ints_io_ao_integrals_kinetic("Read")

  A = 0.d0
  B = 0.d0
  iunit = getunitandopen('S.qp','r')
  do 
    read (iunit,*,end=11) i,j, tmp_re, tmp_im
    A(i,j) = tmp_re
    B(i,j) = tmp_im
    if (i.ne.j) then
      A(j,i) =  tmp_re
      B(j,i) = -tmp_im
    endif
  enddo
  11 continue
  close(iunit)
  call ezfio_set_ao_one_e_ints_ao_integrals_overlap(A(1:ao_num, 1:ao_num))
  call ezfio_set_ao_one_e_ints_ao_integrals_overlap_imag(B(1:ao_num, 1:ao_num))
  call ezfio_set_ao_one_e_ints_io_ao_integrals_overlap("Read")

  A = 0.d0
  B = 0.d0
  iunit = getunitandopen('P.qp','r')
  do 
    read (iunit,*,end=14) i,j, tmp_re, tmp_im
    A(i,j) = tmp_re
    B(i,j) = tmp_im
    if (i.ne.j) then
      A(j,i) =  tmp_re
      B(j,i) = -tmp_im
    endif
  enddo
  14 continue
  close(iunit)
  call ezfio_set_ao_one_e_ints_ao_integrals_pseudo(A(1:ao_num,1:ao_num))
  call ezfio_set_ao_one_e_ints_ao_integrals_pseudo_imag(B(1:ao_num,1:ao_num))
  call ezfio_set_ao_one_e_ints_io_ao_integrals_pseudo("Read")

  A = 0.d0
  B = 0.d0
  iunit = getunitandopen('V.qp','r')
  do 
    read (iunit,*,end=12) i,j, tmp_re, tmp_im
    A(i,j) = tmp_re
    B(i,j) = tmp_im
    if (i.ne.j) then
      A(j,i) =  tmp_re
      B(j,i) = -tmp_im
    endif
  enddo
  12 continue
  close(iunit)
  call ezfio_set_ao_one_e_ints_ao_integrals_n_e(A(1:ao_num, 1:ao_num))
  call ezfio_set_ao_one_e_ints_ao_integrals_n_e_imag(B(1:ao_num, 1:ao_num))
  call ezfio_set_ao_one_e_ints_io_ao_integrals_n_e("Read")

!  allocate(buffer_i(ao_num**3), buffer_values(ao_num**3))
!  iunit = getunitandopen('W.qp','r')
!  n_integrals=0
!  buffer_values = 0.d0
!  do 
!    read (iunit,*,end=13) i,j,k,l, integral
!    n_integrals += 1
!    call two_e_integrals_index(i, j, k, l, buffer_i(n_integrals) )
!    buffer_values(n_integrals) = integral
!    if (n_integrals == size(buffer_i)) then
!      call insert_into_ao_integrals_map(n_integrals,buffer_i,buffer_values)
!      n_integrals = 0
!    endif
!  enddo
!  13 continue
!  close(iunit)
!  
!  if (n_integrals > 0) then
!    call insert_into_ao_integrals_map(n_integrals,buffer_i,buffer_values)
!  endif
!
!  call map_sort(ao_integrals_map)
!  call map_unique(ao_integrals_map)
!
!  call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
!  call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read')
  
end
