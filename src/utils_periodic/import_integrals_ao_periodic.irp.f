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

  integer             :: n_integrals_1, n_integrals_2 
  integer(key_kind), allocatable   :: buffer_i_1(:), buffer_i_2(:) 
  real(integral_kind), allocatable :: buffer_values_1(:), buffer_values_2(:)
  logical :: use_map1
  integer(key_kind) :: idx_tmp
  double precision :: sign


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

  allocate(buffer_i_1(ao_num**3), buffer_values_1(ao_num**3))
  allocate(buffer_i_2(ao_num**3), buffer_values_2(ao_num**3))
  iunit = getunitandopen('W.qp','r')
  n_integrals_1=0
  n_integrals_2=0
  buffer_values_1 = 0.d0
  buffer_values_2 = 0.d0
  do 
    read (iunit,*,end=13) i,j,k,l, tmp_re, tmp_im
    call ao_two_e_integral_periodic_map_idx_sign(i,j,k,l,use_map1,idx_tmp,sign)
    print'(4(I4),(L3),(I6),(F7.1))',i,j,k,l,use_map1,idx_tmp,sign
    if (use_map1) then
      n_integrals_1 += 1
      buffer_i_1(n_integrals_1)=idx_tmp
      buffer_values_1(n_integrals_1)=tmp_re
      print'(A,4(I4),(I6),(E15.7))','map1',i,j,k,l,idx_tmp,tmp_re
      if (sign.ne.0.d0) then 
        n_integrals_1 += 1
        buffer_i_1(n_integrals_1)=idx_tmp+1
        buffer_values_1(n_integrals_1)=tmp_im*sign
        print'(A,4(I4),(I6),(E15.7))','map1',i,j,k,l,idx_tmp+1,tmp_im*sign
      endif 
      if (n_integrals_1 >= size(buffer_i_1)-1) then
        call insert_into_ao_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1)
        n_integrals_1 = 0
      endif
    else
      n_integrals_2 += 1
      buffer_i_2(n_integrals_2)=idx_tmp
      buffer_values_2(n_integrals_2)=tmp_re
      print'(A,4(I4),(I6),(E15.7))','map2',i,j,k,l,idx_tmp,tmp_re
      if (sign.ne.0.d0) then
        n_integrals_2 += 1
        buffer_i_2(n_integrals_2)=idx_tmp+1
        buffer_values_2(n_integrals_2)=tmp_im*sign
        print'(A,4(I4),(I6),(E15.7))','map2',i,j,k,l,idx_tmp+1,tmp_im*sign
      endif 
      if (n_integrals_2 >= size(buffer_i_2)-1) then
        call insert_into_ao_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2)
        n_integrals_2 = 0
      endif
    endif
  enddo
  13 continue
  close(iunit)
  
  if (n_integrals_1 > 0) then
    call insert_into_ao_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1)
  endif
  if (n_integrals_2 > 0) then
    call insert_into_ao_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2)
  endif

  call map_sort(ao_integrals_map)
  call map_unique(ao_integrals_map)
  call map_sort(ao_integrals_map_2)
  call map_unique(ao_integrals_map_2)

  call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_periodic_1',ao_integrals_map)
  call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_periodic_2',ao_integrals_map_2)
  call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read')
  
end
