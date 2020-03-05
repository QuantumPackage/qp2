program print_integrals
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



provide ao_two_e_integrals_in_map
  allocate (A(ao_num,ao_num), B(ao_num,ao_num) )
  
  A(1,1) = huge(1.d0)
  iunit = getunitandopen('E.qp','r')
  read (iunit,*,end=9) A(1,1)
  9 continue
  close(iunit)
  if (A(1,1) /= huge(1.d0)) then
!    call ezfio_set_nuclei_nuclear_repulsion(A(1,1))
!    call ezfio_set_nuclei_io_nuclear_repulsion("Read")
     print*, nuclear_repulsion,A(1,1)
  endif

  A = 0.d0
  B = 0.d0
 ! iunit = getunitandopen('T.qp','r')
 ! do 
 !   read (iunit,*,end=10) i,j, tmp_re, tmp_im
 !   A(i,j) = tmp_re
 !   B(i,j) = tmp_im
 !   print*,ao_kinetic_integrals(i,j),A(i,j)
 !   print*,ao_kinetic_integrals_imag(i,j),B(i,j)
 !   if (i.ne.j) then
 !     A(j,i) =  tmp_re
 !     B(j,i) = -tmp_im
 !   print*,ao_kinetic_integrals(j,i),A(j,i)
 !   print*,ao_kinetic_integrals_imag(j,i),B(j,i)
 !   endif
 ! enddo
 ! 10 continue
 ! close(iunit)
!  call ezfio_set_ao_one_e_ints_ao_integrals_kinetic(A(1:ao_num, 1:ao_num))
!  call ezfio_set_ao_one_e_ints_ao_integrals_kinetic_imag(B(1:ao_num, 1:ao_num))
!  call ezfio_set_ao_one_e_ints_io_ao_integrals_kinetic("Read")

  A = 0.d0
  B = 0.d0
 ! iunit = getunitandopen('S.qp','r')
 ! do 
 !   read (iunit,*,end=11) i,j, tmp_re, tmp_im
 !   A(i,j) = tmp_re
 !   B(i,j) = tmp_im
 !   print*,real(ao_overlap_complex(i,j)),A(i,j)
 !   print*,imag(ao_overlap_complex(i,j)),B(i,j)
 !   print*,ao_overlap_imag(i,j),B(i,j)
 !   if (i.ne.j) then
 !     A(j,i) =  tmp_re
 !     B(j,i) = -tmp_im
 !   print*,real(ao_overlap_complex(j,i)),A(j,i)
 !   print*,imag(ao_overlap_complex(j,i)),B(j,i)
 !   print*,ao_overlap_imag(j,i),B(j,i)
 !   endif
 ! enddo
 ! 11 continue
 ! close(iunit)
!  call ezfio_set_ao_one_e_ints_ao_integrals_overlap(A(1:ao_num, 1:ao_num))
!  call ezfio_set_ao_one_e_ints_ao_integrals_overlap_imag(B(1:ao_num, 1:ao_num))
!  call ezfio_set_ao_one_e_ints_io_ao_integrals_overlap("Read")

  A = 0.d0
  B = 0.d0
!  iunit = getunitandopen('P.qp','r')
!  do 
!    read (iunit,*,end=14) i,j, tmp_re, tmp_im
!    A(i,j) = tmp_re
!    B(i,j) = tmp_im
!    print*,ao_pseudo_integrals(i,j),A(i,j)
!    print*,ao_pseudo_integrals_imag(i,j),B(i,j)
!   ! print*,real(ao_integrals_pseudo(i,j)),A(i,j)
!   ! print*,imag(ao_integrals_pseudo(i,j)),B(i,j)
!    if (i.ne.j) then
!      A(j,i) =  tmp_re
!      B(j,i) = -tmp_im
!    print*,ao_pseudo_integrals(j,i),A(j,i)
!    print*,ao_pseudo_integrals_imag(j,i),B(j,i)
!   ! print*,real(ao_integrals_pseudo(j,i)),A(j,i)
!   ! print*,imag(ao_integrals_pseudo(j,i)),B(j,i)
!    endif
!  enddo
!  14 continue
!  close(iunit)
!  call ezfio_set_ao_one_e_ints_ao_integrals_pseudo(A(1:ao_num,1:ao_num))
!  call ezfio_set_ao_one_e_ints_ao_integrals_pseudo_imag(B(1:ao_num,1:ao_num))
!  call ezfio_set_ao_one_e_ints_io_ao_integrals_pseudo("Read")

  A = 0.d0
  B = 0.d0
!  iunit = getunitandopen('V.qp','r')
!  do 
!    read (iunit,*,end=12) i,j, tmp_re, tmp_im
!    A(i,j) = tmp_re
!    B(i,j) = tmp_im
!    print*,ao_integrals_n_e(i,j),A(i,j)
!    print*,ao_integrals_n_e_imag(i,j),B(i,j)
!    if (i.ne.j) then
!      A(j,i) =  tmp_re
!      B(j,i) = -tmp_im
!      print*,ao_integrals_n_e(j,i),A(j,i)
!      print*,ao_integrals_n_e_imag(j,i),B(j,i)
!    endif
!  enddo
!  12 continue
!  close(iunit)
!  call ezfio_set_ao_one_e_ints_ao_integrals_n_e(A(1:ao_num, 1:ao_num))
!  call ezfio_set_ao_one_e_ints_ao_integrals_n_e_imag(B(1:ao_num, 1:ao_num))
!  call ezfio_set_ao_one_e_ints_io_ao_integrals_n_e("Read")
  complex*16 :: int2e_tmp1,int2e_tmp2,get_ao_two_e_integral_complex_simple,get_ao_two_e_integral_complex, tmp_cmplx
  double precision :: tmp3,tmp4,tmp5,tmp6
  double precision :: thr0
  thr0 = 1.d-10
  allocate(buffer_i_1(ao_num**3), buffer_values_1(ao_num**3))
  allocate(buffer_i_2(ao_num**3), buffer_values_2(ao_num**3))
  iunit = getunitandopen('W.qp','r')
  n_integrals_1=0
  n_integrals_2=0
  buffer_values_1 = 0.d0
  buffer_values_2 = 0.d0
  do 
    read (iunit,*,end=13) i,j,k,l, tmp_re, tmp_im
    tmp_cmplx = dcmplx(tmp_re,tmp_im)
    int2e_tmp1 = get_ao_two_e_integral_complex_simple(i,j,k,l,ao_integrals_map,ao_integrals_map_2)
    int2e_tmp2 = get_ao_two_e_integral_complex(i,j,k,l,ao_integrals_map,ao_integrals_map_2)
  !  print'(4(I4),3(E15.7))',i,j,k,l,tmp_re,real(int2e_tmp1),real(int2e_tmp2)
 !   print'(4(I4),3(E15.7))',i,j,k,l,tmp_im,imag(int2e_tmp1),imag(int2e_tmp2)
    call ao_two_e_integral_complex_map_idx_sign(i,j,k,l,use_map1,idx_tmp,sign)
!    print*,use_map1,idx_tmp,sign
    call map_get(ao_integrals_map,idx_tmp,tmp3)
    call map_get(ao_integrals_map_2,idx_tmp,tmp4)
    call map_get(ao_integrals_map,idx_tmp+1,tmp5)
    call map_get(ao_integrals_map_2,idx_tmp+1,tmp6)
   ! print*,tmp3,tmp4
   ! print*,tmp5,tmp6
   if (cdabs(tmp_cmplx-int2e_tmp1).gt.thr0) then
     print'(4(I4),4(E15.7))',i,j,k,l,tmp_cmplx,int2e_tmp1
   endif
    integer*8 :: ii
          ii = l-ao_integrals_cache_min
          ii = ior( shiftl(ii,6), k-ao_integrals_cache_min)
          ii = ior( shiftl(ii,6), j-ao_integrals_cache_min)
          ii = ior( shiftl(ii,6), i-ao_integrals_cache_min)
!    print*,'cache(pbc)=', ao_integrals_cache_complex(ii)
!    print*,'cache(old)=', ao_integrals_cache(ii)
!    print*
!    if (use_map1) then
!      n_integrals_1 += 1
!      buffer_i_1(n_integrals_1-1)=idx_tmp
!      buffer_values_1(n_integrals_1-1)=tmp_re
!      if (sign.ne.0.d0) then 
!        n_integrals_1 += 1
!        buffer_i_1(n_integrals_2)=idx_tmp+1
!        buffer_values_1(n_integrals_1)=tmp_im*sign
!      endif 
!      if (n_integrals_1 >= size(buffer_i_1)-1) then
!!        call insert_into_ao_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1)
!        n_integrals_1 = 0
!      endif
!    else
!      n_integrals_2 += 1
!      buffer_i_2(n_integrals_2-1)=idx_tmp
!      buffer_values_2(n_integrals_2-1)=tmp_re
!      if (sign.ne.0.d0) then
!        n_integrals_2 += 1
!        buffer_i_2(n_integrals_2)=idx_tmp+1
!        buffer_values_2(n_integrals_2)=tmp_im*sign
!      endif 
!      if (n_integrals_2 >= size(buffer_i_2)-1) then
!!        call insert_into_ao_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2)
!        n_integrals_2 = 0
!      endif
!    endif
  enddo
  13 continue
  close(iunit)
  
!  if (n_integrals_1 > 0) then
!!    call insert_into_ao_integrals_map(n_integrals_1,buffer_i_1,buffer_values_1)
!  endif
!  if (n_integrals_2 > 0) then
!!    call insert_into_ao_integrals_map_2(n_integrals_2,buffer_i_2,buffer_values_2)
!  endif
!
!  call map_sort(ao_integrals_map)
!  call map_unique(ao_integrals_map)
!  call map_sort(ao_integrals_map_2)
!  call map_unique(ao_integrals_map_2)
!
!  call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_complex_1',ao_integrals_map)
!  call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_complex_2',ao_integrals_map_2)
!  call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read'
print*,'map1'
  do i=0,ao_integrals_map%map_size
    print*,i,ao_integrals_map%map(i)%value(:)
    print*,i,ao_integrals_map%map(i)%key(:)
  enddo
print*,'map2'
  do i=0,ao_integrals_map_2%map_size
    print*,i,ao_integrals_map_2%map(i)%value(:)
    print*,i,ao_integrals_map_2%map(i)%key(:)
  enddo
end
