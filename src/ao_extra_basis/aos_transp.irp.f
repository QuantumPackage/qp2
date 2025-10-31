
! ---

BEGIN_PROVIDER [ integer, Nucl_ao_extras_transposed, (N_ao_extras_max,nucl_num)]

  BEGIN_DOC
  ! List of ao_extras attached on each atom
  END_DOC
 
  implicit none
  integer              :: i
  integer, allocatable :: nucl_tmp(:)

  allocate(nucl_tmp(nucl_num))
  nucl_tmp = 0
  do i = 1, ao_extra_num
    nucl_tmp(ao_extra_nucl(i)) += 1
    Nucl_ao_extras_transposed(nucl_tmp(ao_extra_nucl(i)),ao_extra_nucl(i)) = i
  enddo
  deallocate(nucl_tmp)

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_extra_expo_ordered_transp_per_nucl, (ao_extra_prim_num_max,N_ao_extras_max,nucl_num) ]
 implicit none
 integer :: i,j,k,l
 do i = 1, nucl_num
  do j = 1,Nucl_N_ao_extras(i)
   k = Nucl_ao_extras_transposed(j,i)
   do l = 1, ao_extra_prim_num(k)
    ao_extra_expo_ordered_transp_per_nucl(l,j,i) = ao_extra_expo_ordered_transp(l,k)
   enddo
  enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [ integer, ao_extra_power_ordered_transp_per_nucl, (3,N_ao_extras_max,nucl_num) ]
 implicit none
 integer :: i,j,k,l
 do i = 1, nucl_num
  do j = 1,Nucl_N_ao_extras(i)
   k = Nucl_ao_extras_transposed(j,i)
   do l = 1, 3
    ao_extra_power_ordered_transp_per_nucl(l,j,i) = ao_extra_power(k,l)
   enddo
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_extra_coef_normalized_ordered_transp_per_nucl, (ao_extra_prim_num_max,N_ao_extras_max,nucl_num) ]
 implicit none
 integer :: i,j,k,l
 do i = 1, nucl_num
  do j = 1,Nucl_N_ao_extras(i)
   k = Nucl_ao_extras_transposed(j,i)
   do l = 1, ao_extra_prim_num(k)
    ao_extra_coef_normalized_ordered_transp_per_nucl(l,j,i) = ao_extra_coef_normalized_ordered_transp(l,k)
   enddo
  enddo
 enddo

END_PROVIDER

