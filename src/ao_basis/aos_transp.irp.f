
! ---

BEGIN_PROVIDER [ integer, Nucl_Aos_transposed, (N_AOs_max,nucl_num)]

  BEGIN_DOC
  ! List of AOs attached on each atom
  END_DOC
 
  implicit none
  integer              :: i
  integer, allocatable :: nucl_tmp(:)

  allocate(nucl_tmp(nucl_num))
  nucl_tmp = 0
  do i = 1, ao_num
    nucl_tmp(ao_nucl(i)) += 1
    Nucl_Aos_transposed(nucl_tmp(ao_nucl(i)),ao_nucl(i)) = i
  enddo
  deallocate(nucl_tmp)

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_expo_ordered_transp_per_nucl, (ao_prim_num_max,N_AOs_max,nucl_num) ]
 implicit none
 integer :: i,j,k,l
 do i = 1, nucl_num
  do j = 1,Nucl_N_Aos(i)
   k = Nucl_Aos_transposed(j,i)
   do l = 1, ao_prim_num(k)
    ao_expo_ordered_transp_per_nucl(l,j,i) = ao_expo_ordered_transp(l,k)
   enddo
  enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [ integer, ao_power_ordered_transp_per_nucl, (3,N_AOs_max,nucl_num) ]
 implicit none
 integer :: i,j,k,l
 do i = 1, nucl_num
  do j = 1,Nucl_N_Aos(i)
   k = Nucl_Aos_transposed(j,i)
   do l = 1, 3
    ao_power_ordered_transp_per_nucl(l,j,i) = ao_power(k,l)
   enddo
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_coef_normalized_ordered_transp_per_nucl, (ao_prim_num_max,N_AOs_max,nucl_num) ]
 implicit none
 integer :: i,j,k,l
 do i = 1, nucl_num
  do j = 1,Nucl_N_Aos(i)
   k = Nucl_Aos_transposed(j,i)
   do l = 1, ao_prim_num(k)
    ao_coef_normalized_ordered_transp_per_nucl(l,j,i) = ao_coef_normalized_ordered_transp(l,k)
   enddo
  enddo
 enddo

END_PROVIDER

