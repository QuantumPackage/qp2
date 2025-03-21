program test_extra_basis
 implicit none
 integer :: i, j, k, l
 double precision :: ao_two_e_integral_mixed_direct
 print*,'System B nuclei information'
 do i = 1, extra_nucl_num
  print*,'charge',extra_nucl_charge(i)
  print*,'position'
  print*,extra_nucl_coord_transp(1:3,i)
 enddo
 print*,'System A nuclei information'
 do i = 1, nucl_num
  print*,'charge',nucl_charge(i)
  print*,'position'
  print*,nucl_coord_transp(1:3,i)
 enddo
 print*,'vne from B on the basis functions of A'
 do i = 1, ao_num
  do j = 1, ao_num
   print*,j,i,pot_vne_extra_basis(j,i)
  enddo
 enddo
 print*,'vne from A on the basis functions of B'
 do i = 1, ao_extra_num
  do j = 1, ao_extra_num
   print*,j,i,pot_vne_A_extra_basis(j,i)
  enddo
 enddo
 print*,'Density matrix from system A'
 do i = 1, ao_num
  do j = 1, ao_num
   print*,j,i,one_e_dm_ao(j,i)
  enddo
 enddo
 print*,'Density matrix from system B'
 output=trim(ezfio_filename)//'.one_e_dm_b'
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, ao_extra_num
  do j = 1, ao_extra_num
   write(i_unit_output,'(2(I3,X),F16.10)')j,i,ao_extra_one_e_dm(j,i,1)
  enddo
 enddo
 print*,'Two electron integrals between A and B'
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.two_e_ints'
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_extra_num
    do l = 1, ao_extra_num
     write(i_unit_output,'(4(I3,X),F16.10)')i,j,k,l,ao_two_e_integral_mixed_direct(i, j, k, l)
    enddo
   enddo
  enddo
 enddo

end
