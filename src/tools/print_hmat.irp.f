program print_h_mat
  implicit none
  BEGIN_DOC
  ! program that prints out the CI matrix in sparse form 
  END_DOC
  read_wf = .True. 
  touch read_wf 
  call print_wf_dets
  call print_wf_coef
  call sparse_mat
  call full_mat
  call test_sparse_mat
end

subroutine print_wf_dets
 implicit none
 integer :: i,j
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.wf_det'
 i_unit_output = getUnitAndOpen(output,'w')
 write(i_unit_output,*)N_det,N_int
 do i = 1, N_det
  write(i_unit_output,*)psi_det_sorted(1:N_int,1,i)
  write(i_unit_output,*)psi_det_sorted(1:N_int,2,i)
 enddo
end

subroutine print_wf_coef
 implicit none
 integer :: i,j
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.wf_coef'
 i_unit_output = getUnitAndOpen(output,'w')
 write(i_unit_output,*)N_det,N_states
 do i = 1, N_det
  write(i_unit_output,*)psi_coef_sorted(i,1:N_states)
 enddo
end

subroutine sparse_mat
 implicit none
 integer :: i,j
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.hmat_sparse'
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, N_det
  write(i_unit_output,*)i,n_connected_per_det(i)
  do j =1, n_connected_per_det(i)
   write(i_unit_output,*)list_connected_det_per_det(j,i),sparse_h_mat(j,i)
  enddo
 enddo
end


subroutine full_mat
 implicit none
 integer :: i,j
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.hmat_full'
 i_unit_output = getUnitAndOpen(output,'w')
 do i = 1, N_det
  do j = i, N_det
   write(i_unit_output,*)i,j,H_matrix_all_dets(j,i)
  enddo
 enddo
end


subroutine test_sparse_mat
 implicit none
 integer :: i,j
 double precision, allocatable :: eigvec(:,:), eigval(:), hmat(:,:)
 allocate(eigval(N_det), eigvec(N_det,N_det),hmat(N_det,N_det))
 hmat = 0.d0
 do i = 1, N_det
  do j =1, n_connected_per_det(i)
   hmat(list_connected_det_per_det(j,i),i) = sparse_h_mat(j,i)
  enddo
 enddo
 call lapack_diag(eigval,eigvec,hmat,N_det,N_det) 
 print*,'The two energies should be the same '
 print*,'eigval(1) = ',eigval(1)
 print*,'psi_energy= ',CI_electronic_energy(1)


end
