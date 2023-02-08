program test_natorb
  implicit none
  BEGIN_DOC
! TODO : Reads psi_det in the EZFIO folder and prints out the left- and right-eigenvectors together with the energy. Saves the left-right wave functions at the end. 
  END_DOC
  print *, 'Hello world'
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 call routine
! call test

end

subroutine routine
 implicit none
 double precision, allocatable :: fock_diag(:),eigval(:),leigvec(:,:),reigvec(:,:),mat_ref(:,:)
 allocate(eigval(mo_num),leigvec(mo_num,mo_num),reigvec(mo_num,mo_num),fock_diag(mo_num),mat_ref(mo_num, mo_num))
 double precision, allocatable :: eigval_ref(:),leigvec_ref(:,:),reigvec_ref(:,:)
 allocate(eigval_ref(mo_num),leigvec_ref(mo_num,mo_num),reigvec_ref(mo_num,mo_num))
 
 double precision :: thr_deg
 integer :: i,n_real,j
 print*,'fock_matrix'
 do i = 1, mo_num
  fock_diag(i) = Fock_matrix_mo(i,i)
  print*,i,fock_diag(i)
 enddo
 thr_deg = 1.d-6
 mat_ref = -one_e_dm_mo
 print*,'diagonalization by block'
 call diag_mat_per_fock_degen(fock_diag,mat_ref,mo_num,thr_deg,leigvec,reigvec,eigval)
 call non_hrmt_bieig( mo_num, mat_ref&
                     , leigvec_ref, reigvec_ref& 
                     , n_real, eigval_ref)
 print*,'TEST ***********************************'
 double precision :: accu_l, accu_r
 do i = 1, mo_num
  accu_l = 0.d0
  accu_r = 0.d0
  do j = 1, mo_num
   accu_r += reigvec_ref(j,i) * reigvec(j,i)
   accu_l += leigvec_ref(j,i) * leigvec(j,i)
  enddo
  print*,i
  write(*,'(I3,X,100(F16.10,X))')i,eigval(i),eigval_ref(i),accu_l,accu_r
 enddo
end
