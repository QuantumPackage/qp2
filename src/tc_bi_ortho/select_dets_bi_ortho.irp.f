program tc_bi_ortho
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  !!!!!!!!!!!!!!! WARNING NO 3-BODY
  !!!!!!!!!!!!!!! WARNING NO 3-BODY
  three_body_h_tc = .False.
  touch three_body_h_tc 
  !!!!!!!!!!!!!!! WARNING NO 3-BODY
  !!!!!!!!!!!!!!! WARNING NO 3-BODY

  call routine_test
! call test
end

subroutine routine_test
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,n_good,degree
 integer(bit_kind), allocatable :: dets(:,:,:)
 integer, allocatable :: iorder(:)
 double precision, allocatable :: coef(:),coef_new(:,:)
 double precision :: thr
 allocate(coef(N_det), iorder(N_det))
 do i = 1, N_det
  iorder(i) = i
  call get_excitation_degree(HF_bitmask,psi_det(1,1,i),degree,N_int)
  if(degree==1)then
   coef(i) = -0.5d0
  else
   coef(i) = -dabs(coef_pt1_bi_ortho(i))
  endif
 enddo
 call dsort(coef,iorder,N_det)
 !thr = save_threshold
 thr = 1d-15
 n_good = 0
 do i = 1, N_det
  if(dabs(coef(i)).gt.thr)then
   n_good += 1
  endif
 enddo
 print*,'n_good = ',n_good
 allocate(dets(N_int,2,n_good),coef_new(n_good,n_states))
 do i = 1, n_good
  dets(:,:,i) = psi_det(:,:,iorder(i))
  coef_new(i,:) = psi_coef(iorder(i),:)
 enddo
 call save_wavefunction_general(n_good,n_states,dets,n_good,coef_new)


end
