subroutine routine_save_rotated_mos
 implicit none
 integer :: i,j,k,n_degen_list,m,n,n_degen,ilast,ifirst
 double precision, allocatable :: mo_r_coef_good(:,:),mo_l_coef_good(:,:)
 allocate(mo_l_coef_good(ao_num, mo_num), mo_r_coef_good(ao_num,mo_num))
 double precision, allocatable :: mo_r_coef_new(:,:)
 double precision :: norm
 mo_r_coef_good = mo_r_coef
 mo_l_coef_good = mo_l_coef
 allocate(mo_r_coef_new(ao_num, mo_num))
 mo_r_coef_new = mo_r_coef
 do i = 1, mo_num
  norm = 1.d0/dsqrt(overlap_mo_r(i,i))
  do j = 1, ao_num
   mo_r_coef_new(j,i) *= norm
  enddo
 enddo
 double precision, allocatable :: fock_diag(:),s_mat(:,:)
 integer, allocatable :: list_degen(:,:)
 allocate(list_degen(2,mo_num),s_mat(mo_num,mo_num),fock_diag(mo_num))
 do i = 1, mo_num
  fock_diag(i) = fock_matrix_mo(i,i)
 enddo
 ! compute the overlap between the left and rescaled right
 call build_s_matrix(ao_num,mo_num,mo_r_coef_new,mo_r_coef_new,ao_overlap,s_mat)
 call give_degen(fock_diag,mo_num,thr_degen_tc,list_degen,n_degen_list)
 print*,'fock_matrix_mo'
 do i = 1, mo_num
  print*,i,fock_diag(i),angle_left_right(i)
 enddo
 print*,'Overlap '
 do i = 1, mo_num
  write(*,'(I2,X,100(F8.4,X))')i,s_mat(:,i)
 enddo
   
 do i = 1, n_degen_list
  ifirst = list_degen(1,i)
  ilast  = list_degen(2,i)
  n_degen = ilast - ifirst +1
  print*,'ifirst,n_degen = ',ifirst,n_degen
  double precision, allocatable :: stmp(:,:),T(:,:),Snew(:,:),smat2(:,:)
  double precision, allocatable :: mo_l_coef_tmp(:,:),mo_r_coef_tmp(:,:),mo_l_coef_new(:,:)
  allocate(stmp(n_degen,n_degen),smat2(n_degen,n_degen))
  allocate(mo_r_coef_tmp(ao_num,n_degen),mo_l_coef_tmp(ao_num,n_degen),mo_l_coef_new(ao_num,n_degen))
  allocate(T(n_degen,n_degen),Snew(n_degen,n_degen))
  do j = 1, n_degen
   mo_r_coef_tmp(1:ao_num,j) = mo_r_coef_new(1:ao_num,j+ifirst-1)
   mo_l_coef_tmp(1:ao_num,j) = mo_l_coef(1:ao_num,j+ifirst-1)
  enddo
  ! Orthogonalization of right functions
  print*,'Orthogonalization of right functions'
  call orthog_functions(ao_num,n_degen,mo_r_coef_tmp,ao_overlap)
  ! Orthogonalization of left functions
  print*,'Orthogonalization of left functions'
  call orthog_functions(ao_num,n_degen,mo_r_coef_tmp,ao_overlap)
  print*,'Overlap lef-right '
  call build_s_matrix(ao_num,n_degen,mo_r_coef_tmp,mo_l_coef_tmp,ao_overlap,stmp)
  do j = 1, n_degen
   write(*,'(100(F8.4,X))')stmp(:,j)
  enddo
  if(maxovl_tc)then
   T    = 0.d0
   Snew = 0.d0
   call maxovl(n_degen, n_degen, stmp, T, Snew)
   print*,'overlap after'
   do j = 1, n_degen
    write(*,'(100(F16.10,X))')Snew(:,j)
   enddo
   call dgemm( 'N', 'N', ao_num, n_degen, n_degen, 1.d0               &
             , mo_l_coef_tmp, size(mo_l_coef_tmp, 1), T(1,1), size(T, 1) &
             , 0.d0, mo_l_coef_new, size(mo_l_coef_new, 1) )
   call build_s_matrix(ao_num,n_degen,mo_l_coef_new,mo_r_coef_tmp,ao_overlap,stmp)
   print*,'Overlap test'
   do j = 1, n_degen
    write(*,'(100(F16.10,X))')stmp(:,j)
   enddo
  else 
   mo_l_coef_new = mo_l_coef_tmp
  endif
  call impose_biorthog_svd_overlap(ao_num, n_degen, ao_overlap, mo_l_coef_new, mo_r_coef_tmp)
  call build_s_matrix(ao_num,n_degen,mo_l_coef_new,mo_r_coef_tmp,ao_overlap,stmp)
  print*,'LAST OVERLAP '
  do j = 1, n_degen
   write(*,'(100(F16.10,X))')stmp(:,j)
  enddo
  call build_s_matrix(ao_num,n_degen,mo_l_coef_new,mo_l_coef_new,ao_overlap,stmp)
  print*,'LEFT OVERLAP '
  do j = 1, n_degen
   write(*,'(100(F16.10,X))')stmp(:,j)
  enddo
  call build_s_matrix(ao_num,n_degen,mo_r_coef_tmp,mo_r_coef_tmp,ao_overlap,stmp)
  print*,'RIGHT OVERLAP '
  do j = 1, n_degen
   write(*,'(100(F16.10,X))')stmp(:,j)
  enddo
  do j = 1, n_degen
   mo_l_coef_good(1:ao_num,j+ifirst-1) = mo_l_coef_new(1:ao_num,j)
   mo_r_coef_good(1:ao_num,j+ifirst-1) = mo_r_coef_tmp(1:ao_num,j)
  enddo
  deallocate(stmp,smat2)
  deallocate(mo_r_coef_tmp,mo_l_coef_tmp,mo_l_coef_new)
  deallocate(T,Snew)
 enddo

 allocate(stmp(mo_num, mo_num))
 print*,'l coef'
 do i = 1, mo_num
  write(*,'(100(F8.4,X))')mo_l_coef_good(:,i)
 enddo
 print*,'r coef'
 do i = 1, mo_num
  write(*,'(100(F8.4,X))')mo_r_coef_good(:,i)
 enddo
 call build_s_matrix(ao_num,mo_num,mo_l_coef_good,mo_r_coef_good,ao_overlap,stmp)
  print*,'LEFT/RIGHT OVERLAP '
  do j = 1, mo_num
   write(*,'(100(F16.10,X))')stmp(:,j)
  enddo
 call build_s_matrix(ao_num,mo_num,mo_l_coef_good,mo_l_coef_good,ao_overlap,stmp)
  print*,'LEFT/LEFT OVERLAP '
  do j = 1, mo_num
   write(*,'(100(F16.10,X))')stmp(:,j)
  enddo
 call build_s_matrix(ao_num,mo_num,mo_r_coef_good,mo_r_coef_good,ao_overlap,stmp)
  print*,'RIGHT/RIGHT OVERLAP '
  do j = 1, mo_num
   write(*,'(100(F16.10,X))')stmp(:,j)
  enddo
  mo_r_coef = mo_r_coef_good
  mo_l_coef = mo_l_coef_good
  call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
  call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
  TOUCH mo_l_coef mo_r_coef
end

subroutine build_s_matrix(m,n,C1,C2,overlap,smat)
 implicit none
 integer, intent(in) :: m,n
 double precision, intent(in) :: C1(m,n),C2(m,n),overlap(m,m)
 double precision, intent(out):: smat(n,n)
 integer :: i,j,k,l
 smat = 0.D0
  do i = 1, n
   do j = 1, n
    do k = 1, m
     do l = 1, m
      smat(i,j) += C1(k,i) * overlap(l,k) * C2(l,j) 
     enddo
    enddo
   enddo
  enddo
end

subroutine orthog_functions(m,n,coef,overlap)
 implicit none
 integer, intent(in) :: m,n
 double precision, intent(in)    :: overlap(m,m)
 double precision, intent(inout) :: coef(m,n)
 double precision, allocatable :: stmp(:,:)
 integer :: j
 allocate(stmp(n,n))
  call build_s_matrix(m,n,coef,coef,overlap,stmp)
  print*,'overlap before'
  do j = 1, n
   write(*,'(100(F16.10,X))')stmp(:,j)
  enddo
  call impose_orthog_svd_overlap(m, n, coef,overlap)
  call build_s_matrix(m,n,coef,coef,overlap,stmp)
  do j = 1, n
   coef(1,:m) *= 1.d0/dsqrt(stmp(j,j))
  enddo
  print*,'overlap after'
  call build_s_matrix(m,n,coef,coef,overlap,stmp)
  do j = 1, n
   write(*,'(100(F16.10,X))')stmp(:,j)
  enddo
end
