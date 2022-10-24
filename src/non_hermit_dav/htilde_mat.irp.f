BEGIN_PROVIDER [ integer, n_mat]
 implicit none
 n_mat = 2
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, h_non_hermit, (n_mat, n_mat)]
&BEGIN_PROVIDER [ double precision, h_non_hermit_transp, (n_mat, n_mat)]
&BEGIN_PROVIDER [ double precision, reigvec_ht, (n_mat, n_mat)]
&BEGIN_PROVIDER [ double precision, leigvec_ht, (n_mat, n_mat)]
&BEGIN_PROVIDER [ double precision, eigval_ht, (n_mat)]
&BEGIN_PROVIDER [ integer, n_real_ht, (n_mat)]
 implicit none
 integer :: i,j
 do i = 1, n_mat
  read(33,*)h_non_hermit(i,1:n_mat)
 enddo
 print*,''
 print*,'H_mat '
 print*,''
 do i = 1, n_mat
  write(*,'(1000(F16.10,X))')h_non_hermit(i,:)
 enddo
 do i = 1, n_mat
  do j = 1, n_mat
   h_non_hermit_transp(j,i) = h_non_hermit(i,j)
  enddo
 enddo
 call non_hrmt_real_diag(n_mat,h_non_hermit,reigvec_ht,leigvec_ht,n_real_ht,eigval_ht)


END_PROVIDER 


subroutine hcalc_r_tmp(v,u,N_st,sze) ! v = H u
  implicit none
  BEGIN_DOC
  ! Template of routine for the application of H
  !
  ! Here, it is done with the Hamiltonian matrix 
  !
  ! on the set of determinants of psi_det 
  !
  ! Computes $v = H | u \rangle$ 
  !
  END_DOC
  integer, intent(in)              :: N_st,sze
  double precision, intent(in)     :: u(sze,N_st)
  double precision, intent(inout)  :: v(sze,N_st)
  integer :: i,j,istate
  v = 0.d0
  do istate = 1, N_st
   do j = 1, sze
    do i = 1, sze
      v(i,istate) += h_non_hermit(i,j) * u(j,istate)
!      print*,i,j,h_non_hermit(i,j),u(j,istate)
    enddo
   enddo
  enddo
  print*,'HU'
  do i = 1, sze
   print*,v(i,1)
  enddo
end

subroutine hcalc_l_tmp(v,u,N_st,sze) ! v = H^\dagger u
  implicit none
  BEGIN_DOC
  ! Template of routine for the application of H
  !
  ! Here, it is done with the Hamiltonian matrix 
  !
  ! on the set of determinants of psi_det 
  !
  ! Computes $v = H | u \rangle$ 
  !
  END_DOC
  integer, intent(in)              :: N_st,sze
  double precision, intent(in)     :: u(sze,N_st)
  double precision, intent(inout)  :: v(sze,N_st)
  integer :: i,j,istate
  v = 0.d0
  do istate = 1, N_st
   do j = 1, sze
    do i = 1, sze
      v(i,istate) += h_non_hermit_transp(i,j) * u(j,istate)
    enddo
   enddo
  enddo
  print*,'HU'
  do i = 1, sze
   print*,v(i,1)
  enddo
end
