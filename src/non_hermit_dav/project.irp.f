subroutine h_non_hermite(v,u,Hmat,a,N_st,sze) 
 implicit none
 BEGIN_DOC
 ! Template of routine for the application of H
 !
 ! Here, it is done with the Hamiltonian matrix 
 !
 ! on the set of determinants of psi_det 
 !
 ! Computes $v = a * H | u \rangle$ 
 !
 END_DOC
 integer, intent(in)              :: N_st,sze
 double precision, intent(in)     :: u(sze,N_st), Hmat(sze,sze), a
 double precision, intent(inout)  :: v(sze,N_st)
 integer :: i,j,k 
 do k = 1, N_st
  do j = 1, sze
   do i = 1, sze
    v(i,k) += a * u(j,k) * Hmat(i,j)
   enddo
  enddo
 enddo
end


subroutine exp_tau_H(u,v,hmat,tau,et,N_st,sze)
 implicit none
 BEGIN_DOC
! realises v = (1 - tau (H - et)) u
 END_DOC
 integer, intent(in) :: N_st,sze
 double precision, intent(in) :: hmat(sze,sze), u(sze,N_st), tau, et
 double precision, intent(out):: v(sze,N_st)
 double precision :: a
 integer :: i,j
 v = (1.d0 + tau * et) * u 
 a = -1.d0 * tau
 call h_non_hermite(v,u,Hmat,a,N_st,sze)
end

double precision function project_phi0(u,Hmat0,N_st,sze)
 implicit none
 integer, intent(in)              :: N_st,sze
 double precision, intent(in)     :: u(sze,N_st), Hmat0(sze)
 integer :: j
 project_phi0 = 0.d0
 do j = 1, sze
  project_phi0 += u(j,1) * Hmat0(j) 
 enddo
 project_phi0 *= 1.d0 / u(1,1)
end

