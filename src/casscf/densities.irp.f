use bitmasks

BEGIN_PROVIDER [real*8, D0tu, (n_act_orb,n_act_orb) ]
  implicit none
  BEGIN_DOC
  ! the first-order density matrix in the basis of the starting MOs.
  ! matrix is state averaged.
  END_DOC
  integer                        :: t,u
  
  do u=1,n_act_orb
    do t=1,n_act_orb
      D0tu(t,u) = one_e_dm_mo_alpha_average( list_act(t), list_act(u) ) + &
                  one_e_dm_mo_beta_average ( list_act(t), list_act(u) ) 
    enddo
  enddo
  
END_PROVIDER

BEGIN_PROVIDER [real*8, P0tuvx, (n_act_orb,n_act_orb,n_act_orb,n_act_orb) ]
   BEGIN_DOC
   ! The second-order density matrix in the basis of the starting MOs ONLY IN THE RANGE OF ACTIVE MOS
   ! The values are state averaged
   !
   ! We use the spin-free generators of mono-excitations
   ! E_pq destroys q and creates p
   ! D_pq   =     <0|E_pq|0> = D_qp
   ! P_pqrs = 1/2 <0|E_pq E_rs - delta_qr E_ps|0>
   !
   ! P0tuvx(p,q,r,s) = chemist notation : 1/2 <0|E_pq E_rs - delta_qr E_ps|0>
   END_DOC
   implicit none
   integer                        :: t,u,v,x
   integer                        :: tt,uu,vv,xx
   integer                        :: mu,nu,istate,ispin,jspin,ihole,ipart,jhole,jpart
   integer                        :: ierr
   real*8                         :: phase1,phase11,phase12,phase2,phase21,phase22
   integer                        :: nu1,nu2,nu11,nu12,nu21,nu22
   integer                        :: ierr1,ierr2,ierr11,ierr12,ierr21,ierr22
   real*8                         :: cI_mu(N_states),term
   integer(bit_kind), dimension(N_int,2) :: det_mu, det_mu_ex
   integer(bit_kind), dimension(N_int,2) :: det_mu_ex1, det_mu_ex11, det_mu_ex12
   integer(bit_kind), dimension(N_int,2) :: det_mu_ex2, det_mu_ex21, det_mu_ex22
   
  if (bavard) then
    write(6,*) ' providing the 2 body RDM on the active part'
  endif

  P0tuvx= 0.d0
  do istate=1,N_states
   do x = 1, n_act_orb
    xx = list_act(x)
    do v = 1, n_act_orb
     vv = list_act(v)
     do u = 1, n_act_orb
      uu = list_act(u)
      do t = 1, n_act_orb
       tt = list_act(t)
       P0tuvx(t,u,v,x) = state_av_act_2_rdm_spin_trace_mo(t,v,u,x)
      enddo
     enddo 
    enddo
   enddo
  enddo

END_PROVIDER
