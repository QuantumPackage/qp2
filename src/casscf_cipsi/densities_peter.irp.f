use bitmasks

BEGIN_PROVIDER [real*8, P0tuvx_peter, (n_act_orb,n_act_orb,n_act_orb,n_act_orb) ]
   BEGIN_DOC
   ! the second-order density matrix in the basis of the starting MOs
   ! matrices are state averaged
   !
   ! we use the spin-free generators of mono-excitations
   ! E_pq destroys q and creates p
   ! D_pq   =     <0|E_pq|0> = D_qp
   ! P_pqrs = 1/2 <0|E_pq E_rs - delta_qr E_ps|0>
   !
   END_DOC
   implicit none
   integer                        :: t,u,v,x,mu,nu,istate,ispin,jspin,ihole,ipart,jhole,jpart
   integer                        :: ierr
   real*8                         :: phase1,phase11,phase12,phase2,phase21,phase22
   integer                        :: nu1,nu2,nu11,nu12,nu21,nu22
   integer                        :: ierr1,ierr2,ierr11,ierr12,ierr21,ierr22
   real*8                         :: cI_mu(N_states),term
   integer(bit_kind), dimension(N_int,2) :: det_mu, det_mu_ex
   integer(bit_kind), dimension(N_int,2) :: det_mu_ex1, det_mu_ex11, det_mu_ex12
   integer(bit_kind), dimension(N_int,2) :: det_mu_ex2, det_mu_ex21, det_mu_ex22
   
  if (bavard) then
    write(6,*) ' providing density matrix P0'
  endif

   P0tuvx_peter = 0.d0
   
   ! first loop: we apply E_tu, once for D_tu, once for -P_tvvu
   do mu=1,n_det
     call det_extract(det_mu,mu,N_int)
     do istate=1,n_states
       cI_mu(istate)=psi_coef(mu,istate)
     end do
     do t=1,n_act_orb
       ipart=list_act(t)
       do u=1,n_act_orb
         ihole=list_act(u)
         ! apply E_tu
         call det_copy(det_mu,det_mu_ex1,N_int)
         call det_copy(det_mu,det_mu_ex2,N_int)
         call do_spinfree_mono_excitation(det_mu,det_mu_ex1          &
             ,det_mu_ex2,nu1,nu2,ihole,ipart,phase1,phase2,ierr1,ierr2)
         ! det_mu_ex1 is in the list
         if (nu1.ne.-1) then
           do istate=1,n_states
             term=cI_mu(istate)*psi_coef(nu1,istate)*phase1
             ! and we fill P0_tvvu
             do v=1,n_act_orb
               P0tuvx_peter(t,v,v,u)-=term
             end do
           end do
         end if
         ! det_mu_ex2 is in the list
         if (nu2.ne.-1) then
           do istate=1,n_states
             term=cI_mu(istate)*psi_coef(nu2,istate)*phase2
             do v=1,n_act_orb
               P0tuvx_peter(t,v,v,u)-=term
             end do
           end do
         end if
       end do
     end do
   end do
   ! now we do the double excitation E_tu E_vx |0>
   do mu=1,n_det
     call det_extract(det_mu,mu,N_int)
     do istate=1,n_states
       cI_mu(istate)=psi_coef(mu,istate)
     end do
     do v=1,n_act_orb
       ipart=list_act(v)
       do x=1,n_act_orb
         ihole=list_act(x)
         ! apply E_vx
         call det_copy(det_mu,det_mu_ex1,N_int)
         call det_copy(det_mu,det_mu_ex2,N_int)
         call do_spinfree_mono_excitation(det_mu,det_mu_ex1          &
             ,det_mu_ex2,nu1,nu2,ihole,ipart,phase1,phase2,ierr1,ierr2)
         ! we apply E_tu to the first resultant determinant, thus E_tu E_vx |0>
         if (ierr1.eq.1) then
           do t=1,n_act_orb
             jpart=list_act(t)
             do u=1,n_act_orb
               jhole=list_act(u)
               call det_copy(det_mu_ex1,det_mu_ex11,N_int)
               call det_copy(det_mu_ex1,det_mu_ex12,N_int)
               call do_spinfree_mono_excitation(det_mu_ex1,det_mu_ex11&
                   ,det_mu_ex12,nu11,nu12,jhole,jpart,phase11,phase12,ierr11,ierr12)
               if (nu11.ne.-1) then
                 do istate=1,n_states
                   P0tuvx_peter(t,u,v,x)+=cI_mu(istate)*psi_coef(nu11,istate)&
                       *phase11*phase1
                 end do
               end if
               if (nu12.ne.-1) then
                 do istate=1,n_states
                   P0tuvx_peter(t,u,v,x)+=cI_mu(istate)*psi_coef(nu12,istate)&
                       *phase12*phase1
                 end do
               end if
             end do
           end do
         end if
         
         ! we apply E_tu to the second resultant determinant
         if (ierr2.eq.1) then
           do t=1,n_act_orb
             jpart=list_act(t)
             do u=1,n_act_orb
               jhole=list_act(u)
               call det_copy(det_mu_ex2,det_mu_ex21,N_int)
               call det_copy(det_mu_ex2,det_mu_ex22,N_int)
               call do_spinfree_mono_excitation(det_mu_ex2,det_mu_ex21&
                   ,det_mu_ex22,nu21,nu22,jhole,jpart,phase21,phase22,ierr21,ierr22)
               if (nu21.ne.-1) then
                 do istate=1,n_states
                   P0tuvx_peter(t,u,v,x)+=cI_mu(istate)*psi_coef(nu21,istate)&
                       *phase21*phase2
                 end do
               end if
               if (nu22.ne.-1) then
                 do istate=1,n_states
                   P0tuvx_peter(t,u,v,x)+=cI_mu(istate)*psi_coef(nu22,istate)&
                       *phase22*phase2
                 end do
               end if
             end do
           end do
         end if
         
       end do
     end do
   end do
   
   ! we average by just dividing by the number of states
   do x=1,n_act_orb
     do v=1,n_act_orb
       do u=1,n_act_orb
         do t=1,n_act_orb
           P0tuvx_peter(t,u,v,x)*=0.5D0/dble(N_states)
         end do
       end do
     end do
   end do
   
END_PROVIDER
