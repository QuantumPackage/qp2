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

 BEGIN_PROVIDER [double precision, D0tu_alpha_ao, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, D0tu_beta_ao, (ao_num, ao_num)]
 implicit none
 integer :: i,ii,j,u,t,uu,tt
 double precision, allocatable :: D0_tmp_alpha(:,:),D0_tmp_beta(:,:)
 allocate(D0_tmp_alpha(mo_num, mo_num),D0_tmp_beta(mo_num, mo_num))
 D0_tmp_beta = 0.d0
 D0_tmp_alpha = 0.d0
 do i = 1, n_core_inact_orb
  ii = list_core_inact(i)
  D0_tmp_alpha(ii,ii) = 1.d0
  D0_tmp_beta(ii,ii) = 1.d0
 enddo
 print*,'Diagonal elements of the 1RDM in the active space'
 do u=1,n_act_orb
   uu = list_act(u)
   print*,uu,one_e_dm_mo_alpha_average(uu,uu),one_e_dm_mo_beta_average(uu,uu)
   do t=1,n_act_orb
    tt = list_act(t)
    D0_tmp_alpha(tt,uu) = one_e_dm_mo_alpha_average(tt,uu)
    D0_tmp_beta(tt,uu) = one_e_dm_mo_beta_average(tt,uu)
   enddo
 enddo

 call mo_to_ao_no_overlap(D0_tmp_alpha,mo_num,D0tu_alpha_ao,ao_num)
 call mo_to_ao_no_overlap(D0_tmp_beta,mo_num,D0tu_beta_ao,ao_num)

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
  if(fast_2rdm)then
   do istate=1,N_states
    do x = 1, n_act_orb
     do v = 1, n_act_orb
      do u = 1, n_act_orb
       do t = 1, n_act_orb
        !      1 1 2 2                                     1 2 1 2
        P0tuvx(t,u,v,x) = 0.5d0 * state_av_act_2_rdm_spin_trace_mo(t,v,u,x)
       enddo
      enddo 
     enddo
    enddo
   enddo
  else
   P0tuvx = P0tuvx_peter
  endif

END_PROVIDER
