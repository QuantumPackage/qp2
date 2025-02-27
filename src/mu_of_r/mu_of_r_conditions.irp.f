
 BEGIN_PROVIDER [double precision, mu_of_r_prov, (n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
 ! general variable for mu(r)
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018)
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint,istate
 double precision :: wall0,wall1
 print*,'providing mu_of_r ...'
 call wall_time(wall0)

 if (read_mu_of_r) then
   print*,'Reading mu(r) from disk ...'
   call ezfio_get_mu_of_r_mu_of_r_disk(mu_of_r_prov)
   return
 endif

 do istate = 1, N_states
   if(mu_of_r_potential.EQ."hf")then
     do ipoint = 1, n_points_final_grid
       mu_of_r_prov(ipoint,istate) =  mu_of_r_hf(ipoint)
     enddo
   else if(mu_of_r_potential.EQ."hf_old")then
     do ipoint = 1, n_points_final_grid
       mu_of_r_prov(ipoint,istate) =  mu_of_r_hf_old(ipoint)
     enddo
   else if(mu_of_r_potential.EQ."hf_sparse")then
     do ipoint = 1, n_points_final_grid
       mu_of_r_prov(ipoint,istate) =  mu_of_r_hf_sparse(ipoint)
     enddo
   else if(mu_of_r_potential.EQ."cas_full".or.mu_of_r_potential.EQ."cas_truncated".or.mu_of_r_potential.EQ."pure_act")then
     do ipoint = 1, n_points_final_grid
       mu_of_r_prov(ipoint,istate) =  mu_of_r_psi_cas(ipoint,istate)
     enddo
   else if(mu_of_r_potential.EQ."proj")then
     do ipoint = 1, n_points_final_grid
       mu_of_r_prov(ipoint,istate) =  mu_of_r_projector_mo(ipoint)
     enddo
   else
    print*,'you requested the following mu_of_r_potential'
    print*,mu_of_r_potential
    print*,'which does not correspond to any of the options for such keyword'
    stop
   endif
 enddo

 if (write_mu_of_r) then
   print*,'Writing mu(r) on disk ...'
   call ezfio_set_mu_of_r_io_mu_of_r('Read')
   call ezfio_set_mu_of_r_mu_of_r_disk(mu_of_r_prov)
 endif

 call wall_time(wall1)
 print*,'Time to provide mu_of_r = ',wall1-wall0
 END_PROVIDER


 BEGIN_PROVIDER [double precision, mu_of_r_hf, (n_points_final_grid) ]
 implicit none
 BEGIN_DOC
 ! mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint
 double precision :: wall0,wall1,f_hf,on_top,w_hf,sqpi
 print*,'providing mu_of_r_hf ...'
 call wall_time(wall0)
 PROVIDE f_hf_cholesky on_top_hf_grid
 sqpi = dsqrt(dacos(-1.d0))
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,f_hf,on_top,w_hf) &
 !$OMP ShARED (n_points_final_grid,mu_of_r_hf,f_hf_cholesky,on_top_hf_grid,sqpi)
 do ipoint = 1, n_points_final_grid
  f_hf   = f_hf_cholesky(ipoint)
  on_top = on_top_hf_grid(ipoint)
  if(on_top.le.1.d-12.or.f_hf.le.0.d0.or.f_hf * on_top.lt.0.d0)then
    w_hf   = 1.d+10
  else
    w_hf  = f_hf /  on_top
  endif
  mu_of_r_hf(ipoint) =  w_hf * sqpi * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide mu_of_r_hf = ',wall1-wall0
 END_PROVIDER

 BEGIN_PROVIDER [double precision, mu_of_r_hf_sparse, (n_points_final_grid) ]
 implicit none
 BEGIN_DOC
 ! mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint
 double precision :: wall0,wall1,f_hf,on_top,w_hf,sqpi
 print*,'providing mu_of_r_hf_sparse ...'
 call wall_time(wall0)
 sqpi = dsqrt(dacos(-1.d0))
 PROVIDE f_hf_cholesky_sparse on_top_hf_grid
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,f_hf,on_top,w_hf) &
 !$OMP ShARED (n_points_final_grid,mu_of_r_hf_sparse,f_hf_cholesky_sparse,on_top_hf_grid,sqpi)
 do ipoint = 1, n_points_final_grid
  f_hf   = f_hf_cholesky_sparse(ipoint)
  on_top = on_top_hf_grid(ipoint)
  if(on_top.le.1.d-12.or.f_hf.le.0.d0.or.f_hf * on_top.lt.0.d0)then
    w_hf   = 1.d+10
  else
    w_hf  = f_hf /  on_top
  endif
  mu_of_r_hf_sparse(ipoint) =  w_hf * sqpi * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide mu_of_r_hf_sparse = ',wall1-wall0
 END_PROVIDER

 BEGIN_PROVIDER [double precision, mu_of_r_hf_old, (n_points_final_grid) ]
 implicit none
 BEGIN_DOC
 ! mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint
 double precision :: wall0,wall1,f_hf,on_top,w_hf,sqpi
 PROVIDE all_mo_integrals
 print*,'providing mu_of_r_hf_old ...'
 call wall_time(wall0)
 sqpi = dsqrt(dacos(-1.d0))
 provide f_psi_hf_ab
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,f_hf,on_top,w_hf) &
 !$OMP ShARED (n_points_final_grid,mu_of_r_hf_old,f_psi_hf_ab,on_top_hf_mu_r,sqpi)
 do ipoint = 1, n_points_final_grid
  f_hf   = f_psi_hf_ab(ipoint)
  on_top = on_top_hf_mu_r(ipoint)
  if(on_top.le.1.d-12.or.f_hf.le.0.d0.or.f_hf * on_top.lt.0.d0)then
    w_hf   = 1.d+10
  else
    w_hf  = f_hf /  on_top
  endif
  mu_of_r_hf_old(ipoint) =  w_hf * sqpi * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide mu_of_r_hf_old = ',wall1-wall0
 END_PROVIDER


 BEGIN_PROVIDER [double precision, mu_of_r_psi_cas, (n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
 ! mu(r) computed with a wave function developped in an active space
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018)
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
 !
 ! in the one- and two-body density matrix are excluded
 END_DOC
 integer :: ipoint,istate
 double precision :: wall0,wall1,f_psi,on_top,w_psi,sqpi
 print*,'providing mu_of_r_psi_cas ...'
 call wall_time(wall0)
 sqpi = dsqrt(dacos(-1.d0))

 provide f_psi_cas_ab
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,f_psi,on_top,w_psi,istate) &
 !$OMP SHARED (n_points_final_grid,mu_of_r_psi_cas,f_psi_cas_ab,on_top_cas_mu_r,sqpi,N_states)
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   f_psi  = f_psi_cas_ab(ipoint,istate)
   on_top = on_top_cas_mu_r(ipoint,istate)
   if(on_top.le.1.d-12.or.f_psi.le.0.d0.or.f_psi * on_top.lt.0.d0)then
     w_psi   = 1.d+10
   else
     w_psi  = f_psi /  on_top
   endif
   mu_of_r_psi_cas(ipoint,istate) = w_psi * sqpi * 0.5d0
  enddo
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide mu_of_r_psi_cas = ',wall1-wall0
 END_PROVIDER


 BEGIN_PROVIDER [double precision, mu_average_prov, (N_states)]
 implicit none
 BEGIN_DOC
 ! average value of mu(r) weighted with the total one-e density and divided by the number of electrons
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
 !
 ! in the one- and two-body density matrix are excluded
 END_DOC
 integer :: ipoint,istate
 double precision :: weight,density
 mu_average_prov = 0.d0
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   weight =final_weight_at_r_vector(ipoint)
   density = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate) &
           + one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   if(mu_of_r_prov(ipoint,istate).gt.1.d+09)cycle
   mu_average_prov(istate) += mu_of_r_prov(ipoint,istate) * weight * density
  enddo
  mu_average_prov(istate) = mu_average_prov(istate) / elec_num_grid_becke(istate)
 enddo
 END_PROVIDER


BEGIN_PROVIDER [double precision, mu_of_r_projector_mo, (n_points_final_grid) ]
 implicit none
 BEGIN_DOC
 ! mu(r) computed with the projector onto the atomic basis
 !  P_B(\mathbf{r},\mathbf{r}') = \sum_{ij} |
 !  \chi_{i} \rangle \left[S^{-1}\right]_{ij} \langle \chi_{j} |
 !  \] where $i$ and $j$ denote all atomic orbitals.
 END_DOC

 double precision, parameter :: factor = dsqrt(2.d0*dacos(-1.d0))
 double precision, allocatable :: tmp(:,:)
 integer :: ipoint


 do ipoint=1,n_points_final_grid
   mu_of_r_projector_mo(ipoint) = 0.d0
   integer :: i,j
   do j=1,n_inact_act_orb
     i = list_inact_act(j)
     mu_of_r_projector_mo(ipoint) = mu_of_r_projector_mo(ipoint) + &
         mos_in_r_array_omp(i,ipoint) * mos_in_r_array_omp(i,ipoint)
   enddo
   do j=1,n_virt_orb
     i = list_virt(j)
     mu_of_r_projector_mo(ipoint) = mu_of_r_projector_mo(ipoint) + &
         mos_in_r_array_omp(i,ipoint) * mos_in_r_array_omp(i,ipoint)
   enddo
 enddo

 do ipoint=1,n_points_final_grid
   ! epsilon
   mu_of_r_projector_mo(ipoint) = 1.d0/(2.d0*dacos(-1.d0) * mu_of_r_projector_mo(ipoint)**(2.d0/3.d0))
   ! mu
   mu_of_r_projector_mo(ipoint) = 1.d0/dsqrt( 2.d0*mu_of_r_projector_mo(ipoint) )
 enddo
END_PROVIDER



BEGIN_PROVIDER [double precision, mu_average_proj, (N_states)]
 implicit none
 BEGIN_DOC
 ! average value of mu(r) weighted with the total one-e density and divided by the number of electrons
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals
 !
 ! in the one- and two-body density matrix are excluded
 END_DOC
 integer :: ipoint,istate
 double precision :: weight,density
 do istate = 1, N_states
  mu_average_proj(istate) = 0.d0
  do ipoint = 1, n_points_final_grid
   weight =final_weight_at_r_vector(ipoint)
   density = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate) &
           + one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   mu_average_proj(istate) += mu_of_r_projector_mo(ipoint) * weight * density
  enddo
  mu_average_proj(istate) = mu_average_proj(istate) / elec_num_grid_becke(istate)
 enddo
END_PROVIDER

