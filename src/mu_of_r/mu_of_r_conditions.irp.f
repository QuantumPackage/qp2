
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
! PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals 
 call wall_time(wall0)

 if (read_mu_of_r) then
   print*,'Reading mu(r) from disk ...'
   call ezfio_get_mu_of_r_mu_of_r_disk(mu_of_r_prov)
   return
 endif

 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   if(mu_of_r_potential.EQ."hf")then
    mu_of_r_prov(ipoint,istate) =  mu_of_r_hf(ipoint)
   else if(mu_of_r_potential.EQ."cas_ful".or.mu_of_r_potential.EQ."cas_truncated")then
    mu_of_r_prov(ipoint,istate) =  mu_of_r_psi_cas(ipoint,istate)
   else 
    print*,'you requested the following mu_of_r_potential'
    print*,mu_of_r_potential
    print*,'which does not correspond to any of the options for such keyword'
    stop
   endif
  enddo
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
 PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals 
 print*,'providing mu_of_r_hf ...'
 call wall_time(wall0)
 sqpi = dsqrt(dacos(-1.d0))
 provide f_psi_hf_ab 
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,f_hf,on_top,w_hf) & 
 !$OMP ShARED (n_points_final_grid,mu_of_r_hf,f_psi_hf_ab,on_top_hf_mu_r,sqpi) 
 do ipoint = 1, n_points_final_grid
  f_hf   = f_psi_hf_ab(ipoint)
  on_top = on_top_hf_mu_r(ipoint)
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
 ! average value of mu(r) weighted with the total one-e density and divised by the number of electrons 
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

