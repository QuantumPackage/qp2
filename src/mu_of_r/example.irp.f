
subroutine test_f_HF_valence_ab
 implicit none
 BEGIN_DOC
! routine to test the function f_HF(r1,r2) 
!
! the integral over r1,r2 should be equal to the alpha/beta interaction of HF determinant
 END_DOC
 integer :: ipoint,i,j,i_i,j_j,jpoint
 double precision :: accu_val,accu_ful, weight1,weight2, r1(3),integral_psi_val,integral_psi,r2(3),two_bod
 accu_2 = 0.d0
 ! You compute the coulomb repulsion between alpha-beta electrons for HF
 do i = 1, n_occ_val_orb_for_hf(1)
  i_i = list_valence_orb_for_hf(i,1)
  do j = 1, n_occ_val_orb_for_hf(2)
   j_j = list_valence_orb_for_hf(j,2)
   accu_2 += mo_two_e_integrals_jj(j_j,i_i)
  enddo
 enddo
 print*,''
 print*,''
 print*,''
 print*,'**************************'
 print*,'**************************'
 print*,'Routine to test the f_HF(r1,r2) function'
 print*,'**************************'
 print*,''
 print*,''
 print*,''
 print*,'**************************'
 print*,'<HF| We_ee^{ab}|HF>     = ',accu_2
 print*,'**************************'

 print*,'semi analytical form '
 accu_val = 0.d0
 ! You integrate on r2 the analytical integral over r1 of f_HF(r1,r2) 
 do ipoint  = 1, n_points_final_grid
  weight1 =final_weight_at_r_vector(ipoint)
  r2(1)   = final_grid_points(1,ipoint)
  r2(2)   = final_grid_points(2,ipoint)
  r2(3)   = final_grid_points(3,ipoint)
  call integral_f_HF_valence_ab(r2,integral_psi_val)
  accu_val += integral_psi_val * weight1
 enddo
 print*,'**************************'
 ! Should give you the alpha-beta repulsion of HF, excluding core contributions,  
 print*,'int dr1 dr2 f_HF(r1,r2) = ',accu_val
 double precision :: accu_2


 print*,'pure numerical form (might take quite some time as it grows as N_g^2 * N_e^2 * N_b^2 ...)'
 ! You integrate brut force on r1 and r2 
 accu_val = 0.d0
 do jpoint = 1, n_points_final_grid
  weight1 =final_weight_at_r_vector(jpoint)
  r1(1)   = final_grid_points(1,jpoint)
  r1(2)   = final_grid_points(2,jpoint)
  r1(3)   = final_grid_points(3,jpoint)
  do ipoint  = 1, n_points_final_grid
   weight2 =final_weight_at_r_vector(ipoint)
   r2(1)   = final_grid_points(1,ipoint)
   r2(2)   = final_grid_points(2,ipoint)
   r2(3)   = final_grid_points(3,ipoint)
   call f_HF_valence_ab(r1,r2,integral_psi_val,two_bod)
   accu_val += integral_psi_val * weight1 * weight2
  enddo
 enddo
 print*,'int dr1 dr2 f_HF(r1,r2) = ',accu_val


 print*,'**************************'
 print*,'**************************'
 print*,'**************************'
 accu_val = 0.d0
 r1 = 0.d0
 r1(1) = 0.5d0
 print*,'r1 = ',r1
 ! You compute the integral over r2 of f_HF(r1,r2) 
 call integral_f_HF_valence_ab(r1,integral_psi)
 do ipoint  = 1, n_points_final_grid
  weight1 =final_weight_at_r_vector(ipoint)
  r2(1)   = final_grid_points(1,ipoint)
  r2(2)   = final_grid_points(2,ipoint)
  r2(3)   = final_grid_points(3,ipoint)
  call f_HF_valence_ab(r1,r2,integral_psi_val,two_bod)
  accu_val += integral_psi_val * weight1
 enddo
 print*,'int dr2 f_HF(r1,r2)     = ',integral_psi
 print*,'analytical form         = ',accu_val
 print*,'**************************'
end



subroutine test_f_ii_valence_ab
 implicit none
 BEGIN_DOC
! routine to test the function f_ii(r1,r2) 
!
! it should be the same that f_HF(r1,r2) only for inactive orbitals 
 END_DOC
 integer :: ipoint
 double precision :: accu_f, accu_n2, weight, r1(3),r2(3)
 double precision :: accu_f_on_top
 double precision :: f_HF_val_ab,two_bod_dens_hf,f_ii_val_ab,two_bod_dens_ii
 accu_f  = 0.d0
 accu_n2 = 0.d0
 do ipoint  = 1, n_points_final_grid
  weight  = final_weight_at_r_vector(ipoint)
  r1(1)   = final_grid_points(1,ipoint)
  r1(2)   = final_grid_points(2,ipoint)
  r1(3)   = final_grid_points(3,ipoint)
  r2 = r1
  call f_HF_valence_ab(r1,r2,f_HF_val_ab,two_bod_dens_hf)
  call give_f_ii_val_ab(r1,r2,f_ii_val_ab,two_bod_dens_ii)
  accu_f += dabs(f_HF_val_ab - f_ii_val_ab) * weight
  accu_n2+= dabs(two_bod_dens_hf - two_bod_dens_ii) * weight
  accu_f_on_top += dabs(two_bod_dens_hf) * weight
 enddo
 print*,'**************************'
 print*,''
 print*,'accu_f  = ',accu_f
 print*,'accu_n2 = ',accu_n2
 print*,''
 print*,'accu_f_on_top  = ',accu_f_on_top
end


subroutine test_f_ia_valence_ab
 implicit none
 BEGIN_DOC
! routine to test the function f_ii(r1,r2), f_ia(r1,r2) and f_aa(r1,r2)
 END_DOC
 integer :: ipoint,istate
 double precision :: accu_f, accu_n2, weight, r1(3),r2(3)
 double precision :: accu_f_on_top
 double precision :: f_ref,f_comp,on_top_ref,on_top_comp
 double precision :: f_ii_val_ab,two_bod_dens_ii,f_ia_val_ab,two_bod_dens_ia,f_aa_val_ab,two_bod_dens_aa
 double precision :: accu
 accu_f  = 0.d0
 accu_n2 = 0.d0
 accu = 0.d0
 istate = 1
 do ipoint  = 1, n_points_final_grid
  weight  = final_weight_at_r_vector(ipoint)
  r1(1)   = final_grid_points(1,ipoint)
  r1(2)   = final_grid_points(2,ipoint)
  r1(3)   = final_grid_points(3,ipoint)
  r2 = r1
  call give_f_ii_val_ab(r1,r2,f_ii_val_ab,two_bod_dens_ii)
  call give_f_ia_val_ab(r1,r2,f_ia_val_ab,two_bod_dens_ia,istate)
  call give_f_aa_val_ab(r1,r2,f_aa_val_ab,two_bod_dens_aa,istate)
  f_ref = f_psi_cas_ab_old(ipoint,istate)
  f_comp = f_ii_val_ab + f_ia_val_ab + f_aa_val_ab
  on_top_ref = total_cas_on_top_density(ipoint,istate)
  on_top_comp= two_bod_dens_ii + two_bod_dens_ia + two_bod_dens_aa 
  accu_f += dabs(f_ref - f_comp) * weight
  accu_n2+= dabs(on_top_ref - on_top_comp) * weight
  accu += f_ref * weight
 enddo
 print*,'**************************'
 print*,''
 print*,'accu_f  = ',accu_f
 print*,'accu_n2 = ',accu_n2
 print*,''
 print*,'accu    = ',accu

end

subroutine test_f_ii_ia_aa_valence_ab
 implicit none
 BEGIN_DOC
! routine to test the function f_Psi(r1,r2) based on core/inactive/active orbitals 
 END_DOC
 integer :: ipoint,istate
 double precision :: accu_f, accu_n2, weight, r1(3),r2(3)
 double precision :: accu_f_on_top
 double precision :: f_ref,f_comp,on_top_ref,on_top_comp
 double precision :: f_ii_val_ab,two_bod_dens_ii,f_ia_val_ab,two_bod_dens_ia,f_aa_val_ab,two_bod_dens_aa
 double precision :: accu
 accu_f  = 0.d0
 accu_n2 = 0.d0
 accu = 0.d0
 istate = 1
 do ipoint  = 1, n_points_final_grid
  weight  = final_weight_at_r_vector(ipoint)
  f_ref = f_psi_cas_ab(ipoint,istate)
  f_comp = f_psi_cas_ab_old(ipoint,istate)
  on_top_ref = total_cas_on_top_density(ipoint,istate)
  on_top_comp= on_top_cas_mu_r(ipoint,istate)
  accu_f += dabs(f_ref - f_comp) * weight
  accu_n2+= dabs(on_top_ref - on_top_comp) * weight
  accu += f_ref * weight
 enddo
 print*,'**************************'
 print*,''
 print*,'accu_f  = ',accu_f
 print*,'accu_n2 = ',accu_n2
 print*,''
 print*,'accu    = ',accu

end
