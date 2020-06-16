
 BEGIN_PROVIDER [double precision, one_e_act_dm_beta_mo_for_dft, (n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
 ! one_e_act_dm_beta_mo_for_dft = pure ACTIVE part of the ONE ELECTRON REDUCED DENSITY MATRIX for the BETA ELECTRONS 
 END_DOC
 integer :: i,j,ii,jj,istate
 do istate = 1, N_states
  do ii = 1, n_act_orb
   i = list_act(ii)
   do jj = 1, n_act_orb
    j = list_act(jj)
    one_e_act_dm_beta_mo_for_dft(jj,ii,istate) = one_e_dm_mo_beta_for_dft(j,i,istate)
   enddo
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, one_e_act_dm_alpha_mo_for_dft, (n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
 ! one_e_act_dm_alpha_mo_for_dft = pure ACTIVE part of the ONE ELECTRON REDUCED DENSITY MATRIX for the ALPHA ELECTRONS 
 END_DOC
 integer :: i,j,ii,jj,istate
 do istate = 1, N_states
  do ii = 1, n_act_orb
   i = list_act(ii)
   do jj = 1, n_act_orb
    j = list_act(jj)
    one_e_act_dm_alpha_mo_for_dft(jj,ii,istate) = one_e_dm_mo_alpha_for_dft(j,i,istate)
   enddo
  enddo
 enddo

END_PROVIDER 

