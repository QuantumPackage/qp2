
 BEGIN_PROVIDER [ double precision, tc_transition_matrix_mo_beta, (mo_num, mo_num,N_states,N_states) ]
&BEGIN_PROVIDER [ double precision, tc_transition_matrix_mo_alpha, (mo_num, mo_num,N_states,N_states) ]
 implicit none
 BEGIN_DOC
 ! tc_transition_matrix_mo_alpha(p,h,istate,jstate) = <Chi_istate| a^\dagger_p,alpha a_h,alpha |Phi_jstate>
 !
 ! tc_transition_matrix_mo_beta(p,h,istate,jstate) = <Chi_istate| a^\dagger_p,beta a_h,beta |Phi_jstate>
 !
 ! where <Chi_istate| and |Phi_jstate> are the left/right eigenvectors on a bi-ortho basis
 END_DOC
 integer :: i,j,istate,jstate,m,n,p,h
 double precision :: phase
 integer, allocatable           :: occ(:,:)
 integer                        :: n_occ_ab(2),degree,exc(0:2,2,2)
 allocate(occ(N_int*bit_kind_size,2))
 tc_transition_matrix_mo_alpha = 0.d0
 tc_transition_matrix_mo_beta  = 0.d0
   do i = 1, N_det
    do j = 1, N_det
     call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
     if(degree.gt.1)cycle
     do istate = 1, N_states
      do jstate = 1, N_states
       if (degree == 0)then
        call bitstring_to_list_ab(psi_det(1,1,i), occ, n_occ_ab, N_int)
        do p = 1, n_occ_ab(1) ! browsing the alpha electrons
         m = occ(p,1)
         tc_transition_matrix_mo_alpha(m,m,istate,jstate)+= psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,jstate)
        enddo
        do p = 1, n_occ_ab(2) ! browsing the beta electrons
         m = occ(p,2)
         tc_transition_matrix_mo_beta(m,m,istate,jstate)+= psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,jstate)
        enddo
       else
        call get_single_excitation(psi_det(1,1,j),psi_det(1,1,i),exc,phase,N_int)
        if (exc(0,1,1) == 1) then
          ! Single alpha
          h = exc(1,1,1) ! hole in psi_det(1,1,j) 
          p = exc(1,2,1) ! particle in psi_det(1,1,j) 
          tc_transition_matrix_mo_alpha(p,h,istate,jstate)+= &
          phase * psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,jstate)
        else
          ! Single beta
          h = exc(1,1,2) ! hole in psi_det(1,1,j) 
          p = exc(1,2,2) ! particle in psi_det(1,1,j) 
          tc_transition_matrix_mo_beta(p,h,istate,jstate)+=  &
          phase * psi_l_coef_bi_ortho(i,istate) * psi_r_coef_bi_ortho(j,jstate)
        endif
       endif
     enddo
    enddo
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER [double precision, tc_transition_matrix_mo, (mo_num, mo_num,N_states,N_states) ]
 implicit none
 BEGIN_DOC
 ! tc_transition_matrix_mo(p,h,istate,jstate) = \sum_{sigma=alpha,beta} <Chi_istate| a^\dagger_p,sigma a_h,sigma |Phi_jstate>
 !
 ! where <Chi_istate| and |Phi_jstate> are the left/right eigenvectors on a bi-ortho basis
 END_DOC
 tc_transition_matrix_mo = tc_transition_matrix_mo_beta + tc_transition_matrix_mo_alpha
 END_PROVIDER 
 

 BEGIN_PROVIDER [double precision, tc_spin_transition_matrix_mo, (mo_num, mo_num,N_states,N_states) ]
 implicit none
 BEGIN_DOC
 ! tc_spin_transition_matrix_mo = tc_transition_matrix_mo_alpha - tc_transition_matrix_mo_beta
 !
 ! where <Chi_istate| and |Phi_jstate> are the left/right eigenvectors on a bi-ortho basis
 END_DOC
 tc_spin_transition_matrix_mo = tc_transition_matrix_mo_alpha - tc_transition_matrix_mo_beta 
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, tc_bi_ortho_dipole, (3,N_states)]
 implicit none
 integer :: i,j,istate,m
 double precision :: nuclei_part(3)
 tc_bi_ortho_dipole = 0.d0
 do istate = 1, N_states
  do i = 1, mo_num
   do j = 1, mo_num
    tc_bi_ortho_dipole(1,istate) += -(tc_transition_matrix_mo(j,i,istate,istate)) *  mo_bi_orth_bipole_x(j,i)
    tc_bi_ortho_dipole(2,istate) += -(tc_transition_matrix_mo(j,i,istate,istate)) *  mo_bi_orth_bipole_y(j,i)
    tc_bi_ortho_dipole(3,istate) += -(tc_transition_matrix_mo(j,i,istate,istate)) *  mo_bi_orth_bipole_z(j,i)
   enddo
  enddo
 enddo
 print*,'tc_bi_ortho_dipole(3) elec = ',tc_bi_ortho_dipole(3,1)

 nuclei_part = 0.d0
 do m = 1, 3
  do i = 1,nucl_num
   nuclei_part(m) += nucl_charge(i) * nucl_coord(i,m)
  enddo
 enddo
!
 do istate = 1, N_states
  do m = 1, 3
    tc_bi_ortho_dipole(m,istate) += nuclei_part(m)
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER [ double precision, tc_transition_matrix_ao, (ao_num, ao_num,N_states,N_states) ]
 implicit none
 BEGIN_DOC 
! tc_transition_matrix(p,h,istate,jstate) in the AO basis  
 END_DOC 
 integer :: i,j,k,l
 double precision :: dm_mo
 tc_transition_matrix_ao = 0.d0
 integer :: istate,jstate
 do istate = 1, N_states
  do jstate = 1, N_states
    do i = 1, mo_num
      do j = 1, mo_num
        dm_mo = tc_transition_matrix_mo(j,i,jstate,istate)
        do k = 1, ao_num
          do l = 1, ao_num
            tc_transition_matrix_ao(l,k,jstate,istate) += mo_l_coef(l,j) * mo_r_coef(k,i) * dm_mo
          enddo
        enddo
      enddo
    enddo
  enddo
 enddo
 
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, tc_spin_transition_matrix_ao, (ao_num, ao_num,N_states,N_states) ]
 implicit none
 BEGIN_DOC 
! tc_spin_transition_matrix_ao(p,h,istate,jstate) in the AO basis  
 END_DOC 
 integer :: i,j,k,l
 double precision :: dm_mo
 tc_spin_transition_matrix_ao = 0.d0
 integer :: istate,jstate
 do istate = 1, N_states
  do jstate = 1, N_states
    do i = 1, mo_num
      do j = 1, mo_num
        dm_mo = tc_spin_transition_matrix_mo(j,i,jstate,istate)
        do k = 1, ao_num
          do l = 1, ao_num
            tc_spin_transition_matrix_ao(l,k,jstate,istate) += mo_l_coef(l,j) * mo_r_coef(k,i) * dm_mo
          enddo
        enddo
      enddo
    enddo
  enddo
 enddo
 
 END_PROVIDER 
