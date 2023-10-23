! Providers for the dipole moments along x,y,z and the total dipole
! moments.

! The dipole moment along the x axis is:
! \begin{align*}
! \mu_x = < \Psi_m | \sum_i x_i + \sum_A Z_A R_A | \Psi_n >
! \end{align*}
! where $i$ is used for the electrons and $A$ for the nuclei.
! $Z_A$ the charge of the nucleus $A$ and $R_A$ its position in the
! space.

! And it can be computed using the (transition, if n /= m) density
! matrix as a expectation value
! \begin{align*}
! <\Psi_n|x| \Psi_m > = \sum_p \gamma_{pp}^{nm} < \phi_p | x | \phi_p >
!      + \sum_{pq, p \neq q} \gamma_{pq}^{nm} < \phi_p |x | \phi_q > +  < \Psi_m | \sum_A Z_A R_A | \Psi_n >
! \end{align*}



BEGIN_PROVIDER [double precision, multi_s_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_x_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_y_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_z_dipole_moment, (N_states, N_states)]

  implicit none

  BEGIN_DOC
  ! Providers for :
  ! <\Psi_m|\mu_x|\Psi_n>
  ! <\Psi_m|\mu_y|\Psi_n>
  ! <\Psi_m|\mu_z|\Psi_n>
  ! ||\mu|| = \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}
  !
  ! <\Psi_n|x| \Psi_m > = \sum_p \gamma_{pp}^{nm} \bra{\phi_p} x \ket{\phi_p} 
  !   + \sum_{pq, p \neq q} \gamma_{pq}^{nm} \bra{\phi_p} x \ket{\phi_q}
  ! \Psi: wf
  ! n,m indexes for the states
  ! p,q: general spatial MOs 
  ! gamma^{nm}: density matrix \bra{\Psi^n} a^{\dagger}_a a_i \ket{\Psi^m}
  END_DOC

  integer          :: istate,jstate ! States
  integer          :: i,j           ! general spatial MOs
  double precision :: nuclei_part_x, nuclei_part_y, nuclei_part_z
 
  multi_s_x_dipole_moment = 0.d0
  multi_s_y_dipole_moment = 0.d0
  multi_s_z_dipole_moment = 0.d0
 
  do jstate = 1, N_states
    do istate = 1, N_states
 
      do i = 1, mo_num  
        do j = 1, mo_num  
          multi_s_x_dipole_moment(istate,jstate) -= one_e_tr_dm_mo(j,i,istate,jstate) * mo_dipole_x(j,i)  
          multi_s_y_dipole_moment(istate,jstate) -= one_e_tr_dm_mo(j,i,istate,jstate) * mo_dipole_y(j,i) 
          multi_s_z_dipole_moment(istate,jstate) -= one_e_tr_dm_mo(j,i,istate,jstate) * mo_dipole_z(j,i) 
        enddo
      enddo
 
    enddo
  enddo
 
  ! Nuclei part
  nuclei_part_x = 0.d0
  nuclei_part_y = 0.d0
  nuclei_part_z = 0.d0
 
  do i = 1,nucl_num 
    nuclei_part_x += nucl_charge(i) * nucl_coord(i,1) 
    nuclei_part_y += nucl_charge(i) * nucl_coord(i,2) 
    nuclei_part_z += nucl_charge(i) * nucl_coord(i,3) 
  enddo
 
  ! Only if istate = jstate, otherwise 0 by the orthogonality of the states
  do istate = 1, N_states
    multi_s_x_dipole_moment(istate,istate) += nuclei_part_x
    multi_s_y_dipole_moment(istate,istate) += nuclei_part_y
    multi_s_z_dipole_moment(istate,istate) += nuclei_part_z
  enddo
 
  ! d = <Psi|r|Psi>
  do jstate = 1, N_states
    do istate = 1, N_states
      multi_s_dipole_moment(istate,jstate) = &
        dsqrt(multi_s_x_dipole_moment(istate,jstate)**2 & 
            + multi_s_y_dipole_moment(istate,jstate)**2 &
            + multi_s_z_dipole_moment(istate,jstate)**2) 
    enddo
  enddo

END_PROVIDER
