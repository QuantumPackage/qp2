 BEGIN_PROVIDER [double precision, multi_s_deriv_1, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_x_deriv_1, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_y_deriv_1, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_z_deriv_1, (N_states, N_states)]

  implicit none

  BEGIN_DOC
  ! Providers for :
  ! <Psi_m|v_x|Psi_n>
  ! <Psi_m|v_y|Psi_n>
  ! <Psi_m|v_z|Psi_n>
  ! ||v|| = sqrt(v_x^2 + v_y^2 + v_z^2)
  ! v_x = d/dx
  ! Cf. multi_s_dipole_moment for the equations
  END_DOC

  integer          :: istate,jstate ! States
  integer          :: i,j           ! general spatial MOs
  double precision :: nuclei_part_x, nuclei_part_y, nuclei_part_z
 
  multi_s_x_deriv_1 = 0.d0
  multi_s_y_deriv_1 = 0.d0
  multi_s_z_deriv_1 = 0.d0
 
  do jstate = 1, N_states
    do istate = 1, N_states
 
      do i = 1, mo_num  
        ! Diag part
        multi_s_x_deriv_1(istate,jstate) -= one_e_tr_dm_mo(i,i,istate,jstate) * mo_deriv_1_x(i,i)
        multi_s_y_deriv_1(istate,jstate) -= one_e_tr_dm_mo(i,i,istate,jstate) * mo_deriv_1_y(i,i)
        multi_s_z_deriv_1(istate,jstate) -= one_e_tr_dm_mo(i,i,istate,jstate) * mo_deriv_1_z(i,i)
 
        do j = 1, mo_num  
          if (i == j) then
           cycle
          endif
          ! Extra diag part
          multi_s_x_deriv_1(istate,jstate) -= one_e_tr_dm_mo(j,i,istate,jstate) * mo_deriv_1_x(j,i)  
          multi_s_y_deriv_1(istate,jstate) -= one_e_tr_dm_mo(j,i,istate,jstate) * mo_deriv_1_y(j,i) 
          multi_s_z_deriv_1(istate,jstate) -= one_e_tr_dm_mo(j,i,istate,jstate) * mo_deriv_1_z(j,i) 
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
    multi_s_x_deriv_1(istate,istate) += nuclei_part_x
    multi_s_y_deriv_1(istate,istate) += nuclei_part_y
    multi_s_z_deriv_1(istate,istate) += nuclei_part_z
  enddo
 
  ! d = <Psi|r|Psi>
  do jstate = 1, N_states
    do istate = 1, N_states
      multi_s_deriv_1(istate,jstate) = &
        dsqrt(multi_s_x_deriv_1(istate,jstate)**2 & 
            + multi_s_y_deriv_1(istate,jstate)**2 &
            + multi_s_z_deriv_1(istate,jstate)**2) 
    enddo
  enddo

END_PROVIDER

