
BEGIN_PROVIDER [ double precision, ao_integrals_pt_chrg, (ao_num,ao_num)]

  BEGIN_DOC
  !  Point charge-electron interaction, in the |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_charge charge * \frac{1}{|r-R_charge|} | \chi_j \rangle`
  !
  ! Notice the minus sign convention as it is supposed to be for electrons. 
  END_DOC

  implicit none
  integer          :: num_A, num_B, power_A(3), power_B(3)
  integer          :: i, j, k, l, n_pt_in, m
  double precision :: alpha, beta
  double precision :: A_center(3),B_center(3),C_center(3)
  double precision :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult

  ao_integrals_pt_chrg = 0.d0

!  if (read_ao_integrals_pt_chrg) then
!
!    call ezfio_get_ao_one_e_ints_ao_integrals_pt_chrg(ao_integrals_pt_chrg)
!    print *,  'AO N-e integrals read from disk'
!
!  else

!    if(use_cosgtos) then
!      !print *, " use_cosgtos for ao_integrals_pt_chrg ?", use_cosgtos
!
!      do j = 1, ao_num
!        do i = 1, ao_num
!          ao_integrals_pt_chrg(i,j) = ao_integrals_pt_chrg_cosgtos(i,j)
!        enddo
!      enddo
!
!    else

      !$OMP PARALLEL                                                   &
          !$OMP DEFAULT (NONE)                                         &
          !$OMP PRIVATE (i,j,k,l,m,alpha,beta,A_center,B_center,C_center,power_A,power_B,&
          !$OMP          num_A,num_B,Z,c,c1,n_pt_in)                      &
          !$OMP SHARED (ao_num,ao_prim_num,ao_expo_ordered_transp,ao_power,ao_nucl,pts_charge_coord,ao_coef_normalized_ordered_transp,nucl_coord,&
          !$OMP         n_pt_max_integrals,ao_integrals_pt_chrg,n_pts_charge,pts_charge_z)

      n_pt_in = n_pt_max_integrals

      !$OMP DO SCHEDULE (dynamic)

      do j = 1, ao_num
        num_A = ao_nucl(j)
        power_A(1:3)= ao_power(j,1:3)
        A_center(1:3) = nucl_coord(num_A,1:3)

        do i = 1, ao_num

          num_B = ao_nucl(i)
          power_B(1:3)= ao_power(i,1:3)
          B_center(1:3) = nucl_coord(num_B,1:3)

          do l=1,ao_prim_num(j)
            alpha = ao_expo_ordered_transp(l,j)

            do m=1,ao_prim_num(i)
              beta = ao_expo_ordered_transp(m,i)

              double precision               :: c, c1
              c = 0.d0

              do  k = 1, n_pts_charge
                double precision               :: Z
                Z = pts_charge_z(k)

                C_center(1:3) = pts_charge_coord(k,1:3)

                c1 = NAI_pol_mult( A_center, B_center, power_A, power_B &
                                 , alpha, beta, C_center, n_pt_in )

                c = c - Z * c1

              enddo
              ao_integrals_pt_chrg(i,j) = ao_integrals_pt_chrg(i,j)  &
                  + ao_coef_normalized_ordered_transp(l,j)             &
                  * ao_coef_normalized_ordered_transp(m,i) * c
            enddo
          enddo
        enddo
      enddo

    !$OMP END DO
    !$OMP END PARALLEL

!    endif


!    IF(do_pseudo) THEN
!       ao_integrals_pt_chrg += ao_pseudo_integrals
!    ENDIF

!  endif


!  if (write_ao_integrals_pt_chrg) then
!    call ezfio_set_ao_one_e_ints_ao_integrals_pt_chrg(ao_integrals_pt_chrg)
!    print *,  'AO N-e integrals written to disk'
!  endif

END_PROVIDER
