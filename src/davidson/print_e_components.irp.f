subroutine print_energy_components()
  implicit none
  BEGIN_DOC
! Prints the different components of the energy.
  END_DOC
  integer, save :: ifirst = 0
  double precision :: Vee, Ven, Vnn, Vecp, T, f
  complex*16 :: fc
  integer  :: i,j,k,kk

  Vnn = nuclear_repulsion

  print *, 'Energy components'
  print *, '================='
  print *, ''
  do k=1,N_states

    Ven  = 0.d0
    Vecp = 0.d0
    T    = 0.d0
    
    if (is_complex) then
      do kk=1,kpt_num
        do j=1,mo_num_per_kpt
          do i=1,mo_num_per_kpt
            !fc = one_e_dm_mo_alpha_complex(i,j,k) + one_e_dm_mo_beta_complex(i,j,k)
            !Ven  = Ven  + dble(fc * mo_integrals_n_e_complex(j,i))
            !Vecp = Vecp + dble(fc * mo_pseudo_integrals_complex(j,i))
            !T    = T    + dble(fc * mo_kinetic_integrals_complex(j,i))
            fc = one_e_dm_mo_alpha_kpts(i,j,kk,k) + one_e_dm_mo_beta_kpts(i,j,kk,k)
            Ven  = Ven  + dble(fc * mo_integrals_n_e_kpts(j,i,kk))
            Vecp = Vecp + dble(fc * mo_pseudo_integrals_kpts(j,i,kk))
            T    = T    + dble(fc * mo_kinetic_integrals_kpts(j,i,kk))
          enddo
        enddo
      enddo
    else
      do j=1,mo_num
        do i=1,mo_num
          f = one_e_dm_mo_alpha(i,j,k) + one_e_dm_mo_beta(i,j,k)
          Ven  = Ven  + f * mo_integrals_n_e(i,j)
          Vecp = Vecp + f * mo_pseudo_integrals(i,j)
          T    = T    + f * mo_kinetic_integrals(i,j)
        enddo
      enddo
    endif
    Vee = psi_energy(k) - Ven - Vecp - T
    
    if (ifirst == 0) then
      ifirst = 1
      print *, 'Vnn  : Nucleus-Nucleus   potential energy'
      print *, 'Ven  : Electron-Nucleus  potential energy'
      print *, 'Vee  : Electron-Electron potential energy'
      print *, 'Vecp : Potential energy of the pseudo-potentials'
      print *, 'T    : Electronic kinetic energy'
      print *, ''
    endif

    print *, 'State ', k
    print *, '---------'
    print *, ''
    print *, 'Vnn  = ', Vnn
    print *, 'Ven  = ', Ven
    print *, 'Vee  = ', Vee
    print *, 'Vecp = ', Vecp
    print *, 'T    = ', T
    print *, ''
  enddo

  print *, ''

end
