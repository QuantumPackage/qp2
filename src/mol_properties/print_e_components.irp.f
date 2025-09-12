subroutine print_energy_components()
  implicit none
  BEGIN_DOC
! Prints the different components of the energy.
  END_DOC
  integer, save :: ifirst = 0
  double precision :: Vee, Ven, Vnn, Vecp, T, f
  integer  :: i,j,k

  Vnn = nuclear_repulsion

  print *, 'Energy components'
  print *, '================='
  print *, ''

  print *, 'Vnn  : Nucleus-Nucleus   potential energy'
  print *, 'Ven  : Electron-Nucleus  potential energy'
  print *, 'Vee  : Electron-Electron potential energy'
  print *, 'T    : Electronic kinetic energy'
  print *, ''

  do k=1,N_states

    Ven  = 0.d0
    Vecp = 0.d0
    T    = 0.d0

    do j=1,mo_num
      do i=1,mo_num
        f = one_e_dm_mo_alpha(i,j,k) + one_e_dm_mo_beta(i,j,k)
        Ven  = Ven  + f * mo_integrals_n_e(i,j)
        T    = T    + f * mo_kinetic_integrals(i,j)
      enddo
    enddo
    Vee = psi_energy(k) - Ven - T

    print *, 'State ', k
    print *, '---------'
    print *, ''
    print *, 'Vnn  = ', Vnn
    print *, 'Ven  = ', Ven
    print *, 'Vee  = ', Vee
    print *, 'T    = ', T
    print *, ''
  enddo

  print *, ''

end
