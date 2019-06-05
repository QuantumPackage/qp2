BEGIN_PROVIDER [ double precision, H_matrix_all_dets,(N_det,N_det) ]
  use bitmasks
 implicit none
 BEGIN_DOC
 ! |H| matrix on the basis of the Slater determinants defined by psi_det
 END_DOC
 integer :: i,j,k
 double precision :: hij
 integer :: degree(N_det),idx(0:N_det)
 call  i_H_j(psi_det(1,1,1),psi_det(1,1,1),N_int,hij)
 !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) PRIVATE(i,j,hij,degree,idx,k) &
 !$OMP SHARED (N_det, psi_det, N_int,H_matrix_all_dets)
 do i =1,N_det
   do j = i, N_det
    call  i_H_j(psi_det(1,1,i),psi_det(1,1,j),N_int,hij)
    H_matrix_all_dets(i,j) = hij
    H_matrix_all_dets(j,i) = hij
  enddo
 enddo
 !$OMP END PARALLEL DO
END_PROVIDER


BEGIN_PROVIDER [ double precision, S2_matrix_all_dets,(N_det,N_det) ]
  use bitmasks
 implicit none
 BEGIN_DOC
 ! |S^2| matrix on the basis of the Slater determinants defined by psi_det
 END_DOC
 integer :: i,j,k
 double precision :: sij
 integer :: degree(N_det),idx(0:N_det)
 call get_s2(psi_det(1,1,1),psi_det(1,1,1),N_int,sij)
 !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) PRIVATE(i,j,sij,degree,idx,k) &
 !$OMP SHARED (N_det, psi_det, N_int,S2_matrix_all_dets)
 do i =1,N_det
   do j = i, N_det
    call get_s2(psi_det(1,1,i),psi_det(1,1,j),N_int,sij)
    S2_matrix_all_dets(i,j) = sij
    S2_matrix_all_dets(j,i) = sij
  enddo
 enddo
 !$OMP END PARALLEL DO
END_PROVIDER

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
  do k=1,N_states

    Ven  = 0.d0
    Vecp = 0.d0
    T    = 0.d0

    do j=1,mo_num
      do i=1,mo_num
        f = one_e_dm_mo_alpha(i,j,k) + one_e_dm_mo_beta(i,j,k)
        Ven  = Ven  + f * mo_integrals_n_e(i,j)
        Vecp = Vecp + f * mo_pseudo_integrals(i,j)
        T    = T    + f * mo_kinetic_integrals(i,j)
      enddo
    enddo
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
