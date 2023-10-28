
BEGIN_PROVIDER [double precision, M_RPA, (2*nS_exc, 2*nS_exc)]

  BEGIN_DOC
  !
  ! full matrix for direct RPA calculation 
  ! with the TC-Hamiltonian
  !
  END_DOC

  implicit none
  integer                    :: ia, i, a, jb, j, b
  double precision           :: e(mo_num)
  double precision, external :: Kronecker_delta

  PROVIDE mo_tc_effec2e_int
  PROVIDE Fock_matrix_tc_diag_mo_tot

  e(1:mo_num) = Fock_matrix_tc_diag_mo_tot(1:mo_num)


  ! --- --- ---
  ! block A

  ia = 0
  do i = nC_orb+1, nO_orb
    do a = nO_orb+1, mo_num-nR_orb
      ia = ia + 1

      jb = 0
      do j = nC_orb+1, nO_orb
        do b = nO_orb+1, mo_num-nR_orb
          jb = jb + 1
  
          M_RPA(ia,jb) = (e(a) - e(i)) * Kronecker_delta(i,j) * Kronecker_delta(a,b) + 2.d0 * mo_tc_effec2e_int(a,j,i,b)
        enddo
      enddo
    enddo
  enddo

  !
  ! --- --- ---


  ! --- --- ---
  ! block B

  ia = 0
  do i = nC_orb+1, nO_orb
    do a = nO_orb+1, mo_num-nR_orb
      ia = ia + 1

      jb = nS_exc
      do j = nC_orb+1, nO_orb
        do b = nO_orb+1, mo_num-nR_orb
          jb = jb + 1
  
          M_RPA(ia,jb) = 2.d0 * mo_tc_effec2e_int(a,b,i,j)
        enddo
      enddo
    enddo
  enddo

  !
  ! --- --- ---


  ! --- --- ---
  ! block C

  ia = nS_exc
  do i = nC_orb+1, nO_orb
    do a = nO_orb+1, mo_num-nR_orb
      ia = ia + 1

      jb = 0
      do j = nC_orb+1, nO_orb
        do b = nO_orb+1, mo_num-nR_orb
          jb = jb + 1
  
          M_RPA(ia,jb) = 2.d0 * mo_tc_effec2e_int(i,j,a,b)
        enddo
      enddo
    enddo
  enddo

  !
  ! --- --- ---


  ! --- --- ---
  ! block D

  ia = nS_exc
  do i = nC_orb+1, nO_orb
    do a = nO_orb+1, mo_num-nR_orb
      ia = ia + 1

      jb = nS_exc
      do j = nC_orb+1, nO_orb
        do b = nO_orb+1, mo_num-nR_orb
          jb = jb + 1
  
          M_RPA(ia,jb) = (e(a) - e(i)) * Kronecker_delta(i,j) * Kronecker_delta(a,b) + 2.d0 * mo_tc_effec2e_int(i,b,a,j)
        enddo
      enddo
    enddo
  enddo

  !
  ! --- --- ---


END_PROVIDER


