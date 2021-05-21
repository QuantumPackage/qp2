BEGIN_PROVIDER [ double precision, basis_normalization_factor, (shell_num) ]
  implicit none
  BEGIN_DOC
  ! Normalization factors of the shells
  END_DOC
  double precision               :: norm,overlap_x,overlap_y,overlap_z,C_A(3), c
  integer                        :: l, powA(3), nz
  integer                        :: i,j,k
  nz=100
  C_A(1) = 0.d0
  C_A(2) = 0.d0
  C_A(3) = 0.d0

  do i=1,shell_num
    powA(1) = shell_ang_mom(i)
    powA(2) = 0
    powA(3) = 0

    ! Normalization of the contracted basis functions
    norm = 0.d0
    do j=shell_prim_index(i), shell_prim_index(i)+shell_prim_num(i)-1
     do k=shell_prim_index(i), shell_prim_index(i)+shell_prim_num(i)-1
      call overlap_gaussian_xyz(C_A,C_A,shell_prim_expo(j),shell_prim_expo(k), &
          powA,powA,overlap_x,overlap_y,overlap_z,c,nz)
      norm = norm+c*shell_prim_coef(j)*shell_prim_coef(k)
     enddo
    enddo
    basis_normalization_factor(i) = basis_normalization_factor(i) * sqrt(norm)

  enddo

END_PROVIDER
