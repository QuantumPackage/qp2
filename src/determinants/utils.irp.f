BEGIN_PROVIDER [ double precision, H_matrix_all_dets,(N_det,N_det) ]
  use bitmasks
 implicit none
 BEGIN_DOC
 ! |H| matrix on the basis of the Slater determinants defined by psi_det
 END_DOC
 integer :: i,j,k
 double precision :: hij
 call  i_H_j(psi_det(1,1,1),psi_det(1,1,1),N_int,hij)
 print*,'Providing the H_matrix_all_dets ...'
 !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) PRIVATE(i,j,hij,k) &
 !$OMP SHARED (N_det, psi_det, N_int,H_matrix_all_dets)
 do i =1,N_det
   do j = i, N_det
    call  i_H_j(psi_det(1,1,i),psi_det(1,1,j),N_int,hij)
    H_matrix_all_dets(i,j) = hij
    H_matrix_all_dets(j,i) = hij
  enddo
 enddo
 !$OMP END PARALLEL DO
 print*,'H_matrix_all_dets done '
END_PROVIDER

BEGIN_PROVIDER [ double precision, H_matrix_diag_all_dets,(N_det) ]
  use bitmasks
 implicit none
 BEGIN_DOC
 ! |H| matrix on the basis of the Slater determinants defined by psi_det
 END_DOC
 integer :: i
 double precision :: hij

 call  i_H_j(psi_det(1,1,1),psi_det(1,1,1),N_int,hij)
 !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) PRIVATE(i,hij) &
 !$OMP SHARED (N_det, psi_det, N_int,H_matrix_diag_all_dets)
 do i =1,N_det
   call  i_H_j(psi_det(1,1,i),psi_det(1,1,i),N_int,hij)
   H_matrix_diag_all_dets(i) = hij
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
 call get_s2(psi_det(1,1,1),psi_det(1,1,1),N_int,sij)
 !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) PRIVATE(i,j,sij,k) &
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
