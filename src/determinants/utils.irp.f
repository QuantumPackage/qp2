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

BEGIN_PROVIDER [ complex*16, h_matrix_all_dets_complex,(N_det,N_det) ]
  use bitmasks
 implicit none
 BEGIN_DOC
 ! |H| matrix on the basis of the Slater determinants defined by psi_det
 END_DOC
 integer :: i,j,k
 complex*16 :: hij
 integer :: degree(N_det),idx(0:N_det)
 call  i_h_j_complex(psi_det(1,1,1),psi_det(1,1,1),N_int,hij)
 !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) PRIVATE(i,j,hij,degree,idx,k) &
 !$OMP SHARED (N_det, psi_det, N_int,h_matrix_all_dets_complex)
 do i =1,N_det
   do j = i, N_det
    call  i_h_j_complex(psi_det(1,1,i),psi_det(1,1,j),N_int,hij)
    H_matrix_all_dets_complex(i,j) = hij
    H_matrix_all_dets_complex(j,i) = dconjg(hij)
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

