
BEGIN_PROVIDER [ double precision, CI_energy, (N_states_diag) ]
  implicit none
  BEGIN_DOC
  ! :c:data:`n_states` lowest eigenvalues of the |CI| matrix
  END_DOC

  integer                        :: j
  character*(8)                  :: st
  call write_time(6)
  do j=1,min(N_det,N_states_diag)
    CI_energy(j) = CI_electronic_energy(j) + nuclear_repulsion
  enddo
  do j=1,min(N_det,N_states)
    write(st,'(I4)') j
    call write_double(6,CI_energy(j),'Energy of state '//trim(st))
    call write_double(6,CI_s2(j),'S^2 of state '//trim(st))
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, CI_electronic_energy, (N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_s2, (N_states_diag) ]
  implicit none
  if (is_complex) then
    ci_s2(1:N_states_diag) = ci_s2_complex(1:N_states_diag)
    ci_electronic_energy(1:N_states_diag) = ci_electronic_energy_complex(1:N_states_diag)
  else
    ci_s2(1:N_states_diag) = ci_s2_real(1:N_states_diag)
    ci_electronic_energy(1:N_states_diag) = ci_electronic_energy_real(1:N_states_diag)
  endif
END_PROVIDER


 BEGIN_PROVIDER [ double precision, CI_electronic_energy_real, (N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_eigenvectors, (N_det,N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_s2_real, (N_states_diag) ]
   BEGIN_DOC
   ! Eigenvectors/values of the |CI| matrix
   END_DOC
   implicit none
   double precision               :: ovrlp,u_dot_v
   integer                        :: i_good_state
   integer, allocatable           :: index_good_state_array(:)
   logical, allocatable           :: good_state_array(:)
   double precision, allocatable  :: s2_values_tmp(:)
   integer                        :: i_other_state
   double precision, allocatable  :: eigenvectors(:,:), eigenvalues(:), H_prime(:,:)
   integer                        :: i_state
   double precision               :: e_0
   integer                        :: i,j,k
   double precision, allocatable  :: s2_eigvalues(:)
   double precision, allocatable  :: e_array(:)
   integer, allocatable           :: iorder(:)
   logical                        :: converged

   PROVIDE threshold_davidson nthreads_davidson
   ! Guess values for the "N_states" states of the |CI| eigenvectors
   do j=1,min(N_states,N_det)
     do i=1,N_det
       CI_eigenvectors(i,j) = psi_coef(i,j)
     enddo
   enddo

   do j=min(N_states,N_det)+1,N_states_diag
     do i=1,N_det
       CI_eigenvectors(i,j) = 0.d0
     enddo
   enddo

   if (diag_algorithm == "Davidson") then

     call davidson_diag_HS2(psi_det,CI_eigenvectors, CI_s2_real, &
         size(CI_eigenvectors,1),CI_electronic_energy_real,               &
         N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,0,converged)

     integer :: N_states_diag_save
     N_states_diag_save = N_states_diag
     do while (.not.converged)
        double precision, allocatable :: CI_electronic_energy_tmp (:)
        double precision, allocatable :: CI_eigenvectors_tmp (:,:)
        double precision, allocatable :: CI_s2_tmp (:)

        N_states_diag  *= 2
        TOUCH N_states_diag

        allocate (CI_electronic_energy_tmp (N_states_diag) )
        allocate (CI_eigenvectors_tmp (N_det,N_states_diag) )
        allocate (CI_s2_tmp (N_states_diag) )

        CI_electronic_energy_tmp(1:N_states_diag_save) = CI_electronic_energy_real(1:N_states_diag_save)
        CI_eigenvectors_tmp(1:N_det,1:N_states_diag_save) = CI_eigenvectors(1:N_det,1:N_states_diag_save)
        CI_s2_tmp(1:N_states_diag_save) = CI_s2_real(1:N_states_diag_save)

        call davidson_diag_HS2(psi_det,CI_eigenvectors_tmp, CI_s2_tmp, &
            size(CI_eigenvectors_tmp,1),CI_electronic_energy_tmp,               &
            N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,0,converged)

        CI_electronic_energy_real(1:N_states_diag_save) = CI_electronic_energy_tmp(1:N_states_diag_save)
        CI_eigenvectors(1:N_det,1:N_states_diag_save) = CI_eigenvectors_tmp(1:N_det,1:N_states_diag_save)
        CI_s2_real(1:N_states_diag_save) = CI_s2_tmp(1:N_states_diag_save)

        deallocate (CI_electronic_energy_tmp)
        deallocate (CI_eigenvectors_tmp)
        deallocate (CI_s2_tmp)
     enddo
     if (N_states_diag > N_states_diag_save) then
       N_states_diag = N_states_diag_save
       TOUCH N_states_diag
     endif

   else if (diag_algorithm == "Lapack") then

     print *,  'Diagonalization of H using Lapack'
     allocate (eigenvectors(size(H_matrix_all_dets,1),N_det))
     allocate (eigenvalues(N_det))
     if (s2_eig) then
       double precision, parameter :: alpha = 0.1d0
       allocate (H_prime(N_det,N_det) )
       H_prime(1:N_det,1:N_det) = H_matrix_all_dets(1:N_det,1:N_det) +  &
         alpha * S2_matrix_all_dets(1:N_det,1:N_det)
       do j=1,N_det
         H_prime(j,j) = H_prime(j,j) - alpha*expected_s2
       enddo
       call lapack_diag(eigenvalues,eigenvectors,H_prime,size(H_prime,1),N_det)
       CI_electronic_energy_real(:) = 0.d0
       i_state = 0
       allocate (s2_eigvalues(N_det))
       allocate(index_good_state_array(N_det),good_state_array(N_det))
       good_state_array = .False.
       call u_0_S2_u_0(s2_eigvalues,eigenvectors,N_det,psi_det,N_int,&
         N_det,size(eigenvectors,1))
       if (only_expected_s2) then
         do j=1,N_det
           ! Select at least n_states states with S^2 values closed to "expected_s2"
           if(dabs(s2_eigvalues(j)-expected_s2).le.0.5d0)then
             i_state +=1
             index_good_state_array(i_state) = j
             good_state_array(j) = .True.
           endif
           if(i_state.eq.N_states) then
             exit
           endif
         enddo
       else
         do j=1,N_det
           index_good_state_array(j) = j
           good_state_array(j) = .True.
         enddo
       endif
       if(i_state .ne.0)then
         ! Fill the first "i_state" states that have a correct S^2 value
         do j = 1, i_state
           do i=1,N_det
             CI_eigenvectors(i,j) = eigenvectors(i,index_good_state_array(j))
           enddo
           CI_electronic_energy_real(j) = eigenvalues(index_good_state_array(j))
           CI_s2_real(j) = s2_eigvalues(index_good_state_array(j))
         enddo
         i_other_state = 0
         do j = 1, N_det
           if(good_state_array(j))cycle
           i_other_state +=1
           if(i_state+i_other_state.gt.n_states_diag)then
             exit
           endif
           do i=1,N_det
             CI_eigenvectors(i,i_state+i_other_state) = eigenvectors(i,j)
           enddo
           CI_electronic_energy_real(i_state+i_other_state) = eigenvalues(j)
           CI_s2_real(i_state+i_other_state) = s2_eigvalues(i_state+i_other_state)
         enddo

       else
         print*,''
         print*,'!!!!!!!!   WARNING  !!!!!!!!!'
         print*,'  Within the ',N_det,'determinants selected'
         print*,'  and the ',N_states_diag,'states requested'
         print*,'  We did not find only states with S^2 values close to ',expected_s2
         print*,'  We will then set the first N_states eigenvectors of the H matrix'
         print*,'  as the CI_eigenvectors'
         print*,'  You should consider more states and maybe ask for s2_eig to be .True. or just enlarge the CI space'
         print*,''
         do j=1,min(N_states_diag,N_det)
           do i=1,N_det
             CI_eigenvectors(i,j) = eigenvectors(i,j)
           enddo
           CI_electronic_energy_real(j) = eigenvalues(j)
           CI_s2_real(j) = s2_eigvalues(j)
         enddo
       endif
       deallocate(index_good_state_array,good_state_array)
       deallocate(s2_eigvalues)
     else
       call lapack_diag(eigenvalues,eigenvectors,                    &
           H_matrix_all_dets,size(H_matrix_all_dets,1),N_det)
       CI_electronic_energy_real(:) = 0.d0
       call u_0_S2_u_0(CI_s2_real,eigenvectors,N_det,psi_det,N_int,&
          min(N_det,N_states_diag),size(eigenvectors,1))
       ! Select the "N_states_diag" states of lowest energy
       do j=1,min(N_det,N_states_diag)
         do i=1,N_det
           CI_eigenvectors(i,j) = eigenvectors(i,j)
         enddo
         CI_electronic_energy_real(j) = eigenvalues(j)
       enddo
     endif
     do k=1,N_states_diag
       CI_electronic_energy_real(k) = 0.d0
       do j=1,N_det
         do i=1,N_det
           CI_electronic_energy_real(k) +=                                &
               CI_eigenvectors(i,k) * CI_eigenvectors(j,k) *         &
               H_matrix_all_dets(i,j)
         enddo
       enddo
     enddo
     deallocate(eigenvectors,eigenvalues)
   endif

END_PROVIDER

 BEGIN_PROVIDER [ double precision, CI_electronic_energy_complex, (N_states_diag) ]
&BEGIN_PROVIDER [ complex*16, CI_eigenvectors_complex, (N_det,N_states_diag) ]
&BEGIN_PROVIDER [ double precision, CI_s2_complex, (N_states_diag) ]
   BEGIN_DOC
   ! Eigenvectors/values of the |CI| matrix
   END_DOC
   implicit none
   double precision               :: ovrlp
   complex*16                     :: u_dot_v_complex
   integer                        :: i_good_state
   integer, allocatable           :: index_good_state_array(:)
   logical, allocatable           :: good_state_array(:)
   double precision, allocatable  :: s2_values_tmp(:)
   integer                        :: i_other_state
   double precision, allocatable  :: eigenvalues(:)
   complex*16, allocatable        :: eigenvectors(:,:), H_prime(:,:)
   integer                        :: i_state
   double precision               :: e_0
   integer                        :: i,j,k
   double precision, allocatable  :: s2_eigvalues(:)
   double precision, allocatable  :: e_array(:)
   integer, allocatable           :: iorder(:)
   logical                        :: converged

   PROVIDE threshold_davidson nthreads_davidson
   ! Guess values for the "N_states" states of the |CI| eigenvectors
   do j=1,min(N_states,N_det)
     do i=1,N_det
       ci_eigenvectors_complex(i,j) = psi_coef_complex(i,j)
     enddo
   enddo

   do j=min(N_states,N_det)+1,N_states_diag
     do i=1,N_det
       ci_eigenvectors_complex(i,j) = (0.d0,0.d0)
     enddo
   enddo

   if (diag_algorithm == "Davidson") then

     call davidson_diag_hs2_complex(psi_det,ci_eigenvectors_complex, ci_s2_complex, &
         size(ci_eigenvectors_complex,1),ci_electronic_energy_complex,               &
         N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,0,converged)

     integer :: N_states_diag_save
     N_states_diag_save = N_states_diag
     do while (.not.converged)
        double precision, allocatable :: ci_electronic_energy_tmp (:)
        complex*16, allocatable :: ci_eigenvectors_tmp (:,:)
        double precision, allocatable :: ci_s2_tmp (:)

        N_states_diag  *= 2
        TOUCH N_states_diag

        allocate (ci_electronic_energy_tmp (N_states_diag) )
        allocate (ci_eigenvectors_tmp (N_det,N_states_diag) )
        allocate (ci_s2_tmp (N_states_diag) )

        ci_electronic_energy_tmp(1:N_states_diag_save) = ci_electronic_energy_complex(1:N_states_diag_save)
        ci_eigenvectors_tmp(1:N_det,1:N_states_diag_save) = ci_eigenvectors_complex(1:N_det,1:N_states_diag_save)
        ci_s2_tmp(1:N_states_diag_save) = ci_s2_complex(1:N_states_diag_save)

        call davidson_diag_hs2_complex(psi_det,ci_eigenvectors_tmp, ci_s2_tmp, &
            size(ci_eigenvectors_tmp,1),ci_electronic_energy_tmp,               &
            N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,0,converged)

        ci_electronic_energy_complex(1:N_states_diag_save) = ci_electronic_energy_tmp(1:N_states_diag_save)
        ci_eigenvectors_complex(1:N_det,1:N_states_diag_save) = ci_eigenvectors_tmp(1:N_det,1:N_states_diag_save)
        ci_s2_complex(1:N_states_diag_save) = ci_s2_tmp(1:N_states_diag_save)

        deallocate (ci_electronic_energy_tmp)
        deallocate (ci_eigenvectors_tmp)
        deallocate (ci_s2_tmp)
     enddo
     if (N_states_diag > N_states_diag_save) then
       N_states_diag = N_states_diag_save
       TOUCH N_states_diag
     endif

   else if (diag_algorithm == "Lapack") then

     print *,  'Diagonalization of H using Lapack'
     allocate (eigenvectors(size(h_matrix_all_dets_complex,1),N_det))
     allocate (eigenvalues(N_det))
     if (s2_eig) then
       double precision, parameter :: alpha = 0.1d0
       allocate (H_prime(N_det,N_det) )
       H_prime(1:N_det,1:N_det) = h_matrix_all_dets_complex(1:N_det,1:N_det) +  &
         alpha * s2_matrix_all_dets(1:N_det,1:N_det)
       do j=1,N_det
         H_prime(j,j) = H_prime(j,j) - alpha*expected_s2
       enddo
       call lapack_diag_complex(eigenvalues,eigenvectors,H_prime,size(H_prime,1),N_det)
       ci_electronic_energy_complex(:) = 0.d0
       i_state = 0
       allocate (s2_eigvalues(N_det))
       allocate(index_good_state_array(N_det),good_state_array(N_det))
       good_state_array = .False.
       call u_0_s2_u_0_complex(s2_eigvalues,eigenvectors,N_det,psi_det,N_int,&
         N_det,size(eigenvectors,1))
       if (only_expected_s2) then
          do j=1,N_det
            ! Select at least n_states states with S^2 values closed to "expected_s2"
            if(dabs(s2_eigvalues(j)-expected_s2).le.0.5d0)then
              i_state +=1
              index_good_state_array(i_state) = j
              good_state_array(j) = .True.
            endif
            if(i_state.eq.N_states) then
              exit
            endif
          enddo
       else
          do j=1,N_det
             index_good_state_array(j) = j
             good_state_array(j) = .True.
          enddo
       endif
       if(i_state .ne.0)then
         ! Fill the first "i_state" states that have a correct S^2 value
         do j = 1, i_state
           do i=1,N_det
             ci_eigenvectors_complex(i,j) = eigenvectors(i,index_good_state_array(j))
           enddo
           ci_electronic_energy_complex(j) = eigenvalues(index_good_state_array(j))
           ci_s2_complex(j) = s2_eigvalues(index_good_state_array(j))
         enddo
         i_other_state = 0
         do j = 1, N_det
           if(good_state_array(j))cycle
           i_other_state +=1
           if(i_state+i_other_state.gt.n_states_diag)then
             exit
           endif
           do i=1,N_det
             ci_eigenvectors_complex(i,i_state+i_other_state) = eigenvectors(i,j)
           enddo
           ci_electronic_energy_complex(i_state+i_other_state) = eigenvalues(j)
           ci_s2_complex(i_state+i_other_state) = s2_eigvalues(i_state+i_other_state)
         enddo

       else
         print*,''
         print*,'!!!!!!!!   WARNING  !!!!!!!!!'
         print*,'  Within the ',N_det,'determinants selected'
         print*,'  and the ',N_states_diag,'states requested'
         print*,'  We did not find any state with S^2 values close to ',expected_s2
         print*,'  We will then set the first N_states eigenvectors of the H matrix'
         print*,'  as the ci_eigenvectors_complex'
         print*,'  You should consider more states and maybe ask for s2_eig to be .True. or just enlarge the CI space'
         print*,''
         do j=1,min(N_states_diag,N_det)
           do i=1,N_det
             ci_eigenvectors_complex(i,j) = eigenvectors(i,j)
           enddo
           ci_electronic_energy_complex(j) = eigenvalues(j)
           ci_s2_complex(j) = s2_eigvalues(j)
         enddo
       endif
       deallocate(index_good_state_array,good_state_array)
       deallocate(s2_eigvalues)
     else
       call lapack_diag_complex(eigenvalues,eigenvectors,                      &
           H_matrix_all_dets_complex,size(H_matrix_all_dets_complex,1),N_det)
       ci_electronic_energy_complex(:) = 0.d0
       call u_0_S2_u_0_complex(ci_s2_complex,eigenvectors,N_det,psi_det,N_int,&
          min(N_det,N_states_diag),size(eigenvectors,1))
       ! Select the "N_states_diag" states of lowest energy
       do j=1,min(N_det,N_states_diag)
         do i=1,N_det
           ci_eigenvectors_complex(i,j) = eigenvectors(i,j)
         enddo
         ci_electronic_energy_complex(j) = eigenvalues(j)
       enddo
     endif
     do k=1,N_states_diag
       ci_electronic_energy_complex(k) = 0.d0
       do j=1,N_det
         do i=1,N_det
           !todo: accumulate imag parts to test? (should sum to zero)
           ci_electronic_energy_complex(k) +=                                &
               dble(dconjg(ci_eigenvectors_complex(i,k)) * ci_eigenvectors_complex(j,k) *  &
               H_matrix_all_dets_complex(i,j))
         enddo
       enddo
     enddo
     deallocate(eigenvectors,eigenvalues)
   endif

END_PROVIDER

subroutine diagonalize_ci
  implicit none
  if (is_complex) then
    call diagonalize_ci_complex
  else
    call diagonalize_ci_real
  endif
end

subroutine diagonalize_CI_complex
  implicit none
  BEGIN_DOC
!  Replace the coefficients of the |CI| states by the coefficients of the
!  eigenstates of the |CI| matrix.
  END_DOC
  integer :: i,j
  do j=1,N_states
    do i=1,N_det
      psi_coef_complex(i,j) = ci_eigenvectors_complex(i,j)
    enddo
  enddo
  psi_energy(1:N_states) = CI_electronic_energy(1:N_states)
  psi_s2(1:N_states) = CI_s2(1:N_states)
  !todo: touch ci_{s2,electronic_energy}?
  SOFT_TOUCH psi_coef_complex CI_electronic_energy_complex ci_energy CI_eigenvectors_complex CI_s2_complex psi_energy psi_s2
end

subroutine diagonalize_CI_real
  implicit none
  BEGIN_DOC
!  Replace the coefficients of the |CI| states by the coefficients of the
!  eigenstates of the |CI| matrix.
  END_DOC
  integer                        :: i,j
  do j=1,N_states
    do i=1,N_det
      psi_coef(i,j) = CI_eigenvectors(i,j)
    enddo
  enddo
  psi_energy(1:N_states) = CI_electronic_energy(1:N_states)
  psi_s2(1:N_states) = CI_s2(1:N_states)

  !todo: touch ci_{s2,electronic_energy}?
  SOFT_TOUCH psi_coef CI_electronic_energy_real ci_energy CI_eigenvectors CI_s2_real psi_energy psi_s2
end
