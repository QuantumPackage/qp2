subroutine mo_to_ao_complex(A_mo,LDA_mo,A_ao,LDA_ao)
  implicit none
  BEGIN_DOC
  ! Transform A from the MO basis to the AO basis
  !
  ! (S.C).A_mo.(S.C)t
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)         :: A_mo(LDA_mo,mo_num)
  complex*16, intent(out)        :: A_ao(LDA_ao,ao_num)
  complex*16, allocatable        :: T(:,:)
  
  allocate ( T(mo_num,ao_num) )
  
  call zgemm('N','C', mo_num, ao_num, mo_num,                &
      (1.d0,0.d0), A_mo,size(A_mo,1),                                &
      S_mo_coef_complex, size(S_mo_coef_complex,1),                                  &
      (0.d0,0.d0), T, size(T,1))
  
  call zgemm('N','N', ao_num, ao_num, mo_num,                    &
      (1.d0,0.d0), S_mo_coef_complex, size(S_mo_coef_complex,1),                     &
      T, size(T,1),                                                  &
      (0.d0,0.d0), A_ao, size(A_ao,1))
  
  deallocate(T)
end

subroutine mo_to_ao_no_overlap_complex(A_mo,LDA_mo,A_ao,LDA_ao)
  implicit none
  BEGIN_DOC
  ! Transform A from the MO basis to the S^-1 AO basis
  ! Useful for density matrix
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)         :: A_mo(LDA_mo,mo_num)
  complex*16, intent(out)        :: A_ao(LDA_ao,ao_num)
  complex*16, allocatable        :: T(:,:)
  
  allocate ( T(mo_num,ao_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  
  call zgemm('N','C', mo_num, ao_num, mo_num,                &
      (1.d0,0.d0), A_mo,size(A_mo,1),                                &
      mo_coef_complex, size(mo_coef_complex,1),                                      &
      (0.d0,0.d0), T, size(T,1))
  
  call zgemm('N','N', ao_num, ao_num, mo_num,                    &
      (1.d0,0.d0), mo_coef_complex,size(mo_coef_complex,1),                          &
      T, size(T,1),                                                  &
      (0.d0,0.d0), A_ao, size(A_ao,1))
  
  deallocate(T)
end

BEGIN_PROVIDER [ complex*16, S_mo_coef_complex, (ao_num, mo_num) ]
  implicit none
  BEGIN_DOC
  ! Product S.C where S is the overlap matrix in the AO basis and C the mo_coef matrix.
  END_DOC

  call zgemm('N','N',ao_num, mo_num, ao_num, (1.d0,0.d0), &
        ao_overlap_complex, size(ao_overlap_complex,1), &
        mo_coef_complex, size(mo_coef_complex,1), &
        (0.d0,0.d0), &
        S_mo_coef_complex, size(S_mo_coef_complex,1))

END_PROVIDER

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!

subroutine mo_to_ao_kpts(A_mo,LDA_mo,A_ao,LDA_ao)
  implicit none
  BEGIN_DOC
  ! Transform A from the MO basis to the AO basis
  !
  ! (S.C).A_mo.(S.C)t
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)         :: A_mo(LDA_mo,mo_num_per_kpt,kpt_num)
  complex*16, intent(out)        :: A_ao(LDA_ao,ao_num_per_kpt,kpt_num)
  complex*16, allocatable        :: T(:,:)
  
  allocate ( T(mo_num_per_kpt,ao_num_per_kpt) )
  integer :: k
  do k=1,kpt_num
    call zgemm('N','C', mo_num_per_kpt, ao_num_per_kpt, mo_num_per_kpt,                &
        (1.d0,0.d0), A_mo(:,:,k),size(A_mo,1),                                &
        S_mo_coef_kpts(:,:,k), size(S_mo_coef_kpts,1),                                  &
        (0.d0,0.d0), T, size(T,1))
    
    call zgemm('N','N', ao_num_per_kpt, ao_num_per_kpt, mo_num_per_kpt,                    &
        (1.d0,0.d0), S_mo_coef_kpts(:,:,k), size(S_mo_coef_kpts,1),                     &
        T, size(T,1),                                                  &
        (0.d0,0.d0), A_ao(:,:,k), size(A_ao,1))
  enddo
  deallocate(T)
end

subroutine mo_to_ao_no_overlap_kpts(A_mo,LDA_mo,A_ao,LDA_ao)
  implicit none
  BEGIN_DOC
  ! Transform A from the MO basis to the S^-1 AO basis
  ! Useful for density matrix
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)         :: A_mo(LDA_mo,mo_num_per_kpt,kpt_num)
  complex*16, intent(out)        :: A_ao(LDA_ao,ao_num_per_kpt,kpt_num)
  complex*16, allocatable        :: T(:,:)
  
  allocate ( T(mo_num_per_kpt,ao_num_per_kpt) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  integer :: k
  do k=1,kpt_num
    call zgemm('N','C', mo_num_per_kpt, ao_num_per_kpt, mo_num_per_kpt,                &
        (1.d0,0.d0), A_mo(:,:,k),size(A_mo,1),                                &
        mo_coef_kpts(:,:,k), size(mo_coef_kpts,1),                                      &
        (0.d0,0.d0), T, size(T,1))
    
    call zgemm('N','N', ao_num_per_kpt, ao_num_per_kpt, mo_num_per_kpt,                    &
        (1.d0,0.d0), mo_coef_kpts(:,:,k),size(mo_coef_kpts,1),                          &
        T, size(T,1),                                                  &
        (0.d0,0.d0), A_ao(:,:,k), size(A_ao,1))
  enddo
  deallocate(T)
end

BEGIN_PROVIDER [ complex*16, S_mo_coef_kpts, (ao_num_per_kpt, mo_num_per_kpt, kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Product S.C where S is the overlap matrix in the AO basis and C the mo_coef matrix.
  END_DOC

  integer :: k
  do k=1,kpt_num
    call zgemm('N','N',ao_num_per_kpt, mo_num_per_kpt, ao_num_per_kpt, (1.d0,0.d0), &
          ao_overlap_kpts(:,:,k), size(ao_overlap_kpts,1), &
          mo_coef_kpts(:,:,k), size(mo_coef_kpts,1), &
          (0.d0,0.d0), &
          S_mo_coef_kpts(:,:,k), size(S_mo_coef_kpts,1))
  enddo
END_PROVIDER

