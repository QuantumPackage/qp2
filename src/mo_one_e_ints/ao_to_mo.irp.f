use sizet

subroutine mo_to_ao(A_mo,LDA_mo,A_ao,LDA_ao)
 use sizet
  implicit none
  BEGIN_DOC
  ! Transform A from the MO basis to the AO basis
  !
  ! $(S.C).A_{mo}.(S.C)^\dagger$
  END_DOC
  integer(size_t), intent(in)    :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_mo(LDA_mo,mo_num)
  double precision, intent(out)  :: A_ao(LDA_ao,ao_num)
  double precision, allocatable  :: T(:,:)

  allocate ( T(mo_num,ao_num) )

  call dgemm('N','T', int(mo_num,4), int(ao_num,4), int(mo_num,4), &
      1.d0, A_mo,int(size(A_mo,1),4),                               &
      S_mo_coef, int(size(S_mo_coef,1),4),                          &
      0.d0, T, int(size(T,1),4))

  call dgemm('N','N', int(ao_num,4), int(ao_num,4), int(mo_num,4),  &
      1.d0, S_mo_coef, int(size(S_mo_coef,1),4),                    &
      T, int(size(T,1),4),                                          &
      0.d0, A_ao, int(size(A_ao,1),4))

  deallocate(T)
end

subroutine mo_to_ao_no_overlap(A_mo,LDA_mo,A_ao,LDA_ao)
 use sizet
  implicit none
  BEGIN_DOC
  ! $C.A_{mo}.C^\dagger$
  END_DOC
  integer(size_t), intent(in)    :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_mo(LDA_mo,mo_num)
  double precision, intent(out)  :: A_ao(LDA_ao,ao_num)
  double precision, allocatable  :: T(:,:)

  allocate ( T(mo_num,ao_num) )

  call dgemm('N','T', int(mo_num,4), int(ao_num,4), int(mo_num,4), &
      1.d0, A_mo,int(size(A_mo,1),4),                               &
      mo_coef, int(size(mo_coef,1),4),                              &
      0.d0, T, int(size(T,1),4))

  call dgemm('N','N', int(ao_num,4), int(ao_num,4), int(mo_num,4),  &
      1.d0, mo_coef, int(size(mo_coef,1),4),                        &
      T, int(size(T,1),4),                                          &
      0.d0, A_ao, int(size(A_ao,1),4))

  deallocate(T)
end

BEGIN_PROVIDER [ double precision, S_mo_coef, (ao_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Product S.C where S is the overlap matrix in the AO basis and C the mo_coef matrix.
 END_DOC

 call dgemm('N','N', int(ao_num,4), int(mo_num,4), int(ao_num,4), &
     1.d0, ao_overlap,int(size(ao_overlap,1),4),  &
     mo_coef, int(size(mo_coef,1),4),                &
     0.d0, S_mo_coef, int(size(S_mo_coef,1),4))

END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_one_e_integrals_from_mo, (ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! Integrals of the one e hamiltonian obtained from the integrals on the MO basis 
!
! WARNING : this is equal to ao_one_e_integrals only if the AO and MO basis have the same number of functions
 END_DOC
 call mo_to_ao(mo_one_e_integrals, int(mo_num, size_t), &
     ao_one_e_integrals_from_mo, int(ao_num, size_t))
END_PROVIDER 



! ---


subroutine mo_to_ao_sphe(A_mo,LDA_mo,A_ao,LDA_ao)
 use sizet
  implicit none
  BEGIN_DOC
  ! Transform A from the MO basis to the AO basis
  !
  ! $(S.C).A_{mo}.(S.C)^\dagger$
  END_DOC
  integer(size_t), intent(in)    :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_mo(LDA_mo,mo_num)
  double precision, intent(out)  :: A_ao(LDA_ao,ao_sphe_num)
  double precision, allocatable  :: T(:,:)

  allocate ( T(mo_num,ao_sphe_num) )

  call dgemm('N','T', int(mo_num,4), int(ao_sphe_num,4), int(mo_num,4), &
      1.d0, A_mo,int(size(A_mo,1),4),                               &
      S_mo_sphe_coef, int(size(S_mo_sphe_coef,1),4),                &
      0.d0, T, int(size(T,1),4))

  call dgemm('N','N', int(ao_sphe_num,4), int(ao_sphe_num,4), int(mo_num,4), &
      1.d0, S_mo_sphe_coef, int(size(S_mo_sphe_coef,1),4),          &
      T, int(size(T,1),4),                                          &
      0.d0, A_ao, int(size(A_ao,1),4))

  deallocate(T)
end

subroutine mo_to_ao_sphe_no_overlap(A_mo,LDA_mo,A_ao,LDA_ao)
 use sizet
  implicit none
  BEGIN_DOC
  ! $C.A_{mo}.C^\dagger$
  END_DOC
  integer(size_t), intent(in)    :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_mo(LDA_mo,mo_num)
  double precision, intent(out)  :: A_ao(LDA_ao,ao_sphe_num)
  double precision, allocatable  :: T(:,:)

  allocate ( T(mo_num,ao_sphe_num) )

  call dgemm('N','T', int(mo_num,4), int(ao_sphe_num,4), int(mo_num,4), &
      1.d0, A_mo,int(size(A_mo,1),4),                               &
      mo_sphe_coef, int(size(mo_sphe_coef,1),4),                    &
      0.d0, T, int(size(T,1),4))

  call dgemm('N','N', int(ao_sphe_num,4), int(ao_sphe_num,4), int(mo_num,4), &
      1.d0, mo_sphe_coef, int(size(mo_sphe_coef,1),4),              &
      T, int(size(T,1),4),                                          &
      0.d0, A_ao, int(size(A_ao,1),4))

  deallocate(T)
end

BEGIN_PROVIDER [ double precision, S_mo_sphe_coef, (ao_sphe_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Product S.C where S is the overlap matrix in the AO basis and C the mo_sphe_coef matrix.
 END_DOC

 call dgemm('N','N', int(ao_sphe_num,4), int(mo_num,4), int(ao_sphe_num,4), &
     1.d0, ao_overlap,int(size(ao_overlap,1),4),  &
     mo_sphe_coef, int(size(mo_coef,1),4),        &
     0.d0, S_mo_coef, int(size(S_mo_coef,1),4))

END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_sphe_one_e_integrals_from_mo, (ao_sphe_num, ao_sphe_num)]
 implicit none
 BEGIN_DOC
! Integrals of the one e hamiltonian obtained from the integrals on the MO basis 
!
! WARNING : this is equal to ao_one_e_integrals only if the AO and MO basis have the same number of functions
 END_DOC
 call mo_to_ao_sphe(mo_one_e_integrals, int(mo_num, size_t), &
     ao_one_e_integrals_from_mo, int(ao_sphe_num, size_t))
END_PROVIDER 



