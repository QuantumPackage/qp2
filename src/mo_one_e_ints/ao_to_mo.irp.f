subroutine mo_to_ao(A_mo,LDA_mo,A_ao,LDA_ao)
  implicit none
  BEGIN_DOC
  ! Transform A from the MO basis to the AO basis
  !
  ! $(S.C).A_{mo}.(S.C)^\dagger$
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_mo(LDA_mo,mo_num)
  double precision, intent(out)  :: A_ao(LDA_ao,ao_num)
  double precision, allocatable  :: T(:,:)

  allocate ( T(mo_num,ao_num) )

  call dgemm('N','T', mo_num, ao_num, mo_num,                &
      1.d0, A_mo,size(A_mo,1),                                       &
      S_mo_coef, size(S_mo_coef,1),                                  &
      0.d0, T, size(T,1))

  call dgemm('N','N', ao_num, ao_num, mo_num,                    &
      1.d0, S_mo_coef, size(S_mo_coef,1),                            &
      T, size(T,1),                                                  &
      0.d0, A_ao, size(A_ao,1))

  deallocate(T)
end

subroutine mo_to_ao_no_overlap(A_mo,LDA_mo,A_ao,LDA_ao)
  implicit none
  BEGIN_DOC
  ! $C.A_{mo}.C^\dagger$
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_mo(LDA_mo,mo_num)
  double precision, intent(out)  :: A_ao(LDA_ao,ao_num)
  double precision, allocatable  :: T(:,:)

  allocate ( T(mo_num,ao_num) )

  call dgemm('N','T', mo_num, ao_num, mo_num,                &
      1.d0, A_mo,size(A_mo,1),                                       &
      mo_coef, size(mo_coef,1),                                  &
      0.d0, T, size(T,1))

  call dgemm('N','N', ao_num, ao_num, mo_num,                    &
      1.d0, mo_coef, size(mo_coef,1),                            &
      T, size(T,1),                                                  &
      0.d0, A_ao, size(A_ao,1))

  deallocate(T)
end

BEGIN_PROVIDER [ double precision, S_mo_coef, (ao_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Product S.C where S is the overlap matrix in the AO basis and C the mo_coef matrix.
 END_DOC

 call dgemm('N','N', ao_num, mo_num, ao_num,                   &
     1.d0, ao_overlap,size(ao_overlap,1),      &
     mo_coef, size(mo_coef,1),                                     &
     0.d0, S_mo_coef, size(S_mo_coef,1))

END_PROVIDER


