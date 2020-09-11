
BEGIN_PROVIDER [ double precision, mo_overlap,(mo_num,mo_num) ]
  implicit none
  BEGIN_DOC
! Provider to check that the MOs are indeed orthonormal.
  END_DOC
  integer :: i,j,n,l
  double precision :: f
  integer :: lmax


  lmax = (ao_num/4) * 4
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
  !$OMP  PRIVATE(i,j,n,l) &
  !$OMP  SHARED(mo_overlap,mo_coef,ao_overlap, &
  !$OMP    mo_num,ao_num,lmax)
  do j=1,mo_num
   do i= 1,mo_num
    mo_overlap(i,j) = 0.d0
    do n = 1, lmax,4
     do l = 1, ao_num
      mo_overlap(i,j) = mo_overlap(i,j) + mo_coef(l,i) * &
           ( mo_coef(n  ,j) * ao_overlap(l,n  )  &
           + mo_coef(n+1,j) * ao_overlap(l,n+1)  &
           + mo_coef(n+2,j) * ao_overlap(l,n+2)  &
           + mo_coef(n+3,j) * ao_overlap(l,n+3)  )
     enddo
    enddo
    do n = lmax+1, ao_num
     do l = 1, ao_num
      mo_overlap(i,j) = mo_overlap(i,j) + mo_coef(n,j) * mo_coef(l,i) * ao_overlap(l,n)
     enddo
    enddo
   enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_overlap_complex,(mo_num,mo_num) ]
  implicit none
  BEGIN_DOC
! Provider to check that the MOs are indeed orthonormal.
  END_DOC
  integer :: i,j,n,l
  integer :: lmax


  lmax = (ao_num/4) * 4
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
  !$OMP  PRIVATE(i,j,n,l) &
  !$OMP  SHARED(mo_overlap_complex,mo_coef_complex,ao_overlap_complex, &
  !$OMP    mo_num,ao_num,lmax)
  do j=1,mo_num
   do i= 1,mo_num
    mo_overlap_complex(i,j) = (0.d0,0.d0)
    do n = 1, lmax,4
     do l = 1, ao_num
      mo_overlap_complex(i,j) = mo_overlap_complex(i,j) + dconjg(mo_coef_complex(l,i)) * &
           ( mo_coef_complex(n  ,j) * ao_overlap_complex(l,n  )  &
           + mo_coef_complex(n+1,j) * ao_overlap_complex(l,n+1)  &
           + mo_coef_complex(n+2,j) * ao_overlap_complex(l,n+2)  &
           + mo_coef_complex(n+3,j) * ao_overlap_complex(l,n+3)  )
     enddo
    enddo
    do n = lmax+1, ao_num
     do l = 1, ao_num
      mo_overlap_complex(i,j) = mo_overlap_complex(i,j) + mo_coef_complex(n,j) * dconjg(mo_coef_complex(l,i)) * ao_overlap_complex(l,n)
     enddo
    enddo
   enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_overlap_kpts,(mo_num_per_kpt,mo_num_per_kpt,kpt_num) ]
  implicit none
  BEGIN_DOC
! Provider to check that the MOs are indeed orthonormal.
  END_DOC
  integer :: i,j,n,l,k
  integer :: lmax

  print *, 'Providing MO overlap integrals'
  if (read_mo_integrals_overlap) then
    call ezfio_get_mo_one_e_ints_mo_integrals_overlap_kpts(mo_overlap_kpts)
    print *,  'MO overlap integrals read from disk'
  else
  print *, 'Providing MO overlap integrals from AO overlap integrals'
  !  call ao_to_mo_kpts(                                            &
  !      ao_kinetic_integrals_kpts,                                 &
  !      size(ao_kinetic_integrals_kpts,1),                         &
  !      mo_kinetic_integrals_kpts,                                 &
  !      size(mo_kinetic_integrals_kpts,1)                          &
  !      )
  !endif
  

  lmax = (ao_num_per_kpt/4) * 4
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
  !$OMP  PRIVATE(i,j,n,l,k) &
  !$OMP  SHARED(mo_overlap_kpts,mo_coef_kpts,ao_overlap_kpts, &
  !$OMP    mo_num_per_kpt,ao_num_per_kpt,lmax,kpt_num)
  do k=1,kpt_num
    do j=1,mo_num_per_kpt
      do i= 1,mo_num_per_kpt
        mo_overlap_kpts(i,j,k) = (0.d0,0.d0)
        do n = 1, lmax,4
          do l = 1, ao_num_per_kpt
            mo_overlap_kpts(i,j,k) = mo_overlap_kpts(i,j,k) + dconjg(mo_coef_kpts(l,i,k)) * &
                 ( mo_coef_kpts(n  ,j,k) * ao_overlap_kpts(l,n  ,k)  &
                 + mo_coef_kpts(n+1,j,k) * ao_overlap_kpts(l,n+1,k)  &
                 + mo_coef_kpts(n+2,j,k) * ao_overlap_kpts(l,n+2,k)  &
                 + mo_coef_kpts(n+3,j,k) * ao_overlap_kpts(l,n+3,k)  )
          enddo
        enddo
        do n = lmax+1, ao_num_per_kpt
          do l = 1, ao_num_per_kpt
            mo_overlap_kpts(i,j,k) = mo_overlap_kpts(i,j,k) + mo_coef_kpts(n,j,k) * &
              dconjg(mo_coef_kpts(l,i,k)) * ao_overlap_kpts(l,n,k)
          enddo
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO
  endif
END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_overlap_kpts_real, (mo_num_per_kpt, mo_num_per_kpt, kpt_num) ]
  implicit none
  BEGIN_DOC
  ! Overlap for complex MOs
  END_DOC
  integer :: i,j,k
  do k=1,kpt_num
    do j=1,mo_num_per_kpt
      do i=1,mo_num_per_kpt
        mo_overlap_kpts_real(i,j,k) = dble(mo_overlap_kpts(i,j,k))
      enddo
    enddo
  enddo
END_PROVIDER

