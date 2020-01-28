
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
    mo_overlap(i,j) = (0.d0,0.d0)
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

