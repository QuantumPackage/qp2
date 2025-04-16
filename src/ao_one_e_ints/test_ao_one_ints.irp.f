
BEGIN_PROVIDER [double precision, ao_ortho_canonical_overlap, (ao_ortho_canonical_num,ao_ortho_canonical_num)]
  implicit none
  BEGIN_DOC
! overlap matrix of the ao_ortho_canonical.
! Expected to be the Identity
  END_DOC
  integer                        :: i,j,k,l
  double precision               :: c
  do j=1, ao_ortho_canonical_num
    do i=1, ao_ortho_canonical_num
      ao_ortho_canonical_overlap(i,j) = 0.d0
    enddo
  enddo
  do j=1, ao_ortho_canonical_num
    do k=1, ao_num
      c = 0.d0
      do l=1, ao_num
        c +=  ao_ortho_canonical_coef(l,j) * ao_overlap(l,k)
      enddo
      do i=1, ao_ortho_canonical_num
        ao_ortho_canonical_overlap(i,j) += ao_ortho_canonical_coef(k,i) * c
      enddo
    enddo
  enddo
END_PROVIDER

