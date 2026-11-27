use sizet

subroutine orthonormalize_mos
 use sizet
  implicit none
  BEGIN_DOC
  ! Replace mo_coef in place by the Lowdin-orthogonalized
  ! coefficients: $mo\_coef \leftarrow mo\_coef.S^{-1/2}$, $S$
  ! being mo_overlap, so that the MOs are orthonormal with
  ! respect to the AO overlap.
  ! Linear dependencies of $S$ are dropped, controlled by
  ! lin_dep_cutoff (see ortho_lowdin).
  END_DOC
  integer(size_t) :: m,p
  integer :: s,i
  m = size(mo_coef,1)
  p = size(mo_overlap,1)
    call ortho_lowdin(mo_overlap, int(p,size_t), int(mo_num,size_t), mo_coef, &
                    int(m,size_t), int(ao_num,size_t), lin_dep_cutoff)
  SOFT_TOUCH mo_coef
end


