subroutine orthonormalize_mos
  implicit none
  BEGIN_DOC
  ! Replace mo_coef in place by the Lowdin-orthogonalized
  ! coefficients: $mo\_coef \leftarrow mo\_coef.S^{-1/2}$, $S$
  ! being mo_overlap, so that the MOs are orthonormal with
  ! respect to the AO overlap.
  ! Linear dependencies of $S$ are dropped, controlled by
  ! lin_dep_cutoff (see ortho_lowdin).
  END_DOC
  integer :: m,p,s,i
  m = size(mo_coef,1)
  p = size(mo_overlap,1)
  call ortho_lowdin(mo_overlap,p,mo_num,mo_coef,m,ao_num,lin_dep_cutoff)
  SOFT_TOUCH mo_coef
end


