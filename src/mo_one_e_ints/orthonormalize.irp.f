subroutine orthonormalize_mos
  implicit none
  integer :: m,p,s
  m = size(mo_coef,1)
  p = size(mo_overlap,1)
  call ortho_lowdin(mo_overlap,p,mo_num,mo_coef,m,ao_num,lin_dep_cutoff)
  mo_label = 'Orthonormalized'
  SOFT_TOUCH mo_coef mo_label
end


