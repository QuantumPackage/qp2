subroutine orthonormalize_mos
  implicit none
  integer :: m,p,s,i
  m = size(mo_coef,1)
  p = size(mo_overlap,1)
  do i=1,4
    call ortho_lowdin(mo_overlap,p,mo_num,mo_coef,m,ao_num,lin_dep_cutoff)
    call nullify_small_elements(ao_num,mo_num,mo_coef,size(mo_coef,1),1.d-10)
  enddo
  if (restore_symm) then
     call restore_symmetry(ao_num, mo_num, mo_coef, size(mo_coef,1), 1.d-10)
  endif
  SOFT_TOUCH mo_coef
end


