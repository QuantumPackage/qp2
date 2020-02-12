subroutine orthonormalize_mos
  implicit none
  integer :: m,p,s
  if (is_complex) then
    m = size(mo_coef_complex,1)
    p = size(mo_overlap_complex,1)
    call ortho_lowdin_complex(mo_overlap_complex,p,mo_num,mo_coef_complex,m,ao_num)
    mo_label = 'Orthonormalized'
    SOFT_TOUCH mo_coef_complex mo_label
    !TODO: should we do anything with the separate real/imag parts of mo_coef_complex?
  else
    m = size(mo_coef,1)
    p = size(mo_overlap,1)
    call ortho_lowdin(mo_overlap,p,mo_num,mo_coef,m,ao_num)
    mo_label = 'Orthonormalized'
    SOFT_TOUCH mo_coef mo_label
  endif
end


