subroutine orthonormalize_mos
  implicit none
  integer :: m,p,s,k
  if (is_complex) then
    do k=1,kpt_num
      m = size(mo_coef_kpts,1)
      p = size(mo_overlap_kpts,1)
      call ortho_lowdin_complex(mo_overlap_kpts(1,1,k),p,mo_num_per_kpt,mo_coef_kpts(1,1,k),m,ao_num_per_kpt)
    enddo
    mo_label = 'Orthonormalized'
    SOFT_TOUCH mo_coef_kpts mo_label
  else
    m = size(mo_coef,1)
    p = size(mo_overlap,1)
    call ortho_lowdin(mo_overlap,p,mo_num,mo_coef,m,ao_num)
    mo_label = 'Orthonormalized'
    SOFT_TOUCH mo_coef mo_label
  endif
end


