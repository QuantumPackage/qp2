subroutine orthonormalize_mos
  implicit none
  integer :: m,p,s,k
  if (is_complex) then
    do k=1,kpt_num
      m = size(mo_coef_kpts,1)
      p = size(mo_overlap_kpts,1)
      call ortho_lowdin_complex(mo_overlap_kpts(1,1,k),p,mo_num_per_kpt,mo_coef_kpts(1,1,k),m,ao_num_per_kpt,lin_dep_cutoff)
    enddo
    mo_label = 'Orthonormalized'
    SOFT_TOUCH mo_coef_kpts mo_label
  else
    m = size(mo_coef,1)
    p = size(mo_overlap,1)
    call ortho_lowdin(mo_overlap,p,mo_num,mo_coef,m,ao_num,lin_dep_cutoff)
    mo_label = 'Orthonormalized'
    SOFT_TOUCH mo_coef mo_label
  endif
end


subroutine orthonormalize_mos_k_real
  implicit none
  integer :: m,p,s,k
  double precision, allocatable :: mo_coef_tmp(:,:)
  
  allocate(mo_coef_tmp(ao_num_per_kpt,mo_num_per_kpt))
  do k=1,kpt_num
    m = size(mo_coef_kpts,1)
    p = size(mo_overlap_kpts,1)
    mo_coef_tmp = dble(mo_coef_kpts(:,:,k))
    call ortho_lowdin(mo_overlap_kpts_real(1,1,k),p,mo_num_per_kpt,mo_coef_tmp,m,ao_num_per_kpt,lin_dep_cutoff)
    call zlacp2('X',ao_num_per_kpt,mo_num_per_kpt,mo_coef_tmp,size(mo_coef_tmp,1), &
           mo_coef_kpts(1,1,k),size(mo_coef_kpts,1))
  enddo
  deallocate(mo_coef_tmp)
  mo_label = 'Orthonormalized'
  SOFT_TOUCH mo_coef_kpts mo_label
end


