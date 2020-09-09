!BEGIN_PROVIDER [ integer, ao_num_per_kpt ]
! implicit none
! BEGIN_DOC
! ! number of aos per kpt.
! END_DOC
! ao_num_per_kpt = ao_num/kpt_num
!END_PROVIDER

subroutine get_kpt_idx_ao(idx_full,k,i)
  implicit none
  BEGIN_DOC
  ! idx_full is ao index in full range (up to ao_num)
  ! k is index of the k-point for this ao
  ! i is index of this ao within k-point k
  ! this assumes that all kpts have the same number of aos
  END_DOC

  integer, intent(in) :: idx_full
  integer, intent(out) :: i,k
  i = mod(idx_full-1,ao_num_per_kpt)+1
  k = (idx_full-1)/ao_num_per_kpt+1
  ASSERT (k <= kpt_num)
end
