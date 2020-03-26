BEGIN_PROVIDER [integer, kpt_pair_num]
  implicit none
  kpt_pair_num = shiftr(kpt_num*kpt_num+kpt_num,1)
END_PROVIDER

BEGIN_PROVIDER [integer, kconserv, (kpt_num,kpt_num,kpt_num)]
  implicit none
  BEGIN_DOC
  ! Information about k-point symmetry
  !
  ! for k-points I,J,K: kconserv(I,J,K) gives L such that
  ! k_I + k_J = k_K + k_L
  ! two-electron integrals of the form <ij|kx>
  ! (where i,j,k have momentum k_I, k_J, k_K)
  ! will only be nonzero if x has momentum k_L (as described above)
  !
  END_DOC
  integer                        :: i,j,k,l

  if (read_kconserv) then
    call ezfio_get_nuclei_kconserv(kconserv)
    print *,  'kconserv read from disk'
  else
    print*,'kconserv must be provided'
    stop -1
  endif
  if (write_kconserv) then
    call ezfio_set_nuclei_kconserv(kconserv)
    print *,  'kconserv written to disk'
  endif
END_PROVIDER

subroutine double_allowed_kpts(kh1,kh2,kp1,kp2,is_allowed)
  implicit none
  integer, intent(in) :: kh1,kh2,kp1,kp2
  logical, intent(out) :: is_allowed

  is_allowed = (kconserv(kh1,kh2,kp1) == kp2)
end subroutine

