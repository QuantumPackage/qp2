BEGIN_PROVIDER [ double precision, psi_csf_coef, (N_csf, N_states) ]
 implicit none
 BEGIN_DOC
 ! Wafe function in CSF basis
 END_DOC

 double precision, allocatable :: buffer(:,:)
 allocate ( buffer(N_det, N_states) )
 buffer(1:N_det, 1:N_states) = psi_coef(1:N_det, 1:N_states)
 call convertWFfromDETtoCSF(N_states, buffer, psi_csf_coef)
END_PROVIDER

subroutine convertWFfromDETtoCSF(N_st,psi_coef_det_in, psi_coef_cfg_out)
  use cfunctions
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Documentation for DetToCSFTransformationMatrix
  ! Provides the matrix of transformatons for the
  ! conversion between determinant to CSF basis (in BFs)
  END_DOC
  integer, intent(in)            :: N_st
  double precision, intent(in)   :: psi_coef_det_in(N_det,N_st)
  double precision, intent(out)  :: psi_coef_cfg_out(n_CSF,N_st)
  integer*8                      :: Isomo, Idomo, mask
  integer(bit_kind)              :: Ialpha(N_int) ,Ibeta(N_int)
  integer                        :: rows, cols, i, j, k, salpha
  integer                        :: startdet, enddet
  integer                        :: ndetI
  integer                        :: getNSOMO
  double precision,allocatable   :: tempBuffer(:,:)
  double precision,allocatable   :: tempCoeff(:,:)
  double precision               :: phasedet
  integer                        :: idx

  ! initialization
  psi_coef_cfg_out(:,1) = 0.d0

  integer s, bfIcfg
  integer countcsf
  integer MS
  MS = elec_alpha_num-elec_beta_num
  countcsf = 0
  phasedet = 1.0d0
  do i = 1,N_configuration
    startdet = psi_configuration_to_psi_det(1,i)
    enddet = psi_configuration_to_psi_det(2,i)
    ndetI = enddet-startdet+1

    allocate(tempCoeff(ndetI,N_st))
    do j = startdet, enddet
      idx = psi_configuration_to_psi_det_data(j)
      Ialpha(:) = psi_det(:,1,idx)
      Ibeta(:)  = psi_det(:,2,idx)
      call get_phase_qp_to_cfg(Ialpha, Ibeta, phasedet)
      do k=1,N_st
        tempCoeff(j-startdet+1,k) = psi_coef_det_in(idx, k)*phasedet
      enddo
    enddo

    s = 0 ! s == total number of SOMOs
    do k=1,N_int
      if (psi_configuration(k,1,i) == 0_bit_kind) cycle
      s = s + popcnt(psi_configuration(k,1,i))
    enddo

    if(iand(s,1) .EQ. 0) then
      salpha = (s + MS)/2
      bfIcfg = max(1,nint((binom(s,salpha)-binom(s,salpha+1))))
    else
      salpha = (s + MS)/2
      bfIcfg = max(1,nint((binom(s,salpha)-binom(s,salpha+1))))
    endif

    ! perhaps blocking with CFGs of same seniority
    ! can be more efficient
    allocate(tempBuffer(bfIcfg,ndetI))
    tempBuffer = DetToCSFTransformationMatrix(s,:bfIcfg,:ndetI)

    call dgemm('N','N', bfIcfg, N_st, ndetI, 1.d0, tempBuffer, size(tempBuffer,1),&
        tempCoeff, size(tempCoeff,1), 0.d0, psi_coef_cfg_out(countcsf+1,1),&
        size(psi_coef_cfg_out,1))

    deallocate(tempCoeff)
    deallocate(tempBuffer)
    countcsf += bfIcfg
  enddo

end


subroutine convertWFfromCSFtoDET(N_st,psi_coef_cfg_in, psi_coef_det)
  implicit none
  BEGIN_DOC
  ! Documentation for convertCSFtoDET
  ! This function converts the wavefunction
  ! in CFG basis to DET basis using the
  ! transformation matrix provided before.
  END_DOC
  integer, intent(in)            :: N_st
  double precision,intent(in)    :: psi_coef_cfg_in(n_CSF,N_st)
  double precision,intent(out)   :: psi_coef_det(N_det,N_st)
  double precision               :: tmp_psi_coef_det(maxDetDimPerBF,N_st)
  integer                        :: s, bfIcfg, salpha
  integer                        :: countcsf
  integer(bit_kind)              :: Ialpha(N_int), Ibeta(N_int)
  integer                        :: rows, cols, i, j, k
  integer                        :: startdet, enddet
  integer                        :: ndetI
  integer                        :: getNSOMO
  double precision,allocatable   :: tempBuffer(:,:)
  double precision,allocatable   :: tempCoeff (:,:)
  double precision               :: phasedet
  integer                        :: idx
  integer MS
  MS = elec_alpha_num-elec_beta_num
  !print *,"size=",size(tmp_psi_coef_det,1)," ",size(tmp_psi_coef_det,2)

  countcsf = 0

  do i = 1,N_configuration
    startdet = psi_configuration_to_psi_det(1,i)
    enddet = psi_configuration_to_psi_det(2,i)
    ndetI = enddet-startdet+1

    s = 0
    do k=1,N_int
      if (psi_configuration(k,1,i) == 0_bit_kind) cycle
      s = s + popcnt(psi_configuration(k,1,i))
    enddo
    salpha = (s + MS)/2
    bfIcfg = max(1,nint((binom(s,salpha)-binom(s,salpha+1))))

    allocate(tempCoeff(bfIcfg,N_st))

    do k=1,N_st
      do j = 1,bfIcfg
        tempCoeff(j,k) = psi_coef_cfg_in(countcsf+j,k)
      enddo
    enddo

    countcsf += bfIcfg
    ! perhaps blocking with CFGs of same seniority
    ! can be more efficient
    allocate(tempBuffer(bfIcfg,ndetI))
    tempBuffer = DetToCSFTransformationMatrix(s,:bfIcfg,:ndetI)

    call dgemm('T','N', ndetI, N_st, bfIcfg, 1.d0, tempBuffer, size(tempBuffer,1),&
        tempCoeff, size(tempCoeff,1), 0.d0, tmp_psi_coef_det,        &
        size(tmp_psi_coef_det,1))

    do j=startdet,enddet
      idx = psi_configuration_to_psi_det_data(j)
      Ialpha(:) = psi_det(:,1,idx)
      Ibeta(:)  = psi_det(:,2,idx)
      call get_phase_qp_to_cfg(Ialpha, Ibeta, phasedet)
      do k=1,N_st
        psi_coef_det(idx,k) = tmp_psi_coef_det(j-startdet+1,k) * phasedet
      enddo
    enddo

    deallocate(tempCoeff)
    deallocate(tempBuffer)
  enddo

end







