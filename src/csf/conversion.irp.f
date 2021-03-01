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
  integer*8                      :: Isomo, Idomo
  integer(bit_kind)              :: Ialpha(N_int) ,Ibeta(N_int)
  integer                        :: rows, cols, i, j, k
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
    
    s = 0
    do k=1,N_int
      if (psi_configuration(k,1,i) == 0_bit_kind) cycle
      s = s + popcnt(psi_configuration(k,1,i))
    enddo
    bfIcfg = max(1,nint((binom(s,(s+1)/2)-binom(s,((s+1)/2)+1))))
    
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
  
end subroutine convertWFfromDETtoCSF


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
  integer                        :: s, bfIcfg
  integer                        :: countcsf
  integer(bit_kind)              :: Ialpha(N_int), Ibeta(N_int)
  integer                        :: rows, cols, i, j, k
  integer                        :: startdet, enddet
  integer                        :: ndetI
  integer                        :: getNSOMO
  double precision,allocatable   :: tempBuffer(:,:)
  double precision,allocatable   :: tempCoeff(:,:)
  double precision               :: phasedet
  integer                        :: idx
 
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
    bfIcfg = max(1,nint((binom(s,(s+1)/2)-binom(s,((s+1)/2)+1))))
    
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

end subroutine convertCSFtoDET
