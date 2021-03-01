  BEGIN_PROVIDER [ integer, NSOMOMax]
 &BEGIN_PROVIDER [ integer, NCSFMax]
 &BEGIN_PROVIDER [ integer*8, NMO]
 &BEGIN_PROVIDER [ integer, NBFMax]
 &BEGIN_PROVIDER [ integer, dimBasisCSF]
 &BEGIN_PROVIDER [ integer, maxDetDimPerBF]
  implicit none
  BEGIN_DOC
  ! Documentation for NSOMOMax
  ! The maximum number of SOMOs for the current calculation.
  ! required for the calculation of prototype arrays.
  END_DOC
  NSOMOMax = min(elec_num, cfg_nsomo_max + 2)
  ! Note that here we need NSOMOMax + 2 sizes
  NCSFMax = max(1,nint((binom(NSOMOMax,(NSOMOMax+1)/2)-binom(NSOMOMax,((NSOMOMax+1)/2)+1)))) ! TODO: NCSFs for MS=0
  NBFMax = NCSFMax
  maxDetDimPerBF = max(1,nint((binom(NSOMOMax,(NSOMOMax+1)/2))))
  NMO = n_act_orb
  integer i,j,k,l
  integer startdet,enddet
  integer ncfg,ncfgprev
  integer NSOMO
  integer dimcsfpercfg
  integer detDimperBF
  real*8 :: coeff
  integer MS
  integer ncfgpersomo
  detDimperBF = 0
  MS = elec_alpha_num-elec_beta_num
  !print *,"NSOMOMax=",NSOMOMax, cfg_seniority_index(0)
  ! number of cfgs = number of dets for 0 somos
  dimBasisCSF = cfg_seniority_index(0)-1
  ncfgprev = cfg_seniority_index(0)
  do i = 0-iand(MS,1)+2, NSOMOMax,2
     if(cfg_seniority_index(i) .EQ. -1)then
        ncfgpersomo = N_configuration + 1
     else
        ncfgpersomo = cfg_seniority_index(i)
     endif
  ncfg = ncfgpersomo - ncfgprev
  !detDimperBF = max(1,nint((binom(i,(i+1)/2))))
  dimcsfpercfg = max(1,nint((binom(i-2,(i-2+1)/2)-binom(i-2,((i-2+1)/2)+1))))
  dimBasisCSF += ncfg * dimcsfpercfg
  !print *,i,">(",ncfg,ncfgprev,ncfgpersomo,")",",",detDimperBF,">",dimcsfpercfg, " | dimbas= ", dimBasisCSF
  !if(cfg_seniority_index(i+2) == -1) EXIT
  !if(detDimperBF > maxDetDimPerBF) maxDetDimPerBF = detDimperBF
  ncfgprev = cfg_seniority_index(i)
  enddo
  if(NSOMOMax .EQ. elec_num)then
        ncfgpersomo = N_configuration + 1
        ncfg = ncfgpersomo - ncfgprev
        dimcsfpercfg = max(1,nint((binom(i-2,(i-2+1)/2)-binom(i-2,((i-2+1)/2)+1))))
        dimBasisCSF += ncfg * dimcsfpercfg
        !print *,i,">(",ncfg,ncfgprev,ncfgpersomo,")",",",detDimperBF,">",dimcsfpercfg, " | dimbas= ", dimBasisCSF
  endif
  END_PROVIDER

subroutine get_phase_qp_to_cfg(Ialpha, Ibeta, phaseout)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Documentation for get_phase_qp_to_cfg
  !
  ! This function converts from (aaaa)(bbbb)
  ! notation to (ab)(ab)(ab)(ab)
  ! notation.
  ! The cfgCI code works in (ab)(ab)(ab)(ab)
  ! notation throughout.
  END_DOC
  integer(bit_kind),intent(in)   :: Ialpha(N_int)
  integer(bit_kind),intent(in)   :: Ibeta(N_int)
  real*8,intent(out)             :: phaseout
  integer(bit_kind)              :: mask, mask2(N_int), deta(N_int), detb(N_int)
  integer                        :: nbetas
  integer                        :: count, k

  ! Find how many alpha electrons there are in all the N_ints
  integer :: Na(N_int)
  do k=1,N_int
    Na(k) = popcnt(deta(k))
  enddo

  integer :: shift, ipos, nperm
  phaseout = 1.d0
  do k=1,N_int

    do while(detb(k) /= 0_bit_kind)
      ! Find the lowest beta electron and clear it
      ipos = trailz(detb(k))
      detb(k) = ibclr(detb(k),ipos)

      ! Create a mask will all MOs higher than the beta electron
      mask = not(shiftl(1_bit_kind,ipos+1) - 1_bit_kind)

      ! Apply the mask to the alpha string to count how many electrons to cross
      nperm = popcnt( iand(mask, deta(k)) )

      ! Count how many alpha electrons are above the beta electron in the other integers
      nperm = nperm + sum(Na(k+1:N_int))
      if (iand(nperm,1) == 1) then
        phaseout  = -phaseout
      endif

    enddo

  enddo
end subroutine get_phase_qp_to_cfg



  BEGIN_PROVIDER [ real*8, DetToCSFTransformationMatrix, (0:NSOMOMax,NBFMax,maxDetDimPerBF)]
 &BEGIN_PROVIDER [ real*8, psi_coef_config,  (dimBasisCSF,1)]
 &BEGIN_PROVIDER [ integer, psi_config_data, (N_configuration,2)]
 &BEGIN_PROVIDER [ integer, psi_csf_to_config_data, (dimBasisCSF)]
  use cfunctions
  implicit none
  BEGIN_DOC
  ! Documentation for DetToCSFTransformationMatrix
  ! Provides the matrix of transformatons for the
  ! conversion between determinant to CSF basis (in BFs)
  END_DOC
  integer*8 :: Isomo, Idomo, mask, Ialpha,Ibeta
  integer   :: rows, cols, i, j, k
  integer   :: startdet, enddet
  integer*8 MS
  integer ndetI
  integer :: getNSOMO
  real*8,dimension(:,:),allocatable    :: tempBuffer
  real*8,dimension(:),allocatable    :: tempCoeff
  real*8  :: norm_det1, phasedet
  norm_det1 = 0.d0
  MS = elec_alpha_num - elec_beta_num
  print *,"Maxbfdim=",NBFMax
  print *,"Maxdetdim=",maxDetDimPerBF
  print *,"dimBasisCSF=",dimBasisCSF
  print *,"N_configurations=",N_configuration
  print *,"n_core_orb=",n_core_orb
  ! initialization
  psi_coef_config = 0.d0
  DetToCSFTransformationMatrix(0,:,:) = 1.d0
  do i = 2-iand(elec_alpha_num-elec_beta_num,1), NSOMOMax,2
    Isomo = IBSET(0_8, i) - 1_8
    ! rows = Ncsfs
    ! cols = Ndets
    bfIcfg = max(1,nint((binom(i,(i+1)/2)-binom(i,((i+1)/2)+1))))
    ndetI = max(1,nint((binom(i,(i+1)/2))))

    allocate(tempBuffer(bfIcfg,ndetI))
    call getCSFtoDETTransformationMatrix(Isomo, MS, NBFMax, maxDetDimPerBF, tempBuffer)
    DetToCSFTransformationMatrix(i,:bfIcfg,:ndetI) =  tempBuffer
    deallocate(tempBuffer)
  enddo

  integer s, bfIcfg
  integer countcsf
  countcsf = 0
  integer countdet
  countdet = 0
  integer istate
  istate = 1
  phasedet = 1.0d0
  do i = 1,N_configuration
      startdet = psi_configuration_to_psi_det(1,i)
      enddet = psi_configuration_to_psi_det(2,i)
      ndetI = enddet-startdet+1

      allocate(tempCoeff(ndetI))
      countdet = 1
      do j = startdet, enddet
         Ialpha = psi_det(1,1,psi_configuration_to_psi_det_data(j))
         Ibeta  = psi_det(1,2,psi_configuration_to_psi_det_data(j))
         !call debug_spindet(Ialpha,1,1)
         !call debug_spindet(Ibeta ,1,1)
         call get_phase_qp_to_cfg(Ialpha, Ibeta, phasedet)
         !print *,">>",Ialpha,Ibeta,phasedet
         tempCoeff(countdet) = psi_coef(psi_configuration_to_psi_det_data(j), istate)*phasedet
         !tempCoeff(countdet) = psi_coef(psi_configuration_to_psi_det_data(j), istate)
         norm_det1 += tempCoeff(countdet)*tempCoeff(countdet)
         countdet += 1
      enddo

       !print *,"dimcoef=",bfIcfg,norm_det1
       !call printMatrix(tempCoeff,ndetI,1)

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
       !print *,"csftodetdim=",bfIcfg,ndetI
       !call printMatrix(tempBuffer,bfIcfg,ndetI)

       call dgemm('N','N', bfIcfg, 1, ndetI, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, size(tempCoeff,1), 0.d0, psi_coef_config(countcsf+1,1), size(psi_coef_config,1))
       !call dgemv('N', NBFMax, maxDetDimPerBF, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, 1, 0.d0, psi_coef_config(countcsf), 1)

      !call printMatrix(psi_coef_config(countcsf+1,1),bfIcfg,1)
      deallocate(tempCoeff)
      deallocate(tempBuffer)
      psi_config_data(i,1) = countcsf + 1
      countcsf += bfIcfg
      psi_config_data(i,2) = countcsf
      psi_csf_to_config_data(countcsf) = i
  enddo
  print *,"Norm det=",norm_det1, size(psi_coef_config,1), " Dim csf=", countcsf

  END_PROVIDER

  subroutine convertWFfromDETtoCSF(psi_coef_det_in, psi_coef_cfg_out)
  use cfunctions
  implicit none
  BEGIN_DOC
  ! Documentation for DetToCSFTransformationMatrix
  ! Provides the matrix of transformatons for the
  ! conversion between determinant to CSF basis (in BFs)
  END_DOC
  integer*8 :: Isomo, Idomo, mask, Ialpha,Ibeta
  integer   :: rows, cols, i, j, k
  integer   :: startdet, enddet
  integer*8 MS
  integer ndetI
  integer :: getNSOMO
  real*8,intent(in)    :: psi_coef_det_in(n_det,1)
  real*8,intent(out)    :: psi_coef_cfg_out(dimBasisCSF,1)
  real*8,dimension(:,:),allocatable    :: tempBuffer
  real*8,dimension(:),allocatable    :: tempCoeff
  real*8  :: norm_det1, phasedet
  norm_det1 = 0.d0
  MS = elec_alpha_num - elec_beta_num
  print *,"Maxbfdim=",NBFMax
  print *,"Maxdetdim=",maxDetDimPerBF
  print *,"dimBasisCSF=",dimBasisCSF
  print *,"N_configurations=",N_configuration
  print *,"n_core_orb=",n_core_orb
  ! initialization
  psi_coef_cfg_out(:,1) = 0.d0

  integer s, bfIcfg
  integer countcsf
  countcsf = 0
  integer countdet
  countdet = 0
  integer istate
  istate = 1
  phasedet = 1.0d0
  do i = 1,N_configuration
      startdet = psi_configuration_to_psi_det(1,i)
      enddet = psi_configuration_to_psi_det(2,i)
      ndetI = enddet-startdet+1

      allocate(tempCoeff(ndetI))
      countdet = 1
      do j = startdet, enddet
         Ialpha = psi_det(1,1,psi_configuration_to_psi_det_data(j))
         Ibeta  = psi_det(1,2,psi_configuration_to_psi_det_data(j))
         call get_phase_qp_to_cfg(Ialpha, Ibeta, phasedet)
         !print *,">>",Ialpha,Ibeta,phasedet
         tempCoeff(countdet) = psi_coef(psi_configuration_to_psi_det_data(j), istate)*phasedet
         !tempCoeff(countdet) = psi_coef(psi_configuration_to_psi_det_data(j), istate)
         norm_det1 += tempCoeff(countdet)*tempCoeff(countdet)
         countdet += 1
      enddo

       !print *,"dimcoef=",bfIcfg,norm_det1
       !call printMatrix(tempCoeff,ndetI,1)

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
       !print *,"csftodetdim=",bfIcfg,ndetI
       !call printMatrix(tempBuffer,bfIcfg,ndetI)

       call dgemm('N','N', bfIcfg, 1, ndetI, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, size(tempCoeff,1), 0.d0, psi_coef_cfg_out(countcsf+1,1), size(psi_coef_cfg_out,1))

      deallocate(tempCoeff)
      deallocate(tempBuffer)
      psi_config_data(i,1) = countcsf + 1
      countcsf += bfIcfg
      psi_config_data(i,2) = countcsf
  enddo
  print *,"Norm det=",norm_det1, size(psi_coef_cfg_out,1), " Dim csf=", countcsf

  end subroutine convertWFfromDETtoCSF

  subroutine convertWFfromCSFtoDET(psi_coef_cfg_in, psi_coef_det)
    implicit none
    BEGIN_DOC
    ! Documentation for convertCSFtoDET
    ! This function converts the wavefunction
    ! in CFG basis to DET basis using the
    ! transformation matrix provided before.
    END_DOC
    real*8,intent(in)  :: psi_coef_cfg_in(dimBasisCSF,1)
    real*8,intent(out) :: psi_coef_det(N_det,1)
    real*8             :: tmp_psi_coef_det(maxDetDimPerBF)
    integer s, bfIcfg
    integer countcsf
    integer countdet
    integer*8 :: Isomo, Idomo, Ialpha, Ibeta
    integer   :: rows, cols, i, j, k
    integer   :: startdet, enddet
    integer*8 MS
    integer ndetI
    integer :: getNSOMO
    real*8,dimension(:,:),allocatable    :: tempBuffer
    real*8,dimension(:),allocatable    :: tempCoeff
    real*8  :: phasedet
    ! number of states
    integer istate
    istate = 1
    countcsf = 1
    countdet = 1
    print *,"in function convertWFfromCSFtoDET()"


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

       allocate(tempCoeff(bfIcfg))

       do j = 1,bfIcfg
          tempCoeff(j) = psi_coef_cfg_in(countcsf,1)
          countcsf += 1
       enddo
       !print *,"dimcoef=",bfIcfg
       !call printMatrix(tempCoeff,bfIcfg,1)

       ! perhaps blocking with CFGs of same seniority
       ! can be more efficient
       allocate(tempBuffer(bfIcfg,ndetI))
       tempBuffer = DetToCSFTransformationMatrix(s,:bfIcfg,:ndetI)
       !print *,"csftodetdim=",bfIcfg,ndetI
       !call printMatrix(tempBuffer,bfIcfg,ndetI)

       !call dgemm('T','N', ndetI, 1, bfIcfg, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, size(tempCoeff,1), 0.d0, psi_coef_det(countdet,1), size(psi_coef_det,1))
       call dgemm('T','N', ndetI, 1, bfIcfg, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, size(tempCoeff,1), 0.d0, tmp_psi_coef_det, size(tmp_psi_coef_det,1))

       !call dgemv('N', NBFMax, maxDetDimPerBF, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, 1, 0.d0, psi_coef_config(countcsf,1), 1)

       !print *,"result"
       !call printMatrix(tmp_psi_coef_det,ndetI,1)

       countdet = 1
       do j=startdet,enddet
         Ialpha = psi_det(1,1,psi_configuration_to_psi_det_data(j))
         Ibeta  = psi_det(1,2,psi_configuration_to_psi_det_data(j))
         !call debug_spindet(Ialpha,1,1)
         !call debug_spindet(Ibeta ,1,1)
         call get_phase_qp_to_cfg(Ialpha, Ibeta, phasedet)
         !print *,">>",Ialpha,Ibeta,phasedet
         psi_coef_det(psi_configuration_to_psi_det_data(j),1) = tmp_psi_coef_det(countdet)*phasedet
         countdet += 1
       enddo

       deallocate(tempCoeff)
       deallocate(tempBuffer)
       !countdet += ndetI
    enddo

    !countdet = 1
    !tmp_psi_coef_det = psi_coef_det(:,1)
    !do i=1,N_configuration
    !   startdet = psi_configuration_to_psi_det(1,i)
    !   enddet = psi_configuration_to_psi_det(2,i)
    !   ndetI = enddet-startdet+1
    !   print *,i,">>>",startdet,enddet
    !   do k=1,ndetI
    !      !psi_coef_det(startdet+k-1,1) = tmp_psi_coef_det(countdet)
    !      psi_coef_det(countdet,1) = tmp_psi_coef_det(startdet+k-1)
    !      countdet += 1
    !   enddo
    !enddo

    print *,"End ncsfs=",countcsf

  end subroutine convertCSFtoDET

  BEGIN_PROVIDER [ integer, AIJpqMatrixDimsList, (0:NSOMOMax,0:NSOMOMax,4,NSOMOMax+1,NSOMOMax+1,2)]
 &BEGIN_PROVIDER [ integer, rowsmax]
 &BEGIN_PROVIDER [ integer, colsmax]
  use cfunctions
  implicit none
  BEGIN_DOC
  ! Documentation for AIJpqMatrixList
  ! The prototype matrix containing the <I|E_{pq}|J>
  ! matrices for each I,J somo pair and orb ids.
  END_DOC
  integer i,j,k,l
  integer*8 Isomo, Jsomo, tmpsomo
  Isomo = 0
  Jsomo = 0
  integer rows, cols, nsomoi, nsomoj
  rows = -1
  cols = -1
  integer*8 MS
  MS = elec_alpha_num-elec_beta_num
  integer nsomomin
  nsomomin = elec_alpha_num-elec_beta_num
  rowsmax = 0
  colsmax = 0
  print *,"NSOMOMax = ",NSOMOMax
  !allocate(AIJpqMatrixDimsList(NSOMOMax,NSOMOMax,4,NSOMOMax,NSOMOMax,2))
  ! Type
  ! 1. SOMO -> SOMO
  !print *,"Doing SOMO->SOMO"
  AIJpqMatrixDimsList(0,0,1,1,1,1) = 1
  AIJpqMatrixDimsList(0,0,1,1,1,2) = 1
  do i = 2-iand(nsomomin,1), NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i-2,i-2, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LT. 0) then
           cycle
        end if
        do k = 1,i
           do l = 1,i
              ! Define Jsomo
              if(k.NE.l)then
                 Jsomo = IBCLR(Isomo, k-1)
                 Jsomo = IBCLR(Jsomo, l-1)
                 nsomoi = i
                 nsomoj = j
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,i)-1
                 nsomoi = i
                 nsomoj = i
              endif

              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)
              !print *, "SOMO->SOMO \t",i,j,k,l,">",Isomo,Jsomo,">",rows, cols
              if(rowsmax .LT. rows) then
                 rowsmax = rows
              end if
              if(colsmax .LT. cols) then
                 colsmax = cols
              end if
              ! i -> j
              AIJpqMatrixDimsList(nsomoi,nsomoj,1,k,l,1) = rows
              AIJpqMatrixDimsList(nsomoi,nsomoj,1,k,l,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 2. DOMO -> VMO
  !print *,"Doing DOMO->VMO"
  AIJpqMatrixDimsList(0,0,2,1,1,1) = 1
  AIJpqMatrixDimsList(0,0,2,1,1,2) = 1
  do i = 0+iand(nsomomin,1), NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     tmpsomo = ISHFT(1_8,i+2)-1
     do j = i+2,i+2, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LT. 0) then
           cycle
        end if
        do k = 1,j
           do l = 1,j
              if(k .NE. l) then
                 Isomo = IBCLR(tmpsomo,k-1)
                 Isomo = IBCLR(Isomo,l-1)

                 ! Define Jsomo
                 Jsomo = ISHFT(1_8,j)-1;
                 nsomoi = i
                 nsomoj = j
              else
                 Isomo = ISHFT(1_8,j)-1
                 Isomo = IBCLR(Isomo,1-1)
                 Isomo = IBCLR(Isomo,j-1)
                 Jsomo = ISHFT(1_8,j)-1
                 Isomo = ISHFT(1_8,j)-1
                 nsomoi = j
                 nsomoj = j
              endif

              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)
              !print *, i,j,k,l,">",Isomo,Jsomo,">",rows, cols
              if(rowsmax .LT. rows) then
                 rowsmax = rows
              end if
              if(colsmax .LT. cols) then
                 colsmax = cols
              end if
              ! i -> j
              AIJpqMatrixDimsList(nsomoi,nsomoj,2,k,l,1) = rows
              AIJpqMatrixDimsList(nsomoi,nsomoj,2,k,l,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 3. SOMO -> VMO
  !print *,"Doing SOMO->VMO"
  AIJpqMatrixDimsList(0,0,3,1,1,1) = 1
  AIJpqMatrixDimsList(0,0,3,1,1,2) = 1
  do i = 2-iand(nsomomin,1), NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i,i, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,i+1
           do l = 1,i+1
              if(k .NE. l) then
                 Isomo = ISHFT(1_8,i+1)-1
                 Isomo = IBCLR(Isomo,l-1)
                 Jsomo = ISHFT(1_8,j+1)-1
                 Jsomo = IBCLR(Jsomo,k-1)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif
              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)
              !print *, i,j,k,l,">",Isomo,Jsomo,">",rows, cols
              if(rowsmax .LT. rows) then
                 rowsmax = rows
              end if
              if(colsmax .LT. cols) then
                 colsmax = cols
              end if
              ! i -> j
              AIJpqMatrixDimsList(i,j,3,k,l,1) = rows
              AIJpqMatrixDimsList(i,j,3,k,l,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 4. DOMO -> SOMO
  !print *,"Doing DOMO->SOMO"
  AIJpqMatrixDimsList(0,0,4,1,1,1) = 1
  AIJpqMatrixDimsList(0,0,4,1,1,2) = 1
  do i = 2-iand(nsomomin,1), NSOMOMax, 2
     do j = i,i, 2
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,i+1
           do l = 1,i+1
              if(k .NE. l) then
                 Isomo = ISHFT(1_8,i+1)-1
                 Isomo = IBCLR(Isomo,k-1)
                 Jsomo = ISHFT(1_8,j+1)-1
                 Jsomo = IBCLR(Jsomo,l-1)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif
              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)
              !print *, i,j,k,l,">",Isomo,Jsomo,">",rows, cols
              if(rowsmax .LT. rows) then
                 rowsmax = rows
              end if
              if(colsmax .LT. cols) then
                 colsmax = cols
              end if
              ! i -> j
              AIJpqMatrixDimsList(i,j,4,k,l,1) = rows
              AIJpqMatrixDimsList(i,j,4,k,l,2) = cols
           end do
        end do
     end do
  end do
  print *,"Rowsmax=",rowsmax," Colsmax=",colsmax
  END_PROVIDER

  BEGIN_PROVIDER [ real*8, AIJpqContainer, (0:NSOMOMax,0:NSOMOMax,4,NSOMOMax+1,NSOMOMax+1,NBFMax,NBFMax)]
  use cfunctions
  implicit none
  BEGIN_DOC
  ! Documentation for AIJpqMatrixList
  ! The prototype matrix containing the <I|E_{pq}|J>
  ! matrices for each I,J somo pair and orb ids.
  !
  ! Due to the different types of excitations which
  ! include DOMOs and VMOs two prototype DOMOs and two
  ! prototype VMOs are needed. Therefore
  ! the 4th and 5th dimensions are NSOMOMax+4 and NSOMOMax+4
  ! respectively.
  !
  ! The type of excitations are ordered as follows:
  ! Type 1 - SOMO -> SOMO
  ! Type 2 - DOMO -> VMO
  ! Type 3 - SOMO -> VMO
  ! Type 4 - DOMO -> SOMO
  END_DOC
  integer i,j,k,l, orbp, orbq, ri, ci
  orbp = 0
  orbq = 0
  integer*8 Isomo, Jsomo, tmpsomo
  Isomo = 0
  Jsomo = 0
  integer rows, cols, nsomoi, nsomoj
  rows = -1
  cols = -1
  integer*8 MS
  MS = 0
  touch AIJpqMatrixDimsList
  real*8,dimension(:,:),allocatable :: meMatrix
  integer maxdim
  !maxdim = max(rowsmax,colsmax)
  ! allocate matrix
  !print *,"rowsmax =",rowsmax," colsmax=",colsmax
  !print *,"NSOMOMax = ",NSOMOMax
  !allocate(AIJpqMatrixDimsList(NSOMOMax,NSOMOMax,4,NSOMOMax,NSOMOMax,2))
  ! Type
  ! 1. SOMO -> SOMO
  !print *,"Doing SOMO -> SOMO"
  AIJpqContainer(0,0,1,1,1,1,1) = 1.0d0
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i-2,i-2, 2
        if(j .GT. NSOMOMax .OR. j .LT. 0) cycle
        !print *,"i,j=",i,j
        do k = 1,i
           do l = 1,i

              ! Define Jsomo
              if(k .NE. l) then
                 Jsomo = IBCLR(Isomo, k-1)
                 Jsomo = IBCLR(Jsomo, l-1)
                 nsomoi = i
                 nsomoj = j
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,i)-1
                 nsomoi = i
                 nsomoj = i
              endif

              !print *,"k,l=",k,l
              !call debug_spindet(Jsomo,1)
              !call debug_spindet(Isomo,1)

              AIJpqContainer(nsomoi,nsomoj,1,k,l,:,:) = 0.0d0
              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)

              orbp = k
              orbq = l
              allocate(meMatrix(rows,cols))
              meMatrix = 0.0d0
              ! fill matrix
              call getApqIJMatrixDriver(Isomo,           &
                   Jsomo, &
                   orbp,                     &
                   orbq,                     &
                   MS,                       &
                   NMO,                      &
                   meMatrix,                 &
                   rows,                     &
                   cols)
             !print *, i,j,k,l,">",Isomo,Jsomo,">",rows, cols,">",rowsmax,colsmax
             !call printMatrix(meMatrix,rows,cols)
              ! i -> j
             do ri = 1,rows
                 do ci = 1,cols
                    AIJpqContainer(nsomoi,nsomoj,1,k,l,ri,ci) = meMatrix(ri, ci)
                 end do
              end do
              deallocate(meMatrix)
           end do
        end do
     end do
  end do
  ! Type
  ! 2. DOMO -> VMO
  !print *,"Doing DOMO -> VMO"
  AIJpqContainer(0,0,2,1,1,1,1) = 1.0d0
  do i = 0, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     tmpsomo = ISHFT(1_8,i+2)-1
     do j = i+2,i+2, 2
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        Jsomo = ISHFT(1_8,j)-1
        do k = 1,j
           do l = 1,j
              if(k .NE. l) then
                 Isomo = IBCLR(tmpsomo,k-1)
                 Isomo = IBCLR(Isomo,l-1)
                 ! Define Jsomo
                 Jsomo = ISHFT(1_8,j)-1;
                 nsomoi = i
                 nsomoj = j
              else
                 Isomo = ISHFT(1_8,j)-1
                 Isomo = IBCLR(Isomo,1-1)
                 Isomo = IBCLR(Isomo,j-1)
                 Jsomo = ISHFT(1_8,j)-1
                 Isomo = ISHFT(1_8,j)-1
                 nsomoi = j
                 nsomoj = j
              endif

              !print *,"k,l=",k,l
              !call debug_spindet(Jsomo,1)
              !call debug_spindet(Isomo,1)

              AIJpqContainer(nsomoi,nsomoj,2,k,l,:,:) = 0.0d0
              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)

              orbp = k
              orbq = l
              allocate(meMatrix(rows,cols))
              meMatrix = 0.0d0
              ! fill matrix
              call getApqIJMatrixDriver(Isomo,           &
                   Jsomo, &
                   orbp,                     &
                   orbq,                     &
                   MS,                       &
                   NMO,                      &
                   meMatrix,                 &
                   rows,                     &
                   cols)
             !print *, i,j,k,l,">",Isomo,Jsomo,">",rows, cols,">",rowsmax,colsmax
             !call printMatrix(meMatrix,rows,cols)
              ! i -> j
             do ri = 1,rows
                 do ci = 1,cols
                    AIJpqContainer(nsomoi,nsomoj,2,k,l,ri,ci) = meMatrix(ri, ci)
                 end do
              end do
              deallocate(meMatrix)
           end do
        end do
     end do
  end do
  ! Type
  ! 3. SOMO -> VMO
  !print *,"Doing SOMO -> VMO"
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i,i, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = 1,i+1
           do l = 1,i+1
              if(k .NE. l) then
                 Isomo = ISHFT(1_8,i+1)-1
                 Isomo = IBCLR(Isomo,l-1)
                 Jsomo = ISHFT(1_8,j+1)-1
                 Jsomo = IBCLR(Jsomo,k-1)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif

              !print *,"k,l=",k,l
              !call debug_spindet(Jsomo,1)
              !call debug_spindet(Isomo,1)

              AIJpqContainer(i,j,3,k,l,:,:) = 0.0d0
              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)

              orbp = k
              orbq = l
              allocate(meMatrix(rows,cols))
              meMatrix = 0.0d0
              ! fill matrix
              call getApqIJMatrixDriver(Isomo,           &
                   Jsomo, &
                   orbp,                     &
                   orbq,                     &
                   MS,                       &
                   NMO,                      &
                   meMatrix,                 &
                   rows,                     &
                   cols)
             !call printMatrix(meMatrix,rows,cols)
             !print *, i,j,k,l,">",Isomo,Jsomo,">",rows, cols,">",rowsmax,colsmax
              ! i -> j
             do ri = 1,rows
                 do ci = 1,cols
                    AIJpqContainer(i,j,3,k,l,ri,ci) = meMatrix(ri, ci)
                 end do
              end do
              deallocate(meMatrix)
           end do
        end do
     end do
  end do
  ! Type
  ! 4. DOMO -> SOMO
  !print *,"Doing DOMO -> SOMO"
  AIJpqContainer(0,0,4,1,1,1,1) = 1.0d0
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i,i, 2
        Jsomo = ISHFT(1_8,i)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = 1,i+1
           do l = 1,i+1
              if(k .NE. l) then
                 Isomo = ISHFT(1_8,i+1)-1
                 Isomo = IBCLR(Isomo,k-1)
                 Jsomo = ISHFT(1_8,j+1)-1
                 Jsomo = IBCLR(Jsomo,l-1)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif

              AIJpqContainer(i,j,4,k,l,:,:) = 0.0d0
              call getApqIJMatrixDims(Isomo,           &
                   Jsomo, &
                   MS,                       &
                   rows,                     &
                   cols)

              orbp = k
              orbq = l

              allocate(meMatrix(rows,cols))
              meMatrix = 0.0d0
              ! fill matrix
              call getApqIJMatrixDriver(Isomo,           &
                   Jsomo, &
                   orbp,                     &
                   orbq,                     &
                   MS,                       &
                   NMO,                      &
                   meMatrix,                 &
                   rows,                     &
                   cols)
             !call printMatrix(meMatrix,rows,cols)
             !print *, i,j,k,l,">",Isomo,Jsomo,">",rows, cols,">",rowsmax,colsmax
              ! i -> j
             do ri = 1,rows
                 do ci = 1,cols
                    AIJpqContainer(i,j,4,k,l,ri,ci) = meMatrix(ri, ci)
                 end do
              end do
              deallocate(meMatrix)
           end do
        end do
     end do
  end do
  END_PROVIDER
