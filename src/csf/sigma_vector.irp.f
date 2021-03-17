  BEGIN_PROVIDER [ integer, NSOMOMax]
 &BEGIN_PROVIDER [ integer, NSOMOMin]
 &BEGIN_PROVIDER [ integer, NCSFMax]
 &BEGIN_PROVIDER [ integer*8, NMO]
 &BEGIN_PROVIDER [ integer, NBFMax]
 &BEGIN_PROVIDER [ integer, n_CSF]
 &BEGIN_PROVIDER [ integer, maxDetDimPerBF]
  implicit none
  BEGIN_DOC
  ! Documentation for NSOMOMax
  ! The maximum number of SOMOs for the current calculation.
  ! required for the calculation of prototype arrays.
  END_DOC
  NSOMOMax = min(elec_num, cfg_nsomo_max + 2)
  NSOMOMin = max(0,cfg_nsomo_min-2)
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
  print *,"NSOMOMax=",NSOMOMax, cfg_seniority_index(NSOMOMin)
  print *,"NSOMOMin=",NSOMOMin
  do i=0,elec_num
    print *,i," #cfgs=",cfg_seniority_index(i)
  enddo
  ! number of cfgs = number of dets for 0 somos
  n_CSF = cfg_seniority_index(NSOMOMin)-1
  print *,"start=",n_CSF
  ncfgprev = cfg_seniority_index(NSOMOMin)
  !do i = 0-iand(MS,1)+2, NSOMOMax,2
  do i = NSOMOMin+2, NSOMOMax,2
     if(cfg_seniority_index(i) .EQ. -1)then
        ncfgpersomo = N_configuration + 1
     else
        ncfgpersomo = cfg_seniority_index(i)
     endif
  ncfg = ncfgpersomo - ncfgprev
  !detDimperBF = max(1,nint((binom(i,(i+1)/2))))
  dimcsfpercfg = max(1,nint((binom(i-2,(i-2+1)/2)-binom(i-2,((i-2+1)/2)+1))))
  n_CSF += ncfg * dimcsfpercfg
  print *,i,">(",ncfg,ncfgprev,ncfgpersomo,")",",",detDimperBF,">",dimcsfpercfg, " | dimbas= ", n_CSF
  !if(cfg_seniority_index(i+2) == -1) EXIT
  !if(detDimperBF > maxDetDimPerBF) maxDetDimPerBF = detDimperBF
  ncfgprev = cfg_seniority_index(i)
  enddo
  !if(NSOMOMax .EQ. elec_num)then
  !      ncfgpersomo = N_configuration + 1
  !      ncfg = ncfgpersomo - ncfgprev
  !      dimcsfpercfg = max(1,nint((binom(i-2,(i-2+1)/2)-binom(i-2,((i-2+1)/2)+1))))
  !      n_CSF += ncfg * dimcsfpercfg
  !      print *,i,">>(",ncfg,ncfgprev,ncfgpersomo,")",",",detDimperBF,">",dimcsfpercfg, " | dimbas= ", n_CSF
  !endif
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
  integer(bit_kind)              :: mask, deta(N_int), detb(N_int)
  integer                        :: nbetas
  integer                        :: count, k

  ! Initliaze deta and detb
  deta = Ialpha
  detb = Ibeta

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
      mask = not(shiftl(1_bit_kind,ipos + 1) - 1_bit_kind)

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
 &BEGIN_PROVIDER [ real*8, psi_coef_config,  (n_CSF,1)]
 &BEGIN_PROVIDER [ integer, psi_config_data, (N_configuration,2)]
 &BEGIN_PROVIDER [ integer, psi_csf_to_config_data, (n_CSF)]
  use cfunctions
  implicit none
  BEGIN_DOC
  ! Documentation for DetToCSFTransformationMatrix
  ! Provides the matrix of transformatons for the
  ! conversion between determinant to CSF basis (in BFs)
  END_DOC
  integer*8 :: Isomo, Idomo
  integer(bit_kind) :: Ialpha(N_int),Ibeta(N_int)
  integer   :: rows, cols, i, j, k
  integer   :: startdet, enddet
  integer*8 MS
  integer ndetI
  integer :: getNSOMO
  real*8,dimension(:,:),allocatable    :: tempBuffer
  real*8,dimension(:),allocatable    :: tempCoeff
  real*8  :: norm_det1, phasedet

  integer :: nt

  norm_det1 = 0.d0
  MS = elec_alpha_num - elec_beta_num
  print *,"Maxbfdim=",NBFMax
  print *,"Maxdetdim=",maxDetDimPerBF
  print *,"n_CSF=",n_CSF
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
  psi_csf_to_config_data(1) = 1
  phasedet = 1.0d0
  call omp_set_max_active_levels(1)
  !$OMP PARALLEL
  !$OMP MASTER
  do i = 1,N_configuration
      startdet = psi_configuration_to_psi_det(1,i)
      enddet = psi_configuration_to_psi_det(2,i)
      ndetI = enddet-startdet+1

      allocate(tempCoeff(ndetI))
      countdet = 1
      do j = startdet, enddet
         Ialpha = psi_det(:,1,psi_configuration_to_psi_det_data(j))
         Ibeta  = psi_det(:,2,psi_configuration_to_psi_det_data(j))
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
      do k=1,bfIcfg
        psi_csf_to_config_data(countcsf+k) = i
      enddo
      countcsf += bfIcfg
      psi_config_data(i,2) = countcsf
  enddo
  print *,"Norm det=",norm_det1, size(psi_coef_config,1), " Dim csf=", countcsf
  !$OMP END MASTER
  !$OMP END PARALLEL
  call omp_set_max_active_levels(4)

  END_PROVIDER

  BEGIN_PROVIDER [ integer, AIJpqMatrixDimsList, (NSOMOMin:NSOMOMax,4,NSOMOMax+1,NSOMOMax+1,2)]
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
  nsomomin = elec_alpha_num-elec_beta_num
  rowsmax = 0
  colsmax = 0
  print *,"NSOMOMax = ",NSOMOMax
  print *,"NSOMOMin = ",NSOMOMin
  !allocate(AIJpqMatrixDimsList(NSOMOMax,NSOMOMax,4,NSOMOMax,NSOMOMax,2))
  ! Type
  ! 1. SOMO -> SOMO
  !print *,"Doing SOMO->SOMO"
  AIJpqMatrixDimsList(NSOMOMin,1,1,1,1) = 1
  AIJpqMatrixDimsList(NSOMOMin,1,1,1,2) = 1
  do i = NSOMOMin, NSOMOMax, 2
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
              AIJpqMatrixDimsList(nsomoi,1,k,l,1) = rows
              AIJpqMatrixDimsList(nsomoi,1,k,l,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 2. DOMO -> VMO
  !print *,"Doing DOMO->VMO"
  AIJpqMatrixDimsList(NSOMOMin,2,1,1,1) = 1
  AIJpqMatrixDimsList(NSOMOMin,2,1,1,2) = 1
  do i = NSOMOMin, NSOMOMax, 2
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
              AIJpqMatrixDimsList(nsomoi,2,k,l,1) = rows
              AIJpqMatrixDimsList(nsomoi,2,k,l,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 3. SOMO -> VMO
  !print *,"Doing SOMO->VMO"
  AIJpqMatrixDimsList(NSOMOMin,3,1,1,1) = 1
  AIJpqMatrixDimsList(NSOMOMin,3,1,1,2) = 1
  do i = NSOMOMin, NSOMOMax, 2
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
              AIJpqMatrixDimsList(i,3,k,l,1) = rows
              AIJpqMatrixDimsList(i,3,k,l,2) = cols
           end do
        end do
     end do
  end do
  ! Type
  ! 4. DOMO -> SOMO
  !print *,"Doing DOMO->SOMO"
  AIJpqMatrixDimsList(NSOMOMin,4,1,1,1) = 1
  AIJpqMatrixDimsList(NSOMOMin,4,1,1,2) = 1
  do i = NSOMOMin, NSOMOMax, 2
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
              AIJpqMatrixDimsList(i,4,k,l,1) = rows
              AIJpqMatrixDimsList(i,4,k,l,2) = cols
           end do
        end do
     end do
  end do
  print *,"Rowsmax=",rowsmax," Colsmax=",colsmax
  END_PROVIDER

  !BEGIN_PROVIDER [ real*8, AIJpqContainer, (NSOMOMin:NSOMOMax,4,NSOMOMax+1,NSOMOMax+1,NBFMax,NBFMax)]
  BEGIN_PROVIDER [ real*8, AIJpqContainer, (NBFMax,NBFmax,NSOMOMax+1,NSOMOMax+1,4,NSOMOMin:NSOMOMax)]
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
  ! Type
  ! 1. SOMO -> SOMO
  !print *,"Doing SOMO -> SOMO"
  !AIJpqContainer(NSOMOMin,1,1,1,1,1) = 1.0d0
  AIJpqContainer(1,1,1,1,1,NSOMOMin) = 1.0d0
  do i = NSOMOMin+2, NSOMOMax, 2
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

              !AIJpqContainer(nsomoi,1,k,l,:,:) = 0.0d0
              AIJpqContainer(:,:,k,l,1,nsomoi) = 0.0d0
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
                    !AIJpqContainer(nsomoi,1,k,l,ri,ci) = meMatrix(ri, ci)
                    AIJpqContainer(ri,ci,k,l,1,nsomoi) = meMatrix(ri, ci)
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
  !AIJpqContainer(NSOMOMin,2,1,1,1,1) = 1.0d0
  AIJpqContainer(1,1,1,1,2,NSOMOMin) = 1.0d0
  do i = NSOMOMin, NSOMOMax, 2
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

              !AIJpqContainer(nsomoi,2,k,l,:,:) = 0.0d0
              AIJpqContainer(:,:,k,l,2,nsomoi) = 0.0d0
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
                    !AIJpqContainer(nsomoi,2,k,l,ri,ci) = meMatrix(ri, ci)
                    AIJpqContainer(ri,ci,k,l,2,nsomoi) = meMatrix(ri, ci)
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
  !AIJpqContainer(NSOMOMin,3,1,1,1,1) = 1.0d0
  AIJpqContainer(1,1,1,1,3,NSOMOMin) = 1.0d0
  do i = NSOMOMin, NSOMOMax, 2
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

              !AIJpqContainer(i,3,k,l,:,:) = 0.0d0
              AIJpqContainer(:,:,k,l,3,i) = 0.0d0
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
                    !AIJpqContainer(i,3,k,l,ri,ci) = meMatrix(ri, ci)
                    AIJpqContainer(ri,ci,k,l,3,i) = meMatrix(ri, ci)
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
  !AIJpqContainer(NSOMOMin,4,1,1,1,1) = 1.0d0
  AIJpqContainer(1,1,1,1,4,NSOMOMin) = 1.0d0
  do i = NSOMOMin+2, NSOMOMax, 2
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

              !AIJpqContainer(i,4,k,l,:,:) = 0.0d0
              AIJpqContainer(:,:,k,l,4,i) = 0.0d0
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
                    !AIJpqContainer(i,4,k,l,ri,ci) = meMatrix(ri, ci)
                    AIJpqContainer(ri,ci,k,l,4,i) = meMatrix(ri, ci)
                 end do
              end do
              deallocate(meMatrix)
           end do
        end do
     end do
  end do
  END_PROVIDER

subroutine calculate_preconditioner_cfg(diag_energies)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Documentation for calculate_preconditioner
  !
  ! Calculates the diagonal energies of
  ! the configurations in psi_configuration
  ! returns : diag_energies :
  END_DOC
  integer :: i,j,k,l,p,q,noccp,noccq, ii, jj
  real*8,intent(out) :: diag_energies(n_CSF)
  integer                            :: nholes
  integer                            :: nvmos
  integer                            :: listvmos(mo_num)
  integer                            :: vmotype(mo_num) ! 1 -> VMO 2 -> SOMO
  integer                            :: listholes(mo_num)
  integer                            :: holetype(mo_num) ! 1-> SOMO 2->DOMO
  integer*8                          :: Idomo
  integer*8                          :: Isomo
  integer*8                          :: Jdomo
  integer*8                          :: Jsomo
  integer*8                          :: diffSOMO
  integer*8                          :: diffDOMO
  integer                            :: NSOMOI
  integer                            :: NSOMOJ
  integer                            :: ndiffSOMO
  integer                            :: ndiffDOMO
  integer                            :: starti, endi, cnti, cntj, rows,cols
  integer                            :: extype,pmodel,qmodel
  integer(bit_kind) :: Icfg(N_INT,2)
  integer(bit_kind) :: Jcfg(N_INT,2)
  integer,external  :: getNSOMO
  real*8, external  :: mo_two_e_integral
  real*8            :: hpp
  real*8            :: meCC
  real*8            :: ecore

  ! initialize energies
  diag_energies = 0.d0

  ! calculate core energy
  !call get_core_energy(ecore)
  !diag_energies = ecore

  ! calculate the core energy
  !print *,"Core energy=",ref_bitmask_energy

  do i=1,N_configuration

     Isomo = psi_configuration(1,1,i)
     Idomo = psi_configuration(1,2,i)
     Icfg(1,1) = psi_configuration(1,1,i)
     Icfg(1,2) = psi_configuration(1,2,i)
     NSOMOI = getNSOMO(psi_configuration(:,:,i))

     starti = psi_config_data(i,1)
     endi   = psi_config_data(i,2)

     ! find out all pq holes possible
     nholes = 0
     ! holes in SOMO
     !do k = n_core_orb+1,n_core_orb + n_act_orb
     do k = 1,mo_num
        if(POPCNT(IAND(Isomo,IBSET(0_8,k-1))) .EQ. 1) then
           nholes += 1
           listholes(nholes) = k
           holetype(nholes) = 1
        endif
     enddo
     ! holes in DOMO
     !do k = n_core_orb+1,n_core_orb + n_act_orb
     !do k = 1+n_core_inact_orb,n_core_orb+n_core_inact_act_orb
     do k = 1,mo_num
        if(POPCNT(IAND(Idomo,IBSET(0_8,k-1))) .EQ. 1) then
           nholes += 1
           listholes(nholes) = k
           holetype(nholes) = 2
        endif
     enddo

     ! find vmos
     listvmos = -1
     vmotype = -1
     nvmos = 0
     !do k = n_core_orb+1,n_core_orb + n_act_orb
     do k = 1,mo_num
        !print *,i,IBSET(0,i-1),POPCNT(IAND(Isomo,(IBSET(0_8,i-1)))), POPCNT(IAND(Idomo,(IBSET(0_8,i-1))))
        if(POPCNT(IAND(Isomo,(IBSET(0_8,k-1)))) .EQ. 0 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,k-1)))) .EQ. 0) then
           nvmos += 1
           listvmos(nvmos) = k
           vmotype(nvmos) = 0
        else if(POPCNT(IAND(Isomo,(IBSET(0_8,k-1)))) .EQ. 1 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,k-1)))) .EQ. 0 ) then
           nvmos += 1
           listvmos(nvmos) = k
           vmotype(nvmos) = 1
        end if
     enddo
     !print *,"I=",i
     !call debug_spindet(psi_configuration(1,1,i),N_int)
     !call debug_spindet(psi_configuration(1,2,i),N_int)

     do k=1,nholes
        p = listholes(k)
        noccp = holetype(k)

        ! Calculate one-electron
        ! and two-electron coulomb terms
        do l=1,nholes
           q = listholes(l)
           noccq = holetype(l)
           !print *,"--------------- K=",p," L=",q

           ! one-electron term
           if(p.EQ.q) then
              hpp = noccq * h_core_ri(p,q)!mo_one_e_integrals(q,q)
           else
              hpp = 0.d0
           endif


           do j=starti,endi
              ! coulomb term
              ! (pp,qq) = <pq|pq>
              if(p.EQ.q) then
                 diag_energies(j) += hpp !+ 0.5d0 * (noccp * noccq * mo_two_e_integral(p,q,p,q))
                 !print *,"hpp=",hpp,"diga= ",diag_energies(j)
!             else
!                diag_energies(j) +=     !  0.5d0 * noccp * noccq * mo_two_e_integral(p,q,p,q)
!                print *,"diga= ",diag_energies(j)
              endif
           enddo
        enddo

     enddo
  enddo

end subroutine calculate_preconditioner_cfg


subroutine calculate_sigma_vector_cfg_nst_naive_store(psi_out, psi_in, n_st, sze, istart, iend, ishift, istep)
  implicit none
  use bitmasks
  use omp_lib
  BEGIN_DOC
  ! Documentation for sigma-vector calculation
  !
  ! Calculates the result of the
  ! application of the hamiltonian to the
  ! wavefunction in CFG basis once
  ! TODO : Things prepare outside this routine
  !  1. Touch the providers for
  !     a. ApqIJ containers
  !     b. DET to CSF transformation matrices
  !  2. DET to CSF transcormation
  !  2. CSF to DET back transcormation
  ! returns : psi_coef_out_det :
  END_DOC
  integer,intent(in)             :: sze, istart,iend, istep, ishift, n_st
  real*8,intent(in)              :: psi_in(n_st,sze)
  real*8,intent(out)             :: psi_out(n_st,sze)
  integer(bit_kind)              :: Icfg(N_INT,2)
  integer                        :: i,j,k,l,p,q,noccp,noccq, ii, jj, m, n, idxI, kk, nocck,orbk
  integer(bit_kind),dimension(:,:,:),allocatable :: listconnectedJ
  integer(bit_kind),dimension(:,:,:),allocatable :: alphas_Icfg
  integer(bit_kind),dimension(:,:,:),allocatable :: singlesI
  integer(bit_kind),dimension(:,:,:),allocatable :: connectedI_alpha
  integer,dimension(:),allocatable :: idxs_singlesI
  integer,dimension(:),allocatable :: idxs_connectedI_alpha
  integer,dimension(:,:),allocatable :: excitationIds_single
  integer,dimension(:),allocatable :: excitationTypes_single
  integer,dimension(:,:),allocatable :: excitationIds
  integer,dimension(:),allocatable :: excitationTypes
  integer,dimension(:),allocatable :: idslistconnectedJ
  real*8,dimension(:),allocatable :: diagfactors
  integer                        :: nholes
  integer                        :: nvmos
  integer                        :: listvmos(mo_num)
  integer                        :: vmotype(mo_num) ! 1 -> VMO 2 -> SOMO
  integer                        :: listholes(mo_num)
  integer                        :: holetype(mo_num) ! 1-> SOMO 2->DOMO
  integer                        :: Nalphas_Icfg, nconnectedI, rowsikpq, colsikpq, nsinglesI
  integer                        :: extype,NSOMOalpha,NSOMOI,NSOMOJ,pmodel,qmodel
  integer                        :: getNSOMO
  integer                        :: totcolsTKI
  integer                        :: rowsTKI
  integer                        :: noccpp
  integer                        :: istart_cfg, iend_cfg
  integer                        :: nconnectedJ
  integer*8                      :: MS, Isomo, Idomo, Jsomo, Jdomo, Ialpha, Ibeta
  integer                        :: moi, moj, mok, mol, starti, endi, startj, endj, cnti, cntj, cntk
  real*8                         :: norm_coef_cfg, fac2eints
  real*8                         :: norm_coef_det
  real*8                         :: meCC1, meCC2, diagfac
  real*8,dimension(:,:,:),allocatable :: TKI
  real*8,dimension(:,:),allocatable :: GIJpqrs
  real*8,dimension(:,:,:),allocatable :: TKIGIJ
  real*8,dimension(:),allocatable :: psi_out_tmp
  real*8,dimension(:,:),allocatable :: CCmattmp
  real*8, external               :: mo_two_e_integral
  real*8, external               :: get_two_e_integral
  real*8                         :: diag_energies(n_CSF)
  real*8                         :: tmpvar, tmptot

  integer(omp_lock_kind), allocatable :: lock(:)
  call omp_set_max_active_levels(1)

  allocate(lock(sze))
  do i=1,sze
    call omp_init_lock(lock(i))
  enddo

  !print *," sze = ",sze
  call calculate_preconditioner_cfg(diag_energies)

  MS = 0
  norm_coef_cfg=0.d0

  psi_out=0.d0

  istart_cfg = psi_csf_to_config_data(istart)
  iend_cfg   = psi_csf_to_config_data(iend)


  call omp_set_max_active_levels(1)
  !$OMP PARALLEL                                                     &
      !$OMP DEFAULT(NONE)                                            &
      !$OMP private(i,icfg, isomo, idomo, NSOMOI, NSOMOJ, nholes, k, listholes,&
      !$OMP    holetype, vmotype, nvmos, listvmos, starti, endi,     &
      !$OMP    nsinglesI, singlesI,idxs_singlesI,excitationIds_single,&
      !$OMP    excitationTypes_single, idxI, p, q, extype, pmodel, qmodel,&
      !$OMP    Jsomo, Jdomo, startj, endj, kk, jj, ii, cnti, cntj, meCC1,&
      !$OMP   nconnectedJ,listconnectedJ,idslistconnectedJ,             &
      !$OMP   Nalphas_Icfg,alphas_Icfg,connectedI_alpha,             &
      !$OMP   idxs_connectedI_alpha,nconnectedI,excitationIds,excitationTypes,diagfactors,&
      !$OMP   totcolsTKI,rowsTKI,NSOMOalpha,rowsikpq,                &
      !$OMP   colsikpq, GIJpqrs,TKIGIJ,j,l,m,TKI,CCmattmp, moi, moj, mok, mol,&
      !$OMP   diagfac, tmpvar)                                            &
      !$OMP shared(istart_cfg, iend_cfg, psi_configuration, mo_num, psi_config_data,&
      !$OMP    N_int, N_st, psi_out, psi_in, h_core_ri, AIJpqContainer,&
      !$OMP     sze, NalphaIcfg_list,alphasIcfg_list,                &
      !$OMP     AIJpqMatrixDimsList, diag_energies, n_CSF, lock)

  allocate(singlesI(N_INT,2,max(sze,100)))
  allocate(idxs_singlesI(max(sze,100)))
  allocate(excitationIds_single(2,max(sze,100)))
  allocate(excitationTypes_single(max(sze,100)))
!

  !!! Single Excitations !!!

  !$OMP DO SCHEDULE(dynamic,16)
  do i=istart_cfg,iend_cfg

    ! if Seniority_range > 8 then
    ! continue
    ! else
    ! cycle

     Icfg(1,1) = psi_configuration(1,1,i)
     Icfg(1,2) = psi_configuration(1,2,i)
     Isomo = Icfg(1,1)
     Idomo = Icfg(1,2)
     NSOMOI = getNSOMO(Icfg)

     ! find out all pq holes possible
     nholes = 0
     ! holes in SOMO
     ! list_act
     ! list_core
     ! list_core_inact
     ! bitmasks
     !do k = n_core_orb+1,n_core_orb + n_act_orb
     do k = 1,mo_num
       if(POPCNT(IAND(Isomo,IBSET(0_8,k-1))) .EQ. 1) then
         nholes += 1
         listholes(nholes) = k
         holetype(nholes) = 1
       endif
     enddo
     ! holes in DOMO
     !do k = n_core_orb+1,n_core_orb + n_act_orb
     do k = 1,mo_num
       if(POPCNT(IAND(Idomo,IBSET(0_8,k-1))) .EQ. 1) then
         nholes += 1
         listholes(nholes) = k
         holetype(nholes) = 2
       endif
     enddo

     ! find vmos
     listvmos = -1
     vmotype = -1
     nvmos = 0
     !do k = n_core_orb+1,n_core_orb + n_act_orb
     do k = 1,mo_num
       !print *,i,IBSET(0,i-1),POPCNT(IAND(Isomo,(IBSET(0_8,i-1)))), POPCNT(IAND(Idomo,(IBSET(0_8,i-1))))
       if(POPCNT(IAND(Isomo,(IBSET(0_8,k-1)))) .EQ. 0 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,k-1)))) .EQ. 0) then
         nvmos += 1
         listvmos(nvmos) = k
         vmotype(nvmos) = 0
       else if(POPCNT(IAND(Isomo,(IBSET(0_8,k-1)))) .EQ. 1 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,k-1)))) .EQ. 0 ) then
         nvmos += 1
         listvmos(nvmos) = k
         vmotype(nvmos) = 1
       end if
     enddo
 
 
     ! Icsf ids
     starti = psi_config_data(i,1)
     endi   = psi_config_data(i,2)
     NSOMOI = getNSOMO(Icfg)
 
     call generate_all_singles_cfg_with_type(Icfg,singlesI,idxs_singlesI,excitationIds_single,&
         excitationTypes_single,nsinglesI,N_int)
 
     do j = 1,nsinglesI
       idxI = idxs_singlesI(j)
       NSOMOJ = getNSOMO(singlesI(1,1,j))
      p = excitationIds_single(1,j)
      q = excitationIds_single(2,j)
      extype = excitationTypes_single(j)
      ! Off diagonal terms
      call convertOrbIdsToModelSpaceIds(Icfg, singlesI(1,1,j), p, q, extype, pmodel, qmodel)
      Jsomo = singlesI(1,1,j)
      Jdomo = singlesI(1,2,j)

      ! Add the hole on J
      if(POPCNT(IAND(Jsomo,IBSET(0_8,q-1))) .EQ. 1  .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
        nholes += 1
        listholes(nholes) = q
        holetype(nholes) = 1
      endif
      if((POPCNT(IAND(Jdomo,IBSET(0_8,q-1))) .EQ. 1 .AND. POPCNT(IAND(Idomo,IBSET(0_8,q-1))) .EQ. 0) .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
        nholes += 1
        listholes(nholes) = q
        holetype(nholes) = 2
      endif

      startj = psi_config_data(idxI,1)
      endj   = psi_config_data(idxI,2)

      !!! One-electron contribution !!!
      do ii = starti, endi
        cnti = ii-starti+1
        do jj = startj, endj
          cntj = jj-startj+1
          !meCC1 = AIJpqContainer(NSOMOI,extype,pmodel,qmodel,cnti,cntj)
          meCC1 = AIJpqContainer(cnti,cntj,pmodel,qmodel,extype,NSOMOI)
          call omp_set_lock(lock(jj))
          do kk = 1,n_st
            psi_out(kk,jj) += meCC1 * psi_in(kk,ii) * h_core_ri(p,q)
          enddo
          call omp_unset_lock(lock(jj))
        enddo
      enddo

      ! Undo setting in listholes
      if(POPCNT(IAND(Jsomo,IBSET(0_8,q-1))) .EQ. 1  .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
        nholes -= 1
      endif
      if((POPCNT(IAND(Jdomo,IBSET(0_8,q-1))) .EQ. 1 .AND. POPCNT(IAND(Idomo,IBSET(0_8,q-1))) .EQ. 0) .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
        nholes -= 1
      endif
    enddo
  enddo
  !$OMP END DO
  deallocate(singlesI)
  deallocate(idxs_singlesI)
  deallocate(excitationIds_single)
  deallocate(excitationTypes_single)

  allocate(listconnectedJ(N_INT,2,max(sze,100)))
  allocate(alphas_Icfg(N_INT,2,max(sze,100)))
  allocate(connectedI_alpha(N_INT,2,max(sze,100)))
  allocate(idxs_connectedI_alpha(max(sze,100)))
  allocate(excitationIds(2,max(sze,100)))
  allocate(excitationTypes(max(sze,100)))
  allocate(diagfactors(max(sze,100)))
  allocate(idslistconnectedJ(max(sze,100)))

  ! Loop over all selected configurations
  !$OMP DO SCHEDULE(dynamic,16)
  do i = istart_cfg,iend_cfg

     ! if Seniority_range > 8 then
     ! continue
     ! else
     ! cycle

     Icfg(1,1) = psi_configuration(1,1,i)
     Icfg(1,2) = psi_configuration(1,2,i)
     starti = psi_config_data(i,1)
     endi   = psi_config_data(i,2)

     ! Returns all unique (checking the past) singly excited cfgs connected to I
     Nalphas_Icfg = 0
     ! TODO:
     ! test if size(alphas_Icfg,1) < Nmo**2) then deallocate + allocate
     !call obtain_associated_alphaI(i, Icfg, alphas_Icfg, Nalphas_Icfg)
     Nalphas_Icfg = NalphaIcfg_list(i)
     alphas_Icfg(1:n_int,1:2,1:Nalphas_Icfg) = alphasIcfg_list(1:n_int,1:2,i,1:Nalphas_Icfg)
     !print *,"I=",i," Nal=",Nalphas_Icfg
     call obtain_connected_J_givenI(i, Icfg, listconnectedJ, idslistconnectedJ, nconnectedJ)
     !print *,"size listconnected=",size(listconnectedJ)
     !do k=1,nconnectedJ
     !  print *," idJ =",idslistconnectedJ(k)
     !enddo

     ! TODO : remove doubly excited for return
     ! Here we do 2x the loop. One to count for the size of the matrix, then we compute.
     do k = 1,Nalphas_Icfg
        ! Now generate all singly excited with respect to a given alpha CFG
        !call obtain_connected_I_foralpha(i,alphas_Icfg(1,1,k),connectedI_alpha,idxs_connectedI_alpha,nconnectedI,excitationIds,excitationTypes,diagfactors)
        call obtain_connected_I_foralpha_fromfilterdlist(i,nconnectedJ, idslistconnectedJ, listconnectedJ, alphas_Icfg(1,1,k),connectedI_alpha,idxs_connectedI_alpha,nconnectedI,excitationIds,excitationTypes,diagfactors)
        !print *,"\t---Ia=",k," NconI=",nconnectedI

        if(nconnectedI .EQ. 0) then
           cycle
        endif

        totcolsTKI = 0
        rowsTKI = -1
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(1,1,k), connectedI_alpha(1,1,j), p, q, extype, pmodel, qmodel)
           ! for E_pp E_rs and E_ppE_rr case
           if(p.EQ.q) then
              NSOMOalpha = NSOMOI
           endif
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,2)
           totcolsTKI += colsikpq
           if(rowsTKI .LT. rowsikpq .AND. rowsTKI .NE. -1) then
              print *,">",j,"Something is wrong in sigma-vector", rowsTKI, rowsikpq, "(p,q)=",pmodel,qmodel,"ex=",extype,"na=",NSOMOalpha," nI=",NSOMOI
              !rowsTKI = rowsikpq
           else
              rowsTKI = rowsikpq
           endif
        enddo

        allocate(TKI(n_st,rowsTKI,totcolsTKI)) ! coefficients of CSF
        ! Initialize the inegral container
        ! dims : (totcolsTKI, nconnectedI)
        allocate(GIJpqrs(totcolsTKI,nconnectedI))  ! gpqrs
        allocate(TKIGIJ(n_st,rowsTKI,nconnectedI))  ! TKI * gpqrs
        !print *,"\t---rowsTKI=",rowsTKI," totCols=",totcolsTKI

        totcolsTKI = 0
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,2)
           !allocate(CCmattmp(colsikpq,n_st))
           !do kk = 1,n_st
           !do m = 1,colsikpq
           !   CCmattmp(m,kk) = psi_in(kk,idxs_connectedI_alpha(j)+m-1)
           !enddo
           !enddo
           do m = 1,colsikpq
              !   tmpvar = CCmattmp(m,kk)
              do l = 1,rowsTKI
              do kk = 1,n_st
                 !TKI(kk,l,totcolsTKI+m) = AIJpqContainer(NSOMOalpha,extype,pmodel,qmodel,l,m) * tmpvar
                 !TKI(kk,l,totcolsTKI+m) = AIJpqContainer(l,m,pmodel,qmodel,extype,NSOMOalpha) * tmpvar
                 TKI(kk,l,totcolsTKI+m) = AIJpqContainer(l,m,pmodel,qmodel,extype,NSOMOalpha) * psi_in(kk,idxs_connectedI_alpha(j)+m-1)
                 !TKI(kk,l,totcolsTKI+m) = AIJpqContainer(l,m,pmodel,qmodel,extype,NSOMOalpha) * tmpvar
              enddo
           enddo
           enddo
           !deallocate(CCmattmp)
           do m = 1,colsikpq
              do l = 1,nconnectedI
                 ! <ij|kl> = (ik|jl)
                 moi = excitationIds(1,j) ! p
                 mok = excitationIds(2,j) ! q
                 moj = excitationIds(2,l) ! s
                 mol = excitationIds(1,l) ! r
                 if(moi.EQ.mok .AND. moj.EQ.mol)then
                    diagfac = diagfactors(j)
                    diagfac *= diagfactors(l)
                    !print *,"integrals (",totcolsTKI+m,l,")",mok,moi,mol,moj, "|", diagfac
                    GIJpqrs(totcolsTKI+m,l) = diagfac*0.5d0*mo_two_e_integral(mok,mol,moi,moj) ! g(pq,sr) = <ps,qr>
                 else
                       diagfac = diagfactors(j)*diagfactors(l)
                       !print *,"integrals (",totcolsTKI+m,l,")",mok,moi,mol,moj, "|", diagfac
                       GIJpqrs(totcolsTKI+m,l) = diagfac*0.5d0*mo_two_e_integral(mok,mol,moi,moj) ! g(pq,sr) = <ps,qr>
                    !endif
                 endif
              enddo
           enddo
           totcolsTKI += colsikpq
        enddo


        ! Do big BLAS
        ! TODO TKI, size(TKI,1)*size(TKI,2)
        call dgemm('N','N', rowsTKI*n_st, nconnectedI, totcolsTKI, 1.d0,  &
             TKI, size(TKI,1)*size(TKI,2), GIJpqrs, size(GIJpqrs,1), 0.d0, &
             TKIGIJ , size(TKIGIJ,1)*size(TKIGIJ,2) )


        ! Collect the result
        totcolsTKI = 0
        do j = 1,nconnectedI
           NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
           NSOMOI     = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,2)
           allocate(CCmattmp(n_st,colsikpq))
           CCmattmp = 0.d0

           call dgemm('N','N', n_st, colsikpq, rowsTKI, 1.d0,  &
             TKIGIJ(1:n_st,1:rowsTKI,j), n_st,       &
             AIJpqContainer(1:rowsTKI,1:colsikpq,pmodel,qmodel,extype,NSOMOalpha), rowsTKI, 0.d0, &
             CCmattmp, size(CCmattmp,1) )
           !do kk=1,n_st
           !  do m=1,colsikpq
           !    do l=1,rowsTKI
           !      CCmattmp(kk,m) += TKIGIJ(kk,l,j) * AIJpqContainer(l,m,pmodel,qmodel,extype,NSOMOalpha)
           !    enddo
           !  enddo
           !enddo

           do m = 1,colsikpq
              !do l = 1,rowsTKI
              do kk = 1,n_st
            !call omp_set_lock(lock(idxs_connectedI_alpha(j)+m-1))
                 !psi_out(kk,idxs_connectedI_alpha(j)+m-1) += AIJpqContainer(l,m,pmodel,qmodel,extype,NSOMOalpha) * TKIGIJ(kk,l,j)
                 psi_out(kk,idxs_connectedI_alpha(j)+m-1) += CCmattmp(kk,m)
            !call omp_unset_lock(lock(idxs_connectedI_alpha(j)+m-1))
              !enddo
           enddo
           enddo
           deallocate(CCmattmp)
           totcolsTKI += colsikpq
        enddo

        deallocate(TKI) ! coefficients of CSF
        deallocate(GIJpqrs)  ! gpqrs
        deallocate(TKIGIJ)  ! gpqrs

     enddo ! loop over alphas
  enddo ! loop over I
  !$OMP END DO
  call omp_set_max_active_levels(4)
  deallocate(connectedI_alpha)
  deallocate(idxs_connectedI_alpha)
  deallocate(excitationIds)
  deallocate(excitationTypes)
  deallocate(diagfactors)


  ! Add the diagonal contribution
  !$OMP DO
  do kk=1,n_st
  do i = 1,n_CSF
     psi_out(kk,i) += diag_energies(i)*psi_in(kk,i)
  enddo
  enddo
  !$OMP END DO

  !$OMP end parallel


end subroutine calculate_sigma_vector_cfg_nst_naive_store

subroutine calculate_sigma_vector_cfg_nst(psi_out, psi_in, n_st, sze, istart, iend, ishift, istep)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Documentation for sigma-vector calculation
  !
  ! Calculates the result of the
  ! application of the hamiltonian to the
  ! wavefunction in CFG basis once
  ! TODO : Things prepare outside this routine
  !  1. Touch the providers for
  !     a. ApqIJ containers
  !     b. DET to CSF transformation matrices
  !  2. DET to CSF transcormation
  !  2. CSF to DET back transcormation
  ! returns : psi_coef_out_det :
  END_DOC
  integer,intent(in)             :: sze, istart,iend, istep, ishift, n_st
  real*8,intent(in)              :: psi_in(sze,n_st)
  real*8,intent(out)             :: psi_out(sze,n_st)
  integer(bit_kind)              :: Icfg(N_INT,2)
  integer                        :: i,j,k,l,p,q,noccp,noccq, ii, jj, m, n, idxI, kk, nocck,orbk
  integer(bit_kind),dimension(:,:,:),allocatable :: alphas_Icfg
  integer(bit_kind),dimension(:,:,:),allocatable :: singlesI
  integer(bit_kind),dimension(:,:,:),allocatable :: connectedI_alpha
  integer,dimension(:),allocatable :: idxs_singlesI
  integer,dimension(:),allocatable :: idxs_connectedI_alpha
  integer,dimension(:,:),allocatable :: excitationIds_single
  integer,dimension(:),allocatable :: excitationTypes_single
  integer,dimension(:,:),allocatable :: excitationIds
  integer,dimension(:),allocatable :: excitationTypes
  real*8,dimension(:),allocatable :: diagfactors
  integer                        :: nholes
  integer                        :: nvmos
  integer                        :: listvmos(mo_num)
  integer                        :: vmotype(mo_num) ! 1 -> VMO 2 -> SOMO
  integer                        :: listholes(mo_num)
  integer                        :: holetype(mo_num) ! 1-> SOMO 2->DOMO
  integer                        :: Nalphas_Icfg, nconnectedI, rowsikpq, colsikpq, nsinglesI
  integer                        :: extype,NSOMOalpha,NSOMOI,NSOMOJ,pmodel,qmodel
  integer                        :: getNSOMO
  integer                        :: totcolsTKI
  integer                        :: rowsTKI
  integer                        :: noccpp
  integer                        :: istart_cfg, iend_cfg
  integer*8                      :: MS, Isomo, Idomo, Jsomo, Jdomo, Ialpha, Ibeta
  integer                        :: moi, moj, mok, mol, starti, endi, startj, endj, cnti, cntj, cntk
  real*8                         :: norm_coef_cfg, fac2eints
  real*8                         :: norm_coef_det
  real*8                         :: meCC1, meCC2, diagfac
  real*8,dimension(:,:,:),allocatable :: TKI
  real*8,dimension(:,:),allocatable :: GIJpqrs
  real*8,dimension(:,:,:),allocatable :: TKIGIJ
  real*8, external               :: mo_two_e_integral
  real*8, external               :: get_two_e_integral
  real*8                         :: diag_energies(n_CSF)

  ! allocate
  allocate(alphas_Icfg(N_INT,2,max(sze/2,100)))
  allocate(singlesI(N_INT,2,max(sze/2,100)))
  allocate(connectedI_alpha(N_INT,2,max(sze/2,100)))
  allocate(idxs_singlesI(max(sze/2,100)))
  allocate(idxs_connectedI_alpha(max(sze/2,100)))
  allocate(excitationIds_single(2,max(sze/2,100)))
  allocate(excitationTypes_single(max(sze/2,100)))
  allocate(excitationIds(2,max(sze/2,100)))
  allocate(excitationTypes(max(sze/2,100)))
  allocate(diagfactors(max(sze/2,100)))


  !print *," sze = ",sze
  call calculate_preconditioner_cfg(diag_energies)

  MS = 0
  norm_coef_cfg=0.d0

  psi_out=0.d0

  istart_cfg = psi_csf_to_config_data(istart)
  iend_cfg   = psi_csf_to_config_data(iend)


  !!! Single Excitations !!!
  do i=istart_cfg,iend_cfg
    print *,"I=",i

    ! if Seniority_range > 8 then
    ! continue
    ! else
    ! cycle

    Icfg(1,1) = psi_configuration(1,1,i)
    Icfg(1,2) = psi_configuration(1,2,i)
    starti = psi_config_data(i,1)
    endi   = psi_config_data(i,2)

    ! Returns all unique (checking the past) singly excited cfgs connected to I
    Nalphas_Icfg = 0
    ! TODO:
    ! test if size(alphas_Icfg,1) < Nmo**2) then deallocate + allocate
    !call obtain_associated_alphaI(i, Icfg, alphas_Icfg, Nalphas_Icfg)
    Nalphas_Icfg = NalphaIcfg_list(i)
    alphas_Icfg(1:N_int,1:2,1:Nalphas_Icfg) = alphasIcfg_list(1:n_int,1:2,i,1:Nalphas_Icfg)

    ! TODO : remove doubly excited for return
    ! Here we do 2x the loop. One to count for the size of the matrix, then we compute.
    do k = 1,Nalphas_Icfg
      ! Now generate all singly excited with respect to a given alpha CFG
      call obtain_connected_I_foralpha(i,alphas_Icfg(1,1,k),connectedI_alpha,idxs_connectedI_alpha,nconnectedI,excitationIds,excitationTypes,diagfactors)

      totcolsTKI = 0
      rowsTKI = -1
      do j = 1,nconnectedI
        NSOMOalpha = getNSOMO(alphas_Icfg(1,1,k))
        NSOMOI = getNSOMO(connectedI_alpha(1,1,j))
        p = excitationIds(1,j)
        q = excitationIds(2,j)
        extype = excitationTypes(j)
        call convertOrbIdsToModelSpaceIds(alphas_Icfg(1,1,k), connectedI_alpha(1,1,j), p, q, extype, pmodel, qmodel)
        ! for E_pp E_rs and E_ppE_rr case
        if(p.EQ.q) then
          NSOMOalpha = NSOMOI
        endif
        rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,1)
        colsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,2)
        totcolsTKI += colsikpq
!        if(rowsTKI .LT. rowsikpq .AND. rowsTKI .NE. -1) then
!          print *,">",j,"Something is wrong in sigma-vector", rowsTKI, rowsikpq, "(p,q)=",pmodel,qmodel,"ex=",extype,"na=",NSOMOalpha," nI=",NSOMOI
!          !rowsTKI = rowsikpq
!        else
          rowsTKI = rowsikpq
!        endif
      enddo

      allocate(TKI(n_st,rowsTKI,totcolsTKI)) ! coefficients of CSF
      ! Initialize the inegral container
      ! dims : (totcolsTKI, nconnectedI)
      allocate(GIJpqrs(totcolsTKI,nconnectedI))  ! gpqrs
      allocate(TKIGIJ(n_st,rowsTKI,nconnectedI))  ! TKI * gpqrs

      totcolsTKI = 0
      do j = 1,nconnectedI
        NSOMOalpha = getNSOMO(alphas_Icfg(1,1,k))
        NSOMOI = getNSOMO(connectedI_alpha(1,1,j))
        p = excitationIds(1,j)
        q = excitationIds(2,j)
        extype = excitationTypes(j)
        call convertOrbIdsToModelSpaceIds(alphas_Icfg(1,1,k), connectedI_alpha(1,1,j), p, q, extype, pmodel, qmodel)
        rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,1)
        colsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,2)
        do m = 1,colsikpq
          do l = 1,rowsTKI
            do kk = 1,n_st
              TKI(kk,l,totcolsTKI+m) = AIJpqContainer(l,m,pmodel,qmodel,extype,NSOMOalpha) * psi_in(kk,idxs_connectedI_alpha(j)+m-1)
            enddo
          enddo
        enddo
        do m = 1,colsikpq
          do l = 1,nconnectedI
            ! <ij|kl> = (ik|jl)
            moi = excitationIds(1,j) ! p
            mok = excitationIds(2,j) ! q
            moj = excitationIds(2,l) ! s
            mol = excitationIds(1,l) ! r
            if(moi.EQ.mok .AND. moj.EQ.mol)then
              diagfac = diagfactors(j)
              diagfac *= diagfactors(l)
              !print *,"integrals (",totcolsTKI+m,l,")",mok,moi,mol,moj, "|", diagfac
              GIJpqrs(totcolsTKI+m,l) = diagfac*0.5d0*mo_two_e_integral(mok,mol,moi,moj) ! g(pq,sr) = <ps,qr>
            else
              diagfac = diagfactors(j)*diagfactors(l)
              !print *,"integrals (",totcolsTKI+m,l,")",mok,moi,mol,moj, "|", diagfac
              GIJpqrs(totcolsTKI+m,l) = diagfac*0.5d0*mo_two_e_integral(mok,mol,moi,moj) ! g(pq,sr) = <ps,qr>
              !endif
            endif
          enddo
        enddo
        totcolsTKI += colsikpq
      enddo



      ! Do big BLAS
      ! TODO TKI, size(TKI,1)*size(TKI,2)
      call dgemm('N','N', rowsTKI*n_st, nconnectedI, totcolsTKI, 1.d0,&
          TKI, size(TKI,1)*size(TKI,2), GIJpqrs, size(GIJpqrs,1), 0.d0,&
          TKIGIJ , size(TKIGIJ,1)*size(TKIGIJ,2) )


      ! Collect the result
      totcolsTKI = 0
      do j = 1,nconnectedI
        NSOMOalpha = getNSOMO(alphas_Icfg(1,1,k))
        NSOMOI     = getNSOMO(connectedI_alpha(1,1,j))
        p = excitationIds(1,j)
        q = excitationIds(2,j)
        extype = excitationTypes(j)
        call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
        rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,1)
        colsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,2)
        do m = 1,colsikpq
          do l = 1,rowsTKI
            do kk = 1,n_st
              psi_out(kk,idxs_connectedI_alpha(j)+m-1) = psi_out(kk,idxs_connectedI_alpha(j)+m-1) + &
                AIJpqContainer(l,m,pmodel,qmodel,extype,NSOMOalpha) * TKIGIJ(kk,l,j)
            enddo
          enddo
        enddo
        totcolsTKI += colsikpq
      enddo

      deallocate(TKI) ! coefficients of CSF
      ! Initialize the inegral container
      ! dims : (totcolsTKI, nconnectedI)
      deallocate(GIJpqrs)  ! gpqrs
      deallocate(TKIGIJ)  ! gpqrs

    enddo ! loop over alphas
  enddo ! loop over I
  deallocate(connectedI_alpha)
  deallocate(idxs_connectedI_alpha)
  deallocate(excitationIds)
  deallocate(excitationTypes)
  deallocate(diagfactors)


  ! Add the diagonal contribution
  do i = 1,n_CSF
    do kk=1,n_st
      psi_out(kk,i) += diag_energies(i)*psi_in(kk,i)
    enddo
  enddo
  call omp_set_max_active_levels(4)

end subroutine calculate_sigma_vector_cfg_nst_naive_store

