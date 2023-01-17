 real*8 function logabsgamma(x)
 implicit none
 real*8, intent(in) :: x
 logabsgamma = 1.d32  ! Avoid floating point exception
 if (x>0.d0) then
   logabsgamma = log(abs(gamma(x)))
 endif
 end function logabsgamma
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
  integer MS, ialpha
  MS = elec_alpha_num-elec_beta_num
  NSOMOMax = min(elec_num, cfg_nsomo_max + 2)
  if(AND(cfg_nsomo_min , 1) .eq. 0)then
    NSOMOMin = max(0,cfg_nsomo_min-2)
  else
    NSOMOMin = max(1,cfg_nsomo_min-2)
  endif
  ! Note that here we need NSOMOMax + 2 sizes
  ialpha = (NSOMOMax + MS)/2
  NCSFMax = max(1,nint((binom(NSOMOMax,ialpha)-binom(NSOMOMax,ialpha+1)))) ! TODO: NCSFs for MS=0 (CHECK)
  NBFMax = NCSFMax
  maxDetDimPerBF = max(1,nint((binom(NSOMOMax,ialpha))))
  NMO = n_act_orb
  integer i,j,k,l
  integer startdet,enddet
  integer ncfg,ncfgprev
  integer NSOMO
  integer dimcsfpercfg
  integer detDimperBF
  real*8 :: coeff, binom1, binom2
  integer ncfgpersomo
  real*8, external :: logabsgamma
  detDimperBF = 0
  ! number of cfgs = number of dets for 0 somos
  n_CSF = cfg_seniority_index(NSOMOMin)-1
  ncfgprev = cfg_seniority_index(NSOMOMin)
  !do i = 0-iand(MS,1)+2, NSOMOMax,2
  !!print *," i=",0," dimcsf=",1," ncfg=",ncfgprev, " senor=",cfg_seniority_index(0)
  !!do i = NSOMOMin+2, NSOMOMax,2
  !!   if(cfg_seniority_index(i) .EQ. -1)then
  !!      ncfgpersomo = N_configuration + 1
  !!   else
  !!      ncfgpersomo = cfg_seniority_index(i)
  !!   endif
  !!ncfg = ncfgpersomo - ncfgprev
  !!!detDimperBF = max(1,nint((binom(i,(i+1)/2))))
  !!!dimcsfpercfg = max(1,nint((binom(i-2,(i-2+1)/2)-binom(i-2,((i-2+1)/2)+1))))
  !!n_CSF += ncfg * dimcsfpercfg
  !!!if(cfg_seniority_index(i+2) == -1) EXIT
  !!!if(detDimperBF > maxDetDimPerBF) maxDetDimPerBF = detDimperBF
  !!ncfgprev = cfg_seniority_index(i)
  !!print *," i=",i," dimcsf=",dimcsfpercfg," ncfg=",ncfg, " senor=",cfg_seniority_index(i)
  !!enddo
  !!print *," ^^^^^ N_CSF = ",n_CSF," N_CFG=",N_configuration
  n_CSF = 0
  !ncfgprev = cfg_seniority_index(0)
  !ncfgpersomo = ncfgprev
  !do i = iand(MS,1), NSOMOMax-2,2
  !  if(cfg_seniority_index(i) .EQ. -1) then
  !    cycle
  !  endif
  !  if(cfg_seniority_index(i+2) .EQ. -1) then
  !    ncfgpersomo = N_configuration + 1
  !  else
  !    if(cfg_seniority_index(i+2) > ncfgpersomo) then
  !        ncfgpersomo = cfg_seniority_index(i+2)
  !    else
  !      k = 0
  !      do while(cfg_seniority_index(i+2+k) < ncfgpersomo)
  !        k = k + 2
  !        ncfgpersomo = cfg_seniority_index(i+2+k)
  !      enddo
  !    endif
  !  endif
  !  ncfg = ncfgpersomo - ncfgprev
  !  if(i .EQ. 0 .OR. i .EQ. 1) then
  !    dimcsfpercfg = 1
  !  elseif( i .EQ. 3) then
  !    dimcsfpercfg = 2
  !  else
  !    if(iand(MS,1) .EQ. 0) then
  !      dimcsfpercfg = max(1,nint((binom(i,i/2)-binom(i,i/2+1))))
  !    else
  !      dimcsfpercfg = max(1,nint((binom(i,(i+1)/2)-binom(i,(i+3)/2))))
  !    endif
  !  endif
  !  n_CSF += ncfg * dimcsfpercfg
  !  print *," i=",i," dimcsf=",dimcsfpercfg," ncfg=",ncfg, " senor=",cfg_seniority_index(i)
  !  if(cfg_seniority_index(i+2) > ncfgprev) then
  !    ncfgprev = cfg_seniority_index(i+2)
  !  else
  !    k = 0
  !    do while(cfg_seniority_index(i+2+k) < ncfgprev)
  !      k = k + 2
  !      ncfgprev = cfg_seniority_index(i+2+k)
  !    enddo
  !  endif
  !enddo
  n_CSF = 0
  !print *," -9(((((((((((((( NSOMOMin=",NSOMOMin
  ncfgprev = cfg_seniority_index(NSOMOMin) ! can be -1 
  if(ncfgprev.eq.-1)then
    ncfgprev=1
  endif
  do i=NSOMOMin,NSOMOMax+2,2
    !k=0
    !do while((cfg_seniority_index(i+2+k) .eq. -1) .and. (k.le.NSOMOMax))
    !  k=k+2
    !end do
    if(cfg_seniority_index(i).eq.-1)cycle
    if(cfg_seniority_index(i+2).eq.-1)then
      ncfg = N_configuration - ncfgprev + 1
      if(ncfg .eq. 0)then
        ncfg=1
      endif
    else
      ncfg = cfg_seniority_index(i+2) - ncfgprev
    endif
    if(i .EQ. 0 .OR. i .EQ. 1) then
      dimcsfpercfg = 1
    elseif( i .EQ. 3) then
      dimcsfpercfg = 2
    else
      if(iand(MS,1) .EQ. 0) then
        ialpha = (i + MS)/2
        dimcsfpercfg = max(1,nint((binom(i,ialpha)-binom(i,ialpha+1))))
      else
        ialpha = (i + MS)/2
        dimcsfpercfg = max(1,nint((binom(i,ialpha)-binom(i,ialpha+1))))
      endif
    endif
    n_CSF += ncfg*dimcsfpercfg
    !print *," i=",i," dimcsf=",dimcsfpercfg," ncfg=",ncfg, " ncfgprev=",ncfgprev, " senor=",cfg_seniority_index(i)
    ncfgprev = cfg_seniority_index(i+2)
  end do
  !print *," ^^^^^ N_CSF = ",n_CSF," N_CFG=",N_configuration
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

  ! Initialize deta and detb
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
  integer   :: startdet, enddet, idx
  integer*8 MS, salpha
  integer ndetI
  integer :: getNSOMO
  real*8,dimension(:,:),allocatable    :: tempBuffer
  real*8,dimension(:),allocatable    :: tempCoeff
  real*8  :: norm_det1, phasedet

  integer :: nt


  norm_det1 = 0.d0
  MS = elec_alpha_num - elec_beta_num
  ! initialization
  psi_coef_config = 0.d0
  DetToCSFTransformationMatrix(0,:,:) = 1.d0
  do i = 2-iand(MS,1_8), NSOMOMax,2
    Isomo = IBSET(0_8, i) - 1_8
    ! rows = Ncsfs
    ! cols = Ndets
    salpha = (i+MS)/2
    bfIcfg = max(1,nint((binom(i,salpha)-binom(i,salpha+1))))
    ndetI =  max(1,nint((binom(i,salpha))))
    !bfIcfg = max(1,nint((binom(i,(i+1)/2)-binom(i,((i+1)/2)+1))))
    !ndetI =  max(1,nint((binom(i,(i+1)/2))))

    allocate(tempBuffer(bfIcfg,ndetI))
    call getCSFtoDETTransformationMatrix(Isomo, MS, NBFMax, maxDetDimPerBF, tempBuffer)
    DetToCSFTransformationMatrix(i,1:bfIcfg,1:ndetI) =  tempBuffer(1:bfIcfg,1:ndetI)
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
         idx = psi_configuration_to_psi_det_data(j)
         Ialpha(:) = psi_det(:,1,idx)
         Ibeta(:)  = psi_det(:,2,idx)
         call get_phase_qp_to_cfg(Ialpha, Ibeta, phasedet)
         tempCoeff(countdet) = psi_coef(idx, istate)*phasedet
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
      salpha = (s+MS)/2
      bfIcfg = max(1,nint((binom(s,salpha)-binom(s,salpha+1))))
      !bfIcfg = max(1,nint((binom(s,(s+1)/2)-binom(s,((s+1)/2)+1))))

      ! perhaps blocking with CFGs of same seniority
      ! can be more efficient
      allocate(tempBuffer(bfIcfg,ndetI))
      tempBuffer = DetToCSFTransformationMatrix(s,:bfIcfg,:ndetI)

       call dgemm('N','N', bfIcfg, 1, ndetI, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, size(tempCoeff,1), 0.d0, psi_coef_config(countcsf+1,1), size(psi_coef_config,1))
       !call dgemv('N', NBFMax, maxDetDimPerBF, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, 1, 0.d0, psi_coef_config(countcsf), 1)

      deallocate(tempCoeff)
      deallocate(tempBuffer)
      psi_config_data(i,1) = countcsf + 1
      do k=1,bfIcfg
        psi_csf_to_config_data(countcsf+k) = i
      enddo
      countcsf += bfIcfg
      psi_config_data(i,2) = countcsf
  enddo
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
  !print *,"NSOMOMax = ",NSOMOMax
  !print *,"NSOMOMin = ",NSOMOMin
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
  !print *,"Rowsmax=",rowsmax," Colsmax=",colsmax
  END_PROVIDER

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
  MS = elec_alpha_num-elec_beta_num
  real*8,dimension(:,:),allocatable :: meMatrix
  integer maxdim

  ! Type
  ! 1. SOMO -> SOMO
  AIJpqContainer = 0.d0
  AIJpqContainer(1,1,1,1,1,NSOMOMin) = 1.0d0
  integer :: rows_old, cols_old
  rows_old = -1
  cols_old = -1
  allocate(meMatrix(1,1))
  do i = NSOMOMin+2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     j=i-2
     if(j .GT. NSOMOMax .OR. j .LT. 0) cycle
     nsomoi = i
     do k = 1,i
        orbp = k
        do l = 1,i

           ! Define Jsomo
           if(k .NE. l) then
              Jsomo = IBCLR(Isomo, k-1)
              Jsomo = IBCLR(Jsomo, l-1)
              nsomoj = j
           else
              Isomo = ISHFT(1_8,i)-1
              Jsomo = ISHFT(1_8,i)-1
              nsomoj = i
           endif

           call getApqIJMatrixDims(Isomo,           &
                Jsomo, &
                MS,                       &
                rows,                     &
                cols)

           orbq = l
           if ((rows /= rows_old).or.(cols /= cols_old)) then
             deallocate(meMatrix)
             allocate(meMatrix(rows,cols))
             rows_old = rows
             cols_old = cols
           endif
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
           ! i -> j
          do ri = 1,rows
              do ci = 1,cols
                 AIJpqContainer(ri,ci,k,l,1,nsomoi) = meMatrix(ri, ci)
              end do
           end do
        end do
     end do
  end do
  deallocate(meMatrix)

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
  AIJpqContainer(1,1,1:2,1:2,3,NSOMOMin) = 1.0d0
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
  AIJpqContainer(1,1,2,2,4,NSOMOMin) = 1.0d0
  AIJpqContainer(1,1,2,1,4,NSOMOMin) =-1.0d0
  AIJpqContainer(1,1,1,2,4,NSOMOMin) =-1.0d0
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

              !print *,"k,l=",k,l
              !call debug_spindet(Jsomo,1)
              !call debug_spindet(Isomo,1)

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
  integer :: i,j,k,kk,l,p,q,noccp,noccq, ii, jj, iii
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
  real*8            :: core_act_contrib
  integer                        :: listall(N_int*bit_kind_size), nelall

  !PROVIDE h_core_ri
  PROVIDE core_fock_operator
  PROVIDE h_act_ri
  ! initialize energies
  diag_energies = 0.d0
  !print *,"Core energy=",core_energy," nucler rep=",nuclear_repulsion, " n_core_orb=",n_core_orb," n_act_orb=",n_act_orb," mo_num=",mo_num

  ! calculate core energy
  diag_energies = core_energy - nuclear_repulsion

  ! calculate the core energy
  !print *,"Core 2energy=",ref_bitmask_energy

  do i=1,N_configuration

     !Isomo = psi_configuration(1,1,i)
     !Idomo = psi_configuration(1,2,i)
     !Icfg(1,1) = psi_configuration(1,1,i)
     !Icfg(1,2) = psi_configuration(1,2,i)
     !NSOMOI = getNSOMO(psi_configuration(:,:,i))

     starti = psi_config_data(i,1)
     endi   = psi_config_data(i,2)

     core_act_contrib = 0.0d0

     ! find out all pq holes possible
     nholes = 0
     listholes = -1
     ! holes in SOMO
     !do kk = 1,n_act_orb
     !  k = list_act(kk)
     !   if(POPCNT(IAND(Isomo,IBSET(0_8,k-1))) .EQ. 1) then
     !      nholes += 1
     !      listholes(nholes) = k
     !      holetype(nholes) = 1
     !   endif
     !enddo
     call bitstring_to_list(psi_configuration(1,1,i),listall,nelall,N_int)

     do iii=1,nelall
       nholes += 1
       listholes(nholes) = listall(iii)
       holetype(nholes) = 1
     end do

     ! holes in DOMO
     !do kk = 1,n_act_orb
     !  k = list_act(kk)
     !   if(POPCNT(IAND(Idomo,IBSET(0_8,k-1))) .EQ. 1) then
     !      nholes += 1
     !      listholes(nholes) = k
     !      holetype(nholes) = 2
     !   endif
     !enddo
     call bitstring_to_list(psi_configuration(1,2,i),listall,nelall,N_int)

     do iii=1,nelall
       if(listall(iii) .gt. n_core_orb)then
         nholes += 1
         listholes(nholes) = listall(iii)
         holetype(nholes) = 2
       endif
     end do


     !!! find vmos
     !!listvmos = -1
     !!vmotype = -1
     !!nvmos = 0
     !!!do k = n_core_orb+1,n_core_orb + n_act_orb
     !!!do k = 1,mo_num
     !!do kk = 1,n_act_orb
     !!  k = list_act(kk)
     !!   !print *,i,IBSET(0,i-1),POPCNT(IAND(Isomo,(IBSET(0_8,i-1)))), POPCNT(IAND(Idomo,(IBSET(0_8,i-1))))
     !!   if(POPCNT(IAND(Isomo,(IBSET(0_8,k-1)))) .EQ. 0 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,k-1)))) .EQ. 0) then
     !!      nvmos += 1
     !!      listvmos(nvmos) = k
     !!      vmotype(nvmos) = 0
     !!   else if(POPCNT(IAND(Isomo,(IBSET(0_8,k-1)))) .EQ. 1 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,k-1)))) .EQ. 0 ) then
     !!      nvmos += 1
     !!      listvmos(nvmos) = k
     !!      vmotype(nvmos) = 1
     !!   end if
     !!enddo
     !print *,"I=",i
     !call debug_spindet(psi_configuration(1,1,i),N_int)
     !call debug_spindet(psi_configuration(1,2,i),N_int)

     do k=1,nholes
        p = listholes(k)
        noccp = holetype(k)


        ! core-virtual
        do l = 1, n_core_orb
         jj = list_core(l)
         core_act_contrib += noccp * (2.d0 * mo_two_e_integrals_jj(jj,p) - mo_two_e_integrals_jj_exchange(jj,p))
        enddo

        ! Calculate one-electron
        ! and two-electron coulomb terms
        do l=1,nholes
           q = listholes(l)
           noccq = holetype(l)
           !print *,"--------------- K=",p," L=",q

           ! one-electron term
           if(p.EQ.q) then
              hpp = noccq * h_act_ri(p,q)!mo_one_e_integrals(q,q)
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
     !print *,"I=",i," core_act=",core_act_contrib
     do j=starti,endi
       diag_energies(j) += core_act_contrib
     end do
  enddo

end subroutine calculate_preconditioner_cfg

subroutine obtain_connected_I_foralpha_fromfilterdlist(idxI, nconnectedJ, idslistconnectedJ, listconnectedJ, Ialpha, connectedI, idxs_connectedI, nconnectedI, excitationIds, excitationTypes, diagfactors)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Documentation for obtain_connected_I_foralpha
  ! This function returns all those selected configurations
  ! which are connected to the input configuration
  ! Ialpha by a single excitation.
  !
  ! The type of excitations are ordered as follows:
  ! Type 1 - SOMO -> SOMO
  ! Type 2 - DOMO -> VMO
  ! Type 3 - SOMO -> VMO
  ! Type 4 - DOMO -> SOMO
  !
  ! Order of operators
  ! \alpha> = a^\dag_p a_q |I> = E_pq |I>
  END_DOC
  integer          ,intent(in)             :: idxI
  integer          ,intent(in)             :: nconnectedJ
  integer(bit_kind),intent(in)             :: listconnectedJ(N_int,2,*)
  integer(bit_kind),intent(in)             :: Ialpha(N_int,2)
  integer(bit_kind),intent(out)            :: connectedI(N_int,2,*)
  integer          ,intent(in)             :: idslistconnectedJ(*)
  integer          ,intent(out)            :: idxs_connectedI(*)
  integer,intent(out)                      :: nconnectedI
  integer,intent(out)                      :: excitationIds(2,*)
  integer,intent(out)                      :: excitationTypes(*)
  real*8 ,intent(out)                      :: diagfactors(*)
  integer*8                                :: Idomo
  integer*8                                :: Isomo
  integer*8                                :: Jdomo
  integer*8                                :: Jsomo
  integer*8                                :: IJsomo
  integer*8                                :: diffSOMO
  integer*8                                :: diffDOMO
  integer*8                                :: xordiffSOMODOMO
  integer                                  :: ndiffSOMO
  integer                                  :: ndiffDOMO
  integer                                  :: nxordiffSOMODOMO
  integer :: ii,i,j,k,kk,l,p,q,nsomoJ,nsomoalpha,starti,endi,extyp,nholes, idxJ
  integer :: listholes(mo_num)
  integer :: holetype(mo_num)
  integer :: end_index
  integer :: Nsomo_alpha
  logical :: isOKlistJ

  PROVIDE DetToCSFTransformationMatrix

  isOKlistJ = .False.

  nconnectedI = 0
  end_index = N_configuration

  ! Since CFGs are sorted wrt to seniority
  ! we don't have to search the full CFG list
  Isomo = Ialpha(1,1)
  Idomo = Ialpha(1,2)
  Nsomo_alpha = POPCNT(Isomo)
  end_index = min(N_configuration,cfg_seniority_index(min(Nsomo_alpha+4,elec_num))-1)
  if(end_index .LT. 0) end_index= N_configuration
  !end_index = N_configuration


  p = 0
  q = 0
  do i=1,nconnectedJ
     idxJ = idslistconnectedJ(i)
     Isomo = Ialpha(1,1)
     Idomo = Ialpha(1,2)
     Jsomo = listconnectedJ(1,1,i)
     Jdomo = listconnectedJ(1,2,i)
     diffSOMO = IEOR(Isomo,Jsomo)
     ndiffSOMO = POPCNT(diffSOMO)
     diffDOMO = IEOR(Idomo,Jdomo)
     xordiffSOMODOMO = IEOR(diffSOMO,diffDOMO)
     ndiffDOMO = POPCNT(diffDOMO)
     nxordiffSOMODOMO = POPCNT(xordiffSOMODOMO)
     nxordiffSOMODOMO += ndiffSOMO + ndiffDOMO 
     if((nxordiffSOMODOMO .EQ. 4) .AND. ndiffSOMO .EQ. 2) then
        select case(ndiffDOMO)
        case (0)
           ! SOMO -> VMO
           !print *,"obt SOMO -> VMO"
           extyp = 3
           IJsomo = IEOR(Isomo, Jsomo)
IRP_IF WITHOUT_TRAILZ
           p = (popcnt(ieor( IAND(Isomo,IJsomo), IAND(Isomo,IJsomo)-1)) -1) + 1
IRP_ELSE
           p = TRAILZ(IAND(Isomo,IJsomo)) + 1
IRP_ENDIF
           IJsomo = IBCLR(IJsomo,p-1)
IRP_IF WITHOUT_TRAILZ
           q = (popcnt(ieor(IJsomo,IJsomo-1))-1) + 1
IRP_ELSE
           q = TRAILZ(IJsomo) + 1
IRP_ENDIF
        case (1)
           ! DOMO -> VMO
           ! or
           ! SOMO -> SOMO
           nsomoJ = POPCNT(Jsomo)
           nsomoalpha = POPCNT(Isomo)
           if(nsomoJ .GT. nsomoalpha) then
              ! DOMO -> VMO
              !print *,"obt DOMO -> VMO"
              extyp = 2
IRP_IF WITHOUT_TRAILZ
              p = (popcnt(ieor( IEOR(Idomo,Jdomo), IEOR(Idomo,Jdomo)-1))-1) + 1
IRP_ELSE
              p = TRAILZ(IEOR(Idomo,Jdomo)) + 1
IRP_ENDIF
              Isomo = IEOR(Isomo, Jsomo)
              Isomo = IBCLR(Isomo,p-1)
IRP_IF WITHOUT_TRAILZ
              q = (popcnt(ieor(Isomo,Isomo-1))-1) + 1
IRP_ELSE
              q = TRAILZ(Isomo) + 1
IRP_ENDIF
           else
              ! SOMO -> SOMO
              !print *,"obt SOMO -> SOMO"
              extyp = 1
IRP_IF WITHOUT_TRAILZ
              q = (popcnt(ieor( IEOR(Idomo,Jdomo), IEOR(Idomo,Jdomo)-1))-1) + 1
IRP_ELSE
              q = TRAILZ(IEOR(Idomo,Jdomo)) + 1
IRP_ENDIF
              Isomo = IEOR(Isomo, Jsomo)
              Isomo = IBCLR(Isomo,q-1)
IRP_IF WITHOUT_TRAILZ
              p = (popcnt(ieor(Isomo,Isomo-1))-1) + 1
IRP_ELSE
              p = TRAILZ(Isomo) + 1
IRP_ENDIF
           end if
        case (2)
           ! DOMO -> SOMO
           !print *,"obt DOMO -> SOMO"
           extyp = 4
           IJsomo = IEOR(Isomo, Jsomo)
IRP_IF WITHOUT_TRAILZ
           p = (popcnt(ieor(IAND(Jsomo,IJsomo) ,IAND(Jsomo,IJsomo) -1))-1) + 1
IRP_ELSE
           p = TRAILZ(IAND(Jsomo,IJsomo)) + 1
IRP_ENDIF
           IJsomo = IBCLR(IJsomo,p-1)
IRP_IF WITHOUT_TRAILZ
           q = (popcnt(ieor(IJsomo,IJsomo-1))-1) + 1
IRP_ELSE
           q = TRAILZ(IJsomo) + 1
IRP_ENDIF
        case default
           print *,"something went wront in get connectedI"
        end select
        starti = psi_config_data(idxJ,1)
        endi   = psi_config_data(idxJ,2)
        nconnectedI += 1
        connectedI(:,:,nconnectedI) = listconnectedJ(:,:,i)
        idxs_connectedI(nconnectedI)=starti
        excitationIds(1,nconnectedI)=p
        excitationIds(2,nconnectedI)=q
        excitationTypes(nconnectedI) = extyp
        diagfactors(nconnectedI) = 1.0d0
     else if((ndiffSOMO + ndiffDOMO) .EQ. 0) then
        ! find out all pq holes possible
        nholes = 0
        ! holes in SOMO
        Isomo = listconnectedJ(1,1,i)
        Idomo = listconnectedJ(1,2,i)
        do ii = 1,mo_num
           if(POPCNT(IAND(Isomo,IBSET(0_8,ii-1))) .EQ. 1) then
              nholes += 1
              listholes(nholes) = ii
              holetype(nholes) = 1
           endif
        end do
        ! holes in DOMO
        do ii = 1,mo_num
           if(POPCNT(IAND(Idomo,IBSET(0_8,ii-1))) .EQ. 1) then
              nholes += 1
              listholes(nholes) = ii
              holetype(nholes) = 2
           endif
        end do

        do k=1,nholes
           p = listholes(k)
           q = p
           extyp = 1
           if(holetype(k) .EQ. 1) then
              starti = psi_config_data(idxJ,1)
              endi   = psi_config_data(idxJ,2)
              nconnectedI += 1
              connectedI(:,:,nconnectedI) = listconnectedJ(:,:,i)
              idxs_connectedI(nconnectedI)=starti
              excitationIds(1,nconnectedI)=p
              excitationIds(2,nconnectedI)=q
              excitationTypes(nconnectedI) = extyp
              diagfactors(nconnectedI) = 1.0d0
           else
              starti = psi_config_data(idxJ,1)
              endi   = psi_config_data(idxJ,2)
              nconnectedI += 1
              connectedI(:,:,nconnectedI) = listconnectedJ(:,:,i)
              idxs_connectedI(nconnectedI)=starti
              excitationIds(1,nconnectedI)=p
              excitationIds(2,nconnectedI)=q
              excitationTypes(nconnectedI) = extyp
              diagfactors(nconnectedI) = 2.0d0
           endif
        enddo
     endif
  end do

end subroutine obtain_connected_I_foralpha_fromfilterdlist


subroutine convertOrbIdsToModelSpaceIds(Ialpha, Jcfg, p, q, extype, pmodel, qmodel)
  implicit none
  BEGIN_DOC
  ! This function converts the orbital ids
  ! in real space to those used in model space
  ! in order to identify the matrices required
  ! for the calculation of MEs.
  !
  ! The type of excitations are ordered as follows:
  ! Type 1 - SOMO -> SOMO
  ! Type 2 - DOMO -> VMO
  ! Type 3 - SOMO -> VMO
  ! Type 4 - DOMO -> SOMO
  END_DOC
  integer(bit_kind),intent(in)   :: Ialpha(N_int,2)
  integer(bit_kind),intent(in)   :: Jcfg(N_int,2)
  integer,intent(in)             :: p,q
  integer,intent(in)             :: extype
  integer,intent(out)            :: pmodel,qmodel
  integer(bit_kind)              :: Isomo(N_int)
  integer(bit_kind)              :: Idomo(N_int)
  integer(bit_kind)              :: Jsomo(N_int)
  integer(bit_kind)              :: Jdomo(N_int)
  !integer*8                       :: Isomo
  !integer*8                       :: Idomo
  !integer*8                       :: Jsomo
  !integer*8                       :: Jdomo
  integer*8                      :: mask
  integer                        :: iint, ipos, ii
  !integer(bit_kind)              :: Isomotmp(N_int)
  !integer(bit_kind)              :: Jsomotmp(N_int)
  integer*8             :: Isomotmp
  integer*8             :: Jsomotmp
  integer                        :: pos0,pos0prev
  integer                        :: tmpp, tmpq

  ! TODO Flag (print) when model space indices is > 64
  do ii=1,N_int
    Isomo(ii) = Ialpha(ii,1)
    Idomo(ii) = Ialpha(ii,2)
    Jsomo(ii) = Jcfg(ii,1)
    Jdomo(ii) = Jcfg(ii,2)
  end do
  pos0prev = 0
  pmodel = p
  qmodel = q

  if(p .EQ. q) then
     pmodel = 1
     qmodel = 1
  else
     select case(extype)
       case (1)
          ! SOMO -> SOMO
          ! remove all domos
          !print *,"type -> SOMO -> SOMO"
          !mask = ISHFT(1_8,p) - 1
          !Isomotmp = IAND(Isomo,mask)
          !pmodel = POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
          !mask = ISHFT(1_8,q) - 1
          !Isomotmp = IAND(Isomo,mask)
          !qmodel = POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))

          iint = shiftr(p-1,bit_kind_shift) + 1
          ipos = p-shiftl((iint-1),bit_kind_shift)-1
          tmpp = 0
          !print *,"iint=",iint, " p=",p
          do ii=1,iint-1
            !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
            !Isomotmp = IAND(Isomo(ii),mask)
            !tmpp += POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
            tmpp += POPCNT(Isomo(ii))
          end do
          mask = ISHFT(1_bit_kind,ipos+1) - 1
          Isomotmp = IAND(Isomo(iint),mask)
          !pmodel = tmpp + POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
          pmodel = tmpp + POPCNT(Isomotmp)
          !print *,"iint=",iint, " ipos=",ipos,"pmodel=",pmodel, XOR(Isomotmp,mask),Isomo(iint)

          iint = shiftr(q-1,bit_kind_shift) + 1
          ipos = q-shiftl((iint-1),bit_kind_shift)-1
          tmpq = 0
          do ii=1,iint-1
            !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
            !Isomotmp = IAND(Isomo(ii),mask)
            !tmpq += POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
            tmpq += POPCNT(Isomo(ii))
          end do
          mask = ISHFT(1_bit_kind,ipos+1) - 1
          Isomotmp = IAND(Isomo(iint),mask)
          !qmodel = tmpq + POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
          qmodel = tmpq + POPCNT(Isomotmp)
          !print *,"iint=",iint, " ipos=",ipos,"qmodel=",qmodel
       case (2)
          ! DOMO -> VMO
          ! remove all domos except one at p
          !print *,"type -> DOMO -> VMO"
          !mask = ISHFT(1_8,p) - 1
          !Jsomotmp = IAND(Jsomo,mask)
          !pmodel = POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
          !mask = ISHFT(1_8,q) - 1
          !Jsomotmp = IAND(Jsomo,mask)
          !qmodel = POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))

          iint = shiftr(p-1,bit_kind_shift) + 1
          ipos = p-shiftl((iint-1),bit_kind_shift)-1
          tmpp = 0
          do ii=1,iint-1
            !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
            !Jsomotmp = IAND(Jsomo(ii),mask)
            !tmpp += POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
            tmpp += POPCNT(Jsomo(ii))
          end do
          mask = ISHFT(1_bit_kind,ipos+1) - 1
          Jsomotmp = IAND(Jsomo(iint),mask)
          !pmodel = tmpp + POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
          pmodel = tmpp + POPCNT(Jsomotmp)

          iint = shiftr(q-1,bit_kind_shift) + 1
          ipos = q-shiftl((iint-1),bit_kind_shift)-1
          tmpq = 0
          do ii=1,iint-1
            !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
            !Jsomotmp = IAND(Jsomo(ii),mask)
            !tmpq += POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
            tmpq += POPCNT(Jsomo(ii))
          end do
          mask = ISHFT(1_bit_kind,ipos+1) - 1
          Jsomotmp = IAND(Jsomo(iint),mask)
          !qmodel = tmpq + POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
          qmodel = tmpq + POPCNT(Jsomotmp)
       case (3)
          ! SOMO -> VMO
          !print *,"type -> SOMO -> VMO"
          !Isomo = IEOR(Isomo,Jsomo)
          if(p.LT.q) then
             !mask = ISHFT(1_8,p) - 1
             !Isomo = IAND(Isomo,mask)
             !pmodel = POPCNT(mask) - POPCNT(XOR(Isomo,mask))
             !mask = ISHFT(1_8,q) - 1
             !Jsomo = IAND(Jsomo,mask)
             !qmodel = POPCNT(mask) - POPCNT(XOR(Jsomo,mask)) + 1

             iint = shiftr(p-1,bit_kind_shift) + 1
             ipos = p-shiftl((iint-1),bit_kind_shift)-1
             tmpp = 0
             do ii=1,iint-1
               !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
               !Isomotmp = IAND(Isomo(ii),mask)
               !tmpp += POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
               tmpp += POPCNT(Isomo(ii))
             end do
             mask = ISHFT(1_bit_kind,ipos+1) - 1
             Isomotmp = IAND(Isomo(iint),mask)
             !pmodel = tmpp + POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
             pmodel = tmpp + POPCNT(Isomotmp)

             iint = shiftr(q-1,bit_kind_shift) + 1
             ipos = q-shiftl((iint-1),bit_kind_shift)-1
             tmpq = 0
             do ii=1,iint-1
               !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
               !Jsomotmp = IAND(Jsomo(ii),mask)
               !tmpq += POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
               tmpq += POPCNT(Jsomo(ii))
             end do
             mask = ISHFT(1_bit_kind,ipos+1) - 1
             Jsomotmp = IAND(Jsomo(iint),mask)
             !qmodel = tmpq + POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask)) + 1
             qmodel = tmpq + POPCNT(Jsomotmp) + 1
          else
             !mask = ISHFT(1_8,p) - 1
             !Isomo = IAND(Isomo,mask)
             !pmodel = POPCNT(mask) - POPCNT(XOR(Isomo,mask)) + 1
             !mask = ISHFT(1_8,q) - 1
             !Jsomo = IAND(Jsomo,mask)
             !qmodel = POPCNT(mask) - POPCNT(XOR(Jsomo,mask))

             iint = shiftr(p-1,bit_kind_shift) + 1
             ipos = p-shiftl((iint-1),bit_kind_shift)-1
             tmpp = 0
             do ii=1,iint-1
               !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
               !Isomotmp = IAND(Isomo(ii),mask)
               !tmpp += POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
               tmpp += POPCNT(Isomo(ii))
             end do
             mask = ISHFT(1_bit_kind,ipos+1) - 1
             Isomotmp = IAND(Isomo(iint),mask)
             !pmodel = tmpp + POPCNT(mask) - POPCNT(XOR(Isomotmp,mask)) + 1
             pmodel = tmpp + POPCNT(Isomotmp) + 1

             iint = shiftr(q-1,bit_kind_shift) + 1
             ipos = q-shiftl((iint-1),bit_kind_shift)-1
             tmpq = 0
             do ii=1,iint-1
               !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
               !Jsomotmp = IAND(Jsomo(ii),mask)
               !tmpq += POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
               tmpq += POPCNT(Jsomo(ii))
             end do
             mask = ISHFT(1_bit_kind,ipos+1) - 1
             Jsomotmp = IAND(Jsomo(iint),mask)
             !qmodel = tmpq + POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
             qmodel = tmpq + POPCNT(Jsomotmp)
          endif
       case (4)
          ! DOMO -> SOMO
          ! remove all domos except one at p
          !print *,"type -> DOMO -> SOMO"
          !Isomo = IEOR(Isomo,Jsomo)
          if(p.LT.q) then
             !mask = ISHFT(1_8,p) - 1
             !Jsomo = IAND(Jsomo,mask)
             !pmodel = POPCNT(mask) - POPCNT(XOR(Jsomo,mask))
             !mask = ISHFT(1_8,q) - 1
             !Isomo = IAND(Isomo,mask)
             !qmodel = POPCNT(mask) - POPCNT(XOR(Isomo,mask)) + 1

             iint = shiftr(p-1,bit_kind_shift) + 1
             ipos = p-shiftl((iint-1),bit_kind_shift)-1
             tmpp = 0
             do ii=1,iint-1
               !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
               !Jsomotmp = IAND(Jsomo(ii),mask)
               !tmpp += POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
               tmpp += POPCNT(Jsomo(ii))
             end do
             mask = ISHFT(1_bit_kind,ipos+1) - 1
             Jsomotmp = IAND(Jsomo(iint),mask)
             !pmodel = tmpp + POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
             pmodel = tmpp + POPCNT(Jsomotmp)

             iint = shiftr(q-1,bit_kind_shift) + 1
             ipos = q-shiftl((iint-1),bit_kind_shift)-1
             tmpq = 0
             do ii=1,iint-1
               !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
               !Isomotmp = IAND(Isomo(ii),mask)
               !tmpq += POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
               tmpq += POPCNT(Isomo(ii))
             end do
             mask = ISHFT(1_bit_kind,ipos+1) - 1
             Isomotmp = IAND(Isomo(iint),mask)
             !qmodel = tmpq + POPCNT(mask) - POPCNT(XOR(Isomotmp,mask)) + 1
             qmodel = tmpq + POPCNT(Isomotmp) + 1
          else
             !mask = ISHFT(1_8,p) - 1
             !Jsomo = IAND(Jsomo,mask)
             !pmodel = POPCNT(mask) - POPCNT(XOR(Jsomo,mask)) + 1
             !mask = ISHFT(1_8,q) - 1
             !Isomo = IAND(Isomo,mask)
             !qmodel = POPCNT(mask) - POPCNT(XOR(Isomo,mask))

             iint = shiftr(p-1,bit_kind_shift) + 1
             ipos = p-shiftl((iint-1),bit_kind_shift)-1
             tmpp = 0
             do ii=1,iint-1
               !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
               !Jsomotmp = IAND(Jsomo(ii),mask)
               !tmpp += POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask))
               tmpp += POPCNT(Jsomo(ii))
             end do
             mask = ISHFT(1_bit_kind,ipos+1) - 1
             Jsomotmp = IAND(Jsomo(iint),mask)
             !pmodel = tmpp + POPCNT(mask) - POPCNT(XOR(Jsomotmp,mask)) + 1
             pmodel = tmpp + POPCNT(Jsomotmp) + 1

             iint = shiftr(q-1,bit_kind_shift) + 1
             ipos = q-shiftl((iint-1),bit_kind_shift)-1
             tmpq = 0
             do ii=1,iint-1
               !mask = ISHFT(1_bit_kind,-1)-1_bit_kind
               !Isomotmp = IAND(Isomo(ii),mask)
               !tmpq += POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
               tmpq += POPCNT(Isomo(ii))
             end do
             mask = ISHFT(1_bit_kind,ipos+1) - 1
             Isomotmp = IAND(Isomo(iint),mask)
             !qmodel = tmpq + POPCNT(mask) - POPCNT(XOR(Isomotmp,mask))
             qmodel = tmpq + POPCNT(Isomotmp)
          endif
       case default
          print *,"something is wrong in convertOrbIdsToModelSpaceIds"
     end select
  endif
  !print *,p,q,"model ids=",pmodel,qmodel
end subroutine convertOrbIdsToModelSpaceIds

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
  integer                        :: i,j,k,l,p,q,noccp,noccq, m, n, idxI, nocck,orbk
  integer                        :: ii,jj,kk,ll,pp,qq
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
  integer                        :: istart_cfg, iend_cfg, num_threads_max
  integer                        :: iint, jint, ipos, jpos, Nsomo_I, iii
  integer                        :: nconnectedJ,nconnectedtotalmax,nconnectedmaxJ,maxnalphas,ntotJ
  integer*8                      :: MS,Ialpha, Ibeta
  integer(bit_kind)              :: Isomo(N_INT)
  integer(bit_kind)              :: Idomo(N_INT)
  integer(bit_kind)              :: Jsomo(N_INT)
  integer(bit_kind)              :: Jdomo(N_INT)
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
  real*8,dimension(:),allocatable:: diag_energies
  real*8                         :: tmpvar, tmptot
  real*8                         :: core_act_contrib
  integer :: listall(N_int*bit_kind_size), nelall
  integer :: countelec

  integer(omp_lock_kind), allocatable :: lock(:)
  call omp_set_max_active_levels(1)

  !print *," sze = ",sze
  allocate(lock(sze))
  do i=1,sze
    call omp_init_lock(lock(i))
  enddo
  !do i=1,size(psi_config_data,1)
  !  print *,"i=",i," psi_cfg_data_1=",psi_config_data(i,1)," psi_cfg_data_2=",psi_config_data(i,2)
  !end do

  allocate(diag_energies(n_CSF))
  call calculate_preconditioner_cfg(diag_energies)
  !print *," diag energy =",diag_energies(1)

  MS = 0
  norm_coef_cfg=0.d0

  psi_out=0.d0

  istart_cfg = psi_csf_to_config_data(istart)
  iend_cfg   = psi_csf_to_config_data(iend)

  !nconnectedtotalmax = 1000
  !nconnectedmaxJ = 1000
  maxnalphas = elec_num*mo_num
  Icfg(:,1) = psi_configuration(:,1,1)
  Icfg(:,2) = psi_configuration(:,2,1)
  allocate(listconnectedJ(N_INT,2,max(sze,10000)))
  allocate(idslistconnectedJ(max(sze,10000)))
  call obtain_connected_J_givenI(1, Icfg, listconnectedJ, idslistconnectedJ, nconnectedmaxJ, nconnectedtotalmax)
  deallocate(listconnectedJ)
  deallocate(idslistconnectedJ)

  integer*8, allocatable :: bit_tmp(:)
  integer*8, external :: configuration_search_key
  double precision :: diagfactors_0
  allocate( bit_tmp(0:N_configuration+1))
  do j=1,N_configuration
    bit_tmp(j) = configuration_search_key(psi_configuration(1,1,j),N_int)
  enddo

  call omp_set_max_active_levels(1)
  !$OMP PARALLEL                                              & 
      !$OMP DEFAULT(NONE)                               &
      !$OMP private(i,icfg, isomo, idomo, NSOMOI, NSOMOJ, nholes, k, listholes,&
      !$OMP    holetype, vmotype, nvmos, listvmos, starti, endi,     &
      !$OMP    nsinglesI, singlesI,idxs_singlesI,excitationIds_single,&
      !$OMP    excitationTypes_single, idxI, p, q, extype, pmodel, qmodel,&
      !$OMP    Jsomo, Jdomo, startj, endj, kk, jj, ii, cnti, cntj, meCC1,&
      !$OMP   nconnectedJ,listconnectedJ,idslistconnectedJ,ntotJ,       &
      !$OMP   Nalphas_Icfg,alphas_Icfg,connectedI_alpha,             &
      !$OMP   idxs_connectedI_alpha,nconnectedI,excitationIds,excitationTypes,diagfactors,&
      !$OMP   totcolsTKI,rowsTKI,NSOMOalpha,rowsikpq,                &
      !$OMP   colsikpq, GIJpqrs,TKIGIJ,j,l,m,TKI,CCmattmp, moi, moj, mok, mol,&
      !$OMP   diagfac, tmpvar, diagfactors_0)                                            &
      !$OMP shared(istart_cfg, iend_cfg, psi_configuration, mo_num, psi_config_data,&
      !$OMP    N_int, N_st, psi_out, psi_in, h_core_ri, core_energy, h_act_ri, AIJpqContainer,&
      !$OMP     pp, sze, NalphaIcfg_list,alphasIcfg_list, bit_tmp,       &
      !$OMP     qq, iint, jint, ipos, jpos, nelall, listall, Nsomo_I, countelec,&
      !$OMP     AIJpqMatrixDimsList, diag_energies, n_CSF, lock, NBFmax,nconnectedtotalmax, nconnectedmaxJ,maxnalphas,&
      !$OMP     n_core_orb, n_act_orb, list_act, n, list_core,  list_core_is_built,core_act_contrib, num_threads_max,&
      !$OMP     n_core_orb_is_built, mo_integrals_map, mo_integrals_map_is_built)

  allocate(singlesI(N_INT,2,max(sze,10000)))
  allocate(idxs_singlesI(max(sze,10000)))
  allocate(excitationIds_single(2,max(sze,10000)))
  allocate(excitationTypes_single(max(sze,10000)))
!

  !!!====================!!!
  !!! Single Excitations !!!
  !!!====================!!!

  !$OMP DO SCHEDULE(dynamic,16)
  do i=istart_cfg,iend_cfg

    ! if Seniority_range > 8 then
    ! continue
    ! else
    ! cycle

    do ii=1,N_INT
     Icfg(ii,1) = psi_configuration(ii,1,i)
     Icfg(ii,2) = psi_configuration(ii,2,i)
     Isomo(ii) = Icfg(ii,1)
     Idomo(ii) = Icfg(ii,2)
    enddo
    NSOMOI = getNSOMO(Icfg)

     ! find out all pq holes possible
     nholes = 0
     ! holes in SOMO
     ! list_act
     ! list_core
     ! list_core_inact
     ! bitmasks
     !do k = 1,mo_num
    ! do kk = 1,n_act_orb
    !   k = list_act(kk)
    !   if(POPCNT(IAND(Isomo,IBSET(0_8,k-1))) .EQ. 1) then
    !     nholes += 1
    !     listholes(nholes) = k
    !     holetype(nholes) = 1
    !   endif
    ! enddo
    ! ! holes in DOMO
    ! !do k = 1,mo_num
    ! do kk = 1,n_act_orb
    !   k = list_act(kk)
    !   if(POPCNT(IAND(Idomo,IBSET(0_8,k-1))) .EQ. 1) then
    !     nholes += 1
    !     listholes(nholes) = k
    !     holetype(nholes) = 2
    !   endif
    ! enddo

    ! ! find vmos
    ! do kk = 1,n_act_orb
    !   k = list_act(kk)
    !   !print *,i,IBSET(0,i-1),POPCNT(IAND(Isomo,(IBSET(0_8,i-1)))), POPCNT(IAND(Idomo,(IBSET(0_8,i-1))))
    !   if(POPCNT(IAND(Isomo,(IBSET(0_8,k-1)))) .EQ. 0 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,k-1)))) .EQ. 0) then
    !     nvmos += 1
    !     listvmos(nvmos) = k
    !     vmotype(nvmos) = 0
    !   else if(POPCNT(IAND(Isomo,(IBSET(0_8,k-1)))) .EQ. 1 .AND. POPCNT(IAND(Idomo,(IBSET(0_8,k-1)))) .EQ. 0 ) then
    !     nvmos += 1
    !     listvmos(nvmos) = k
    !     vmotype(nvmos) = 1
    !   end if
    ! enddo

  ! find out all pq holes possible
    nholes = 0
        call bitstring_to_list(Isomo,listall,nelall,N_int)

        do iii=1,nelall
          nholes += 1
          listholes(nholes) = listall(iii)
          holetype(nholes) = 1
        end do

        Nsomo_I = nelall

        call bitstring_to_list(Idomo,listall,nelall,N_int)

        do iii=1,nelall
          if(listall(iii) .gt. n_core_orb)then
            nholes += 1
            listholes(nholes) = listall(iii)
            holetype(nholes) = 2
          endif
        end do


     listvmos = -1
     vmotype = -1
     nvmos = 0
  ! find vmos
    ! Take into account N_int
    do ii = 1, n_act_orb
      iii = list_act(ii)
      iint = shiftr(iii-1,bit_kind_shift) + 1
      ipos = iii-shiftl((iint-1),bit_kind_shift)-1

      if(IAND(Idomo(iint),(IBSET(0_8,ipos))) .EQ. 0) then
        if(IAND(Isomo(iint),(IBSET(0_8,ipos))) .EQ. 0) then
          nvmos += 1
          listvmos(nvmos) = iii
          vmotype(nvmos) = 1
        else if(POPCNT(IAND(Isomo(iint),(IBSET(0_8,ipos)))) .EQ. 1) then
          nvmos += 1
          listvmos(nvmos) = iii
          vmotype(nvmos) = 2
        end if
      end if
    end do



     ! Icsf ids
     starti = psi_config_data(i,1)
     endi   = psi_config_data(i,2)
     NSOMOI = getNSOMO(Icfg)

     call generate_all_singles_cfg_with_type(bit_tmp,Icfg,singlesI,idxs_singlesI,excitationIds_single,&
         excitationTypes_single,nsinglesI,N_int)

     do j = 1,nsinglesI
       idxI = idxs_singlesI(j)
       NSOMOJ = getNSOMO(singlesI(1,1,j))
      p = excitationIds_single(1,j)
      q = excitationIds_single(2,j)
      extype = excitationTypes_single(j)
      ! Off diagonal terms
      call convertOrbIdsToModelSpaceIds(Icfg, singlesI(1,1,j), p, q, extype, pmodel, qmodel)
      do ii=1,N_INT
        Jsomo(ii) = singlesI(1,1,j)
        Jdomo(ii) = singlesI(1,2,j)
      enddo

      ! Get actual p pos
      pp  = p
      iint = shiftr(pp-1,bit_kind_shift) + 1
      ipos = pp-shiftl((iint-1),bit_kind_shift)-1

      ! Get actual q pos
      qq  = q
      jint = shiftr(qq-1,bit_kind_shift) + 1
      jpos = qq-shiftl((jint-1),bit_kind_shift)-1

      ! Add the hole on J
      !if(POPCNT(IAND(Jsomo,IBSET(0_8,q-1))) .EQ. 1  .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
      if(POPCNT(IAND(Jsomo(jint),IBSET(0_8,jpos))) .EQ. 1  .AND. POPCNT(IAND(Isomo(jint),IBSET(0_8,jpos))) .EQ. 0) then
        nholes += 1
        listholes(nholes) = q
        holetype(nholes) = 1
      endif
      !if((POPCNT(IAND(Jdomo,IBSET(0_8,q-1))) .EQ. 1 .AND. POPCNT(IAND(Idomo,IBSET(0_8,q-1))) .EQ. 0) .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
      if((POPCNT(IAND(Jdomo(jint),IBSET(0_8,jpos))) .EQ. 1 .AND. POPCNT(IAND(Idomo(jint),IBSET(0_8,jpos))) .EQ. 0) .AND.&
      POPCNT(IAND(Isomo(jint),IBSET(0_8,jpos))) .EQ. 0) then
        nholes += 1
        listholes(nholes) = q
        holetype(nholes) = 2
      endif

      startj = psi_config_data(idxI,1)
      endj   = psi_config_data(idxI,2)
      !print *,"i=",i," idxI=",idxI," startj=",startj," endj=",endj," sze=",sze

      !!! One-electron contribution !!!
      do ii = starti, endi
        cnti = ii-starti+1
        do jj = startj, endj
          cntj = jj-startj+1
          !meCC1 = AIJpqContainer(cnti,cntj,pmodel,qmodel,extype,NSOMOI)* h_core_ri(p,q)
          core_act_contrib = 0.0d0
          if(p.ne.q)then
            do pp=1,n_core_orb
              n=list_core(pp)
              core_act_contrib += 2.d0 * get_two_e_integral(p,n,q,n,mo_integrals_map) - get_two_e_integral(p,n,n,q,mo_integrals_map)
            end do
          endif
          meCC1 = AIJpqContainer(cnti,cntj,pmodel,qmodel,extype,NSOMOI)* (h_act_ri(p,q) + core_act_contrib)
          !if(jj.eq.1.and.ii.eq.1)then
          !  print *,"CC=",AIJpqContainer(cnti,cntj,pmodel,qmodel,extype,NSOMOI), " p=",p," q=",q
          !endif
          call omp_set_lock(lock(jj))
          do kk = 1,n_st
            psi_out(kk,jj) = psi_out(kk,jj) + meCC1 * psi_in(kk,ii)
          enddo
          call omp_unset_lock(lock(jj))
        enddo
      enddo

      ! Undo setting in listholes
      !if(POPCNT(IAND(Jsomo,IBSET(0_8,q-1))) .EQ. 1  .AND. POPCNT(IAND(Isomo,IBSET(0_8,q-1))) .EQ. 0) then
      if(POPCNT(IAND(Jsomo(jint),IBSET(0_8,jpos))) .EQ. 1  .AND. POPCNT(IAND(Isomo(jint),IBSET(0_8,jpos))) .EQ. 0) then
        nholes -= 1
      endif
      if((POPCNT(IAND(Jdomo(jint),IBSET(0_8,jpos))) .EQ. 1 .AND. POPCNT(IAND(Idomo(jint),IBSET(0_8,jpos))) .EQ. 0) .AND.&
      POPCNT(IAND(Isomo(jint),IBSET(0_8,jpos))) .EQ. 0) then
        nholes -= 1
      endif
    enddo
  enddo
  !$OMP END DO
  deallocate(singlesI)
  deallocate(idxs_singlesI)
  deallocate(excitationIds_single)
  deallocate(excitationTypes_single)

  !print *," singles part psi(1,1)=",psi_out(1,1)
  !do i=1,n_CSF
  !  print *,"i=",i," psi(i)=",psi_out(1,i)
  !enddo
  
  allocate(listconnectedJ(N_INT,2,max(sze,10000)))
  allocate(alphas_Icfg(N_INT,2,max(sze,10000)))
  allocate(connectedI_alpha(N_INT,2,max(sze,10000)))
  allocate(idxs_connectedI_alpha(max(sze,10000)))
  allocate(excitationIds(2,max(sze,10000)))
  allocate(excitationTypes(max(sze,10000)))
  allocate(diagfactors(max(sze,10000)))
  allocate(idslistconnectedJ(max(sze,10000)))
  allocate(CCmattmp(n_st,NBFmax))

  !!!====================!!!
  !!! Double Excitations !!!
  !!!====================!!!
  ! Loop over all selected configurations
  !$OMP DO SCHEDULE(static)
  do i = istart_cfg,iend_cfg

     ! if Seniority_range > 8 then
     ! continue
     ! else
     ! cycle

     do ii=1,N_INT
       Icfg(ii,1) = psi_configuration(ii,1,i)
       Icfg(ii,2) = psi_configuration(ii,2,i)
     enddo
     starti = psi_config_data(i,1)
     endi   = psi_config_data(i,2)

     ! Returns all unique (checking the past) singly excited cfgs connected to I
     Nalphas_Icfg = 0
     ! TODO:
     ! test if size(alphas_Icfg,1) < Nmo**2) then deallocate + allocate

     Nalphas_Icfg = NalphaIcfg_list(i)
     alphas_Icfg(1:n_int,1:2,1:Nalphas_Icfg) = alphasIcfg_list(1:n_int,1:2,i,1:Nalphas_Icfg)
     !if(Nalphas_Icfg .GT. maxnalphas) then
     !  print *,"Nalpha > maxnalpha"
     !endif

     !call obtain_connected_J_givenI(i, Icfg, listconnectedJ, idslistconnectedJ, nconnectedJ, ntotJ)

     ! TODO : remove doubly excited for return
     !print *,"I=",i,"isomo=",psi_configuration(1,1,i),psi_configuration(2,1,i),POPCNT(psi_configuration(1,1,i)),POPCNT(psi_configuration(2,1,i)),&
     !"idomo=",psi_configuration(1,2,i),psi_configuration(2,2,i),POPCNT(psi_configuration(1,2,i)),POPCNT(psi_configuration(2,2,i)), "Nalphas_Icfg=",Nalphas_Icfg
     do k = 1,Nalphas_Icfg
        ! Now generate all singly excited with respect to a given alpha CFG

        !call obtain_connected_I_foralpha_fromfilterdlist(i,nconnectedJ, idslistconnectedJ, &
        !  listconnectedJ, alphas_Icfg(1,1,k),connectedI_alpha,idxs_connectedI_alpha,nconnectedI, &
        !  excitationIds,excitationTypes,diagfactors)

        call obtain_connected_I_foralpha(i, alphas_Icfg(1,1,k), connectedI_alpha, idxs_connectedI_alpha, &
                                         nconnectedI, excitationIds, excitationTypes, diagfactors)

        !if(i .EQ. 218) then
        !   print *,'k=',k,' kcfgSOMO=',alphas_Icfg(1,1,k),alphas_Icfg(2,1,k),' ',POPCNT(alphas_Icfg(1,1,k)),' &
        !   kcfgDOMO=',alphas_Icfg(1,2,k),alphas_Icfg(2,2,k),' ',POPCNT(alphas_Icfg(1,2,k)), " NconnectedI=",nconnectedI
        !   !print *,'k=',k,' kcfgSOMO=',alphas_Icfg(1,1,k),' ',POPCNT(alphas_Icfg(1,1,k)),' &
        !   !kcfgDOMO=',alphas_Icfg(1,2,k),' ',POPCNT(alphas_Icfg(1,2,k)), " NconnectedI=",nconnectedI
        !endif

        
        if(nconnectedI .EQ. 0) then
           cycle
        endif

        ! Here we do 2x the loop. One to count for the size of the matrix, then we compute.
        totcolsTKI = 0
        rowsTKI = -1
        NSOMOalpha = getNSOMO(alphas_Icfg(:,:,k))
        do j = 1,nconnectedI
           NSOMOI = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           !print *,"K=",k,"j=",j, "countelec=",countelec," p=",p," q=",q, " extype=",extype, "NSOMOalpha=",NSOMOalpha," NSOMOI=",NSOMOI, "alphas_Icfg(1,1,k)=",alphas_Icfg(1,1,k), &
           !alphas_Icfg(2,1,k), " domo=",alphas_Icfg(1,2,k), alphas_Icfg(2,2,k), " connected somo=",connectedI_alpha(1,1,j), &
           !connectedI_alpha(2,1,j), " domo=",connectedI_alpha(1,2,j), connectedI_alpha(2,2,j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(1,1,k), connectedI_alpha(1,1,j), p, q, extype, pmodel, qmodel)
           ! for E_pp E_rs and E_ppE_rr case
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,2)
           !if(i.eq.218)then
           !  print *,"j=",j," k=",k,"p=",p,"q=",q,"NSOMOalpha=",NSOMOalpha, "pmodel=",pmodel,"qmodel=",qmodel, "extype=",extype,&
           !  "conn somo=",connectedI_alpha(1,1,j),connectedI_alpha(2,1,j),&
           !  "conn domo=",connectedI_alpha(1,2,j),connectedI_alpha(2,2,j)
           !  do m=1,colsikpq
           !    print *,idxs_connectedI_alpha(j)+m-1
           !  enddo
           !endif
           !print *,"j=",j," Nsomo=",NSOMOalpha," rowsikpq=",rowsikpq," colsikpq=",colsikpq, " p=",pmodel," q=",qmodel, " extyp=",extype
           totcolsTKI += colsikpq
           rowsTKI = rowsikpq
        enddo

        !if(i.eq.1)then
        !  print *,"n_st=",n_st,"rowsTKI=",rowsTKI, " nconnectedI=",nconnectedI, &
        !  "totcolsTKI=",totcolsTKI
        !endif
        allocate(TKI(n_st,rowsTKI,totcolsTKI)) ! coefficients of CSF
        ! Initialize the integral container
        ! dims : (totcolsTKI, nconnectedI)
        allocate(GIJpqrs(totcolsTKI,nconnectedI))  ! gpqrs
        allocate(TKIGIJ(n_st,rowsTKI,nconnectedI))  ! TKI * gpqrs
        !print *,"\t---rowsTKI=",rowsTKI," totCols=",totcolsTKI
        TKI = 0.d0
        GIJpqrs = 0.d0
        TKIGIJ = 0.d0

        totcolsTKI = 0
        do j = 1,nconnectedI
           NSOMOI = getNSOMO(connectedI_alpha(:,:,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(:,:,k), connectedI_alpha(:,:,j), p, q, extype, pmodel, qmodel)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,2)
           rowsTKI = rowsikpq
           !if(i.eq.1) then
           !  print *,rowsTKI,colsikpq," | ",pmodel,qmodel,extype,NSOMOalpha
           !endif
           do m = 1,colsikpq
              do l = 1,rowsTKI
              do kk = 1,n_st
                 TKI(kk,l,totcolsTKI+m) = AIJpqContainer(l,m,pmodel,qmodel,extype,NSOMOalpha) &
                    * psi_in(kk,idxs_connectedI_alpha(j)+m-1)
              enddo
           enddo
           !if(i.eq.1) then
           !      print *,"j=",j,"psi_in=",psi_in(1,idxs_connectedI_alpha(j)+m-1)
           !endif
           enddo

           diagfactors_0 = diagfactors(j)*0.5d0
           moi = excitationIds(1,j) ! p
           mok = excitationIds(2,j) ! q
           do l=1,nconnectedI
            moj = excitationIds(2,l) ! s
            mol = excitationIds(1,l) ! r
            diagfac =  diagfactors_0 * diagfactors(l)* mo_two_e_integral(mok,mol,moi,moj)! g(pq,sr) = <ps,qr>
            !print *,"p=",mok,"q=",mol,"r=",moi,"s=",moj
            do m = 1,colsikpq
               ! <ij|kl> = (ik|jl)
               GIJpqrs(totcolsTKI+m,l) = diagfac
            enddo
           enddo
           totcolsTKI += colsikpq
        enddo


        ! Do big BLAS
        call dgemm('N','N', rowsTKI*n_st, nconnectedI, totcolsTKI, 1.d0,  &
             TKI, size(TKI,1)*size(TKI,2), GIJpqrs, size(GIJpqrs,1), 0.d0, &
             TKIGIJ , size(TKIGIJ,1)*size(TKIGIJ,2) )


        ! Collect the result
        totcolsTKI = 0
        do j = 1,nconnectedI
           NSOMOI     = getNSOMO(connectedI_alpha(1,1,j))
           p = excitationIds(1,j)
           q = excitationIds(2,j)
           extype = excitationTypes(j)
           call convertOrbIdsToModelSpaceIds(alphas_Icfg(1,1,k), connectedI_alpha(1,1,j), p, q, extype, pmodel, qmodel)
           rowsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,1)
           colsikpq = AIJpqMatrixDimsList(NSOMOalpha,extype,pmodel,qmodel,2)
           rowsTKI = rowsikpq
           CCmattmp = 0.d0

        !if(i.eq.1)then
        !  print *,"\t n_st=",n_st," colsikpq=",colsikpq," rowsTKI=",rowsTKI,&
        !    " | ",size(TKIGIJ,1),size(AIJpqContainer,1),size(CCmattmp,1)
        !endif
           call dgemm('N','N', n_st, colsikpq, rowsTKI, 1.d0,        &
               TKIGIJ(1,1,j), size(TKIGIJ,1),                        &
               AIJpqContainer(1,1,pmodel,qmodel,extype,NSOMOalpha),  &
               size(AIJpqContainer,1), 0.d0,                         &
               CCmattmp, size(CCmattmp,1) )

           !print *,"j=",j,"colsikpq=",colsikpq, "sizeTIG=",size(TKIGIJ,1),"sizeaijpq=",size(AIJpqContainer,1)
           do m = 1,colsikpq
              call omp_set_lock(lock(idxs_connectedI_alpha(j)+m-1))
              do kk = 1,n_st
                 psi_out(kk,idxs_connectedI_alpha(j)+m-1) += CCmattmp(kk,m)
                 !if(dabs(CCmattmp(kk,m)).gt.1e-10)then
                 !  print *, CCmattmp(kk,m), " | ",idxs_connectedI_alpha(j)+m-1
                 !end if
              enddo
              call omp_unset_lock(lock(idxs_connectedI_alpha(j)+m-1))
           enddo
           totcolsTKI += colsikpq
        enddo

        deallocate(TKI) ! coefficients of CSF
        deallocate(GIJpqrs)  ! gpqrs
        deallocate(TKIGIJ)  ! gpqrs

     enddo ! loop over alphas
  enddo ! loop over I
  !$OMP END DO
  call omp_set_max_active_levels(4)
  deallocate(CCmattmp)
  deallocate(connectedI_alpha)
  deallocate(idxs_connectedI_alpha)
  deallocate(excitationIds)
  deallocate(excitationTypes)
  deallocate(diagfactors)

  !print *," psi(1,823)=",psi_out(1,823), " g(1 8, 3 15)=",mo_two_e_integral(1,8,3,15), " ncore=",n_core_orb
  !print *," psi(1,1)=",psi_out(1,1)

  ! Add the diagonal contribution
  !$OMP DO
  do i = 1,n_CSF
    do kk=1,n_st
     psi_out(kk,i) += diag_energies(i)*psi_in(kk,i)
    enddo
  enddo
  !$OMP END DO

  !$OMP END PARALLEL
  !print *," ----- "
  !do i=1,sze
  !  print *,"i=",i," psi_out(i)=",psi_out(1,i)
  !end do
  call omp_set_max_active_levels(4)

  deallocate(diag_energies)
  deallocate(bit_tmp)

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
