 real*8 function lgamma(x)
 implicit none
 real*8, intent(in) :: x
 lgamma = log(abs(gamma(x)))
 end function lgamma
  
  BEGIN_PROVIDER [ integer, NSOMOMax]
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
  real*8 :: coeff, binom1, binom2
  integer MS
  integer ncfgpersomo
  real*8, external :: lgamma
  detDimperBF = 0
  MS = elec_alpha_num-elec_beta_num
  ! number of cfgs = number of dets for 0 somos
  n_CSF = 0
  ncfgprev = cfg_seniority_index(0)
  ncfgpersomo = ncfgprev
  do i = 1, elec_num
    print *,"i=",i," Ncfg= ",cfg_seniority_index(i)
  enddo
  do i = iand(MS,1), NSOMOMax-2,2
    if(cfg_seniority_index(i) .EQ. -1) then
      cycle
    endif
    if(cfg_seniority_index(i+2) .EQ. -1) then
      ncfgpersomo = N_configuration + 1
    else
      if(cfg_seniority_index(i+2) > ncfgpersomo) then
          ncfgpersomo = cfg_seniority_index(i+2)
      else
        k = 0
        do while(cfg_seniority_index(i+2+k) < ncfgpersomo)
          k = k + 2
          ncfgpersomo = cfg_seniority_index(i+2+k)
        enddo
      endif
    endif
    ncfg = ncfgpersomo - ncfgprev
    if(iand(MS,1) .EQ. 0) then
      !dimcsfpercfg = max(1,nint((binom(i,i/2)-binom(i,i/2+1))))
      binom1 = dexp(lgamma(1.0d0*(i+1))                            &
                  - lgamma(1.0d0*((i/2)+1))                        &
                  - lgamma(1.0d0*(i-((i/2))+1)));
      binom2 = dexp(lgamma(1.0d0*(i+1))                            &
                  - lgamma(1.0d0*(((i/2)+1)+1))                    &
                  - lgamma(1.0d0*(i-((i/2)+1)+1)));
      dimcsfpercfg = max(1,nint(binom1 - binom2))
    else
      !dimcsfpercfg = max(1,nint((binom(i,(i+1)/2)-binom(i,(i+3)/2))))
      binom1 = dexp(lgamma(1.0d0*(i+1))                            &
                  - lgamma(1.0d0*(((i+1)/2)+1))                    &
                  - lgamma(1.0d0*(i-(((i+1)/2))+1)));
      binom2 = dexp(lgamma(1.0d0*(i+1))                            &
                  - lgamma(1.0d0*((((i+3)/2)+1)+1))                &
                  - lgamma(1.0d0*(i-(((i+3)/2)+1)+1)));
      dimcsfpercfg = max(1,nint(binom1 - binom2))
    endif
    n_CSF += ncfg * dimcsfpercfg
    print *,"i=",i," ncfg= ", ncfg, " dims=", dimcsfpercfg, " n_csf=", n_CSF, ncfgpersomo, ncfgprev
    if(cfg_seniority_index(i+2) > ncfgprev) then
      ncfgprev = cfg_seniority_index(i+2)
    else
      k = 0
      do while(cfg_seniority_index(i+2+k) < ncfgprev)
        k = k + 2
        ncfgprev = cfg_seniority_index(i+2+k)
      enddo
    endif
  enddo
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
  integer                        :: k

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



  BEGIN_PROVIDER [ integer, AIJpqMatrixDimsList, (0:NSOMOMax,0:NSOMOMax,4,NSOMOMax,NSOMOMax,2)]
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
  !allocate(AIJpqMatrixDimsList(NSOMOMax,NSOMOMax,4,NSOMOMax,NSOMOMax,2))
  ! Type
  ! 1. SOMO -> SOMO
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
  do i = 2-iand(nsomomin,1), NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i,i, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,i
           do l = 1,i
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
  ! 4. DOMO -> VMO
  !print *,"Doing DOMO->SOMO"
  do i = 2-iand(nsomomin,1), NSOMOMax, 2
     do j = i,i, 2
        if(j .GT. NSOMOMax .OR. j .LE. 0) then
           cycle
        end if
        do k = 1,i
           do l = 1,i
              if(k .NE. l) then
                 Isomo = ISHFT(1_8,i+1)-1
                 Isomo = IBCLR(Isomo,k+1-1)
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
  END_PROVIDER

  BEGIN_PROVIDER [ real*8, AIJpqContainer, (0:NSOMOMax,0:NSOMOMax,4,NSOMOMax,NSOMOMax,NBFMax,NBFMax)]
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
  !allocate(AIJpqMatrixDimsList(NSOMOMax,NSOMOMax,4,NSOMOMax,NSOMOMax,2))
  ! Type
  ! 1. SOMO -> SOMO
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i-2,i-2, 2
        if(j .GT. NSOMOMax .OR. j .LT. 0) cycle
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
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i,i, 2
        Jsomo = ISHFT(1_8,j)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = 1,i
           do l = 1,i
              if(k .NE. l) then
                 Isomo = ISHFT(1_8,i+1)-1
                 Isomo = IBCLR(Isomo,l-1)
                 Jsomo = ISHFT(1_8,j+1)-1
                 Jsomo = IBCLR(Jsomo,k-1)
              else
                 Isomo = ISHFT(1_8,i)-1
                 Jsomo = ISHFT(1_8,j)-1
              endif

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
  do i = 2, NSOMOMax, 2
     Isomo = ISHFT(1_8,i)-1
     do j = i,i, 2
        Jsomo = ISHFT(1_8,i)-1
        if(j .GT. NSOMOMax .OR. j .LE. 0) cycle
        do k = 1,i
           do l = 1,i
              if(k .NE. l) then
                 Isomo = ISHFT(1_8,i+1)-1
                 Isomo = IBCLR(Isomo,k-1)
                 Jsomo = ISHFT(1_8,j+1)-1
                 Jsomo = IBCLR(Jsomo,l+1-1)
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


!!!!!!

  BEGIN_PROVIDER [ real*8, DetToCSFTransformationMatrix, (0:NSOMOMax,NBFMax,maxDetDimPerBF)]
 &BEGIN_PROVIDER [ real*8, psi_coef_config,  (n_CSF)]
 &BEGIN_PROVIDER [ integer, psi_config_data, (N_configuration,2)]
  use cfunctions
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Documentation for DetToCSFTransformationMatrix
  ! Provides the matrix of transformatons for the
  ! conversion between determinant to CSF basis (in BFs)
  END_DOC
  integer(bit_kind) :: mask(N_int), Ialpha(N_int),Ibeta(N_int)
  integer   :: rows, cols, i, j, k
  integer   :: startdet, enddet
  integer*8 MS, Isomo, Idomo
  integer ndetI
  integer :: getNSOMO
  real*8,dimension(:,:),allocatable    :: tempBuffer
  real*8,dimension(:),allocatable    :: tempCoeff
  real*8  :: norm_det1, phasedet
  norm_det1 = 0.d0
  MS = elec_alpha_num - elec_beta_num
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
    DetToCSFTransformationMatrix(i,1:bfIcfg,1:ndetI) =  tempBuffer(1:bfIcfg,1:ndetI)
    deallocate(tempBuffer)
  enddo

  integer s, bfIcfg
  integer countcsf
  countcsf = 0
  integer countdet
  countdet = 0
  integer idx
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
         idx = psi_configuration_to_psi_det_data(j)
         Ialpha(:) = psi_det(:,1,idx)
         Ibeta(:)  = psi_det(:,2,idx)
         call get_phase_qp_to_cfg(Ialpha, Ibeta, phasedet)
         tempCoeff(countdet) = psi_coef(idx, istate)*phasedet
         norm_det1 += tempCoeff(countdet)*tempCoeff(countdet)
         countdet += 1
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

       call dgemm('N','N', bfIcfg, 1, ndetI, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, size(tempCoeff,1), 0.d0, psi_coef_config(countcsf+1), size(psi_coef_config,1))
       !call dgemv('N', NBFMax, maxDetDimPerBF, 1.d0, tempBuffer, size(tempBuffer,1), tempCoeff, 1, 0.d0, psi_coef_config(countcsf), 1)

      deallocate(tempCoeff)
      deallocate(tempBuffer)
      psi_config_data(i,1) = countcsf + 1
      countcsf += bfIcfg
      psi_config_data(i,2) = countcsf
  enddo

  END_PROVIDER
