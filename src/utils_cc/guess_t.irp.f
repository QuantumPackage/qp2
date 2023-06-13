! T1

subroutine guess_t1(nO,nV,f_o,f_v,f_ov,t1)

  implicit none

  BEGIN_DOC
  ! Update the T1 amplitudes for CC
  END_DOC

  ! in
  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: f_o(nO), f_v(nV), f_ov(nO,nV)

  ! inout
  double precision, intent(out) :: t1(nO, nV)

  ! internal
  integer                       :: i,a

  if (trim(cc_guess_t1) == 'none') then
     t1 = 0d0
  else if (trim(cc_guess_t1) == 'MP') then
    do a = 1, nV
      do i = 1, nO
        t1(i,a) = f_ov(i,a) / (f_o(i) - f_v(a) - cc_level_shift_guess)
      enddo
    enddo
  else if (trim(cc_guess_t1) == 'read') then
    call read_t1(nO,nV,t1)
  else
    print*, 'Unknown cc_guess_t1 type: '//trim(cc_guess_t1)
    call abort
  endif
  
end

! T2

subroutine guess_t2(nO,nV,f_o,f_v,v_oovv,t2)

  implicit none

  BEGIN_DOC
  ! Update the T2 amplitudes for CC
  END_DOC

  ! in
  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: f_o(nO), f_v(nV), v_oovv(nO, nO, nV, nV)

  ! inout
  double precision, intent(out) :: t2(nO, nO, nV, nV)

  ! internal
  integer                       :: i,j,a,b

  if (trim(cc_guess_t2) == 'none') then
    t2 = 0d0
  else if (trim(cc_guess_t2) == 'MP') then
    do b = 1, nV
      do a = 1, nV
        do j = 1, nO
          do i = 1, nO
            t2(i,j,a,b) = v_oovv(i,j,a,b) / (f_o(i) + f_o(j) - f_v(a) - f_v(b) - cc_level_shift_guess)
          enddo
        enddo
      enddo
    enddo
  else if (trim(cc_guess_t2) == 'read') then
    call read_t2(nO,nV,t2)
  else
    print*, 'Unknown cc_guess_t1 type: '//trim(cc_guess_t2)
    call abort
  endif
  
end

! T1

subroutine write_t1(nO,nV,t1)

  implicit none

  BEGIN_DOC
  ! Write the T1 amplitudes for CC
  END_DOC

  ! in
  integer, intent(in)          :: nO, nV
  double precision, intent(in) :: t1(nO, nV)

  ! internal
  integer                      :: i,a, iunit
  integer, external :: getunitandopen

  if (write_amplitudes) then
    iunit = getUnitAndOpen(trim(ezfio_filename)//'/work/T1','w')
    do a = 1, nV
      do i = 1, nO
         write(iunit,'(F20.12)') t1(i,a)
      enddo
    enddo
    close(iunit)
  endif
  
end

! T2

subroutine write_t2(nO,nV,t2)

  implicit none

  BEGIN_DOC
  ! Write the T2 amplitudes for CC
  END_DOC

  ! in
  integer, intent(in)          :: nO, nV
  double precision, intent(in) :: t2(nO, nO, nV, nV)

  ! internal
  integer                      :: i,j,a,b, iunit
  integer, external :: getunitandopen

  if (write_amplitudes) then
    iunit = getUnitAndOpen(trim(ezfio_filename)//'/work/T2','w')
    do b = 1, nV
      do a = 1, nV
        do j = 1, nO
          do i = 1, nO
             write(iunit,'(F20.12)') t2(i,j,a,b)
          enddo
        enddo
      enddo
    enddo
    close(iunit)
  endif
  
end

! T1

subroutine read_t1(nO,nV,t1)

  implicit none

  BEGIN_DOC
  ! Read the T1 amplitudes for CC
  END_DOC

  ! in
  integer, intent(in)           :: nO, nV
  double precision, intent(out) :: t1(nO, nV)

  ! internal
  integer                       :: i,a, iunit
  logical                       :: ok
  integer, external :: getunitandopen

  if (read_amplitudes) then
    iunit = getUnitAndOpen(trim(ezfio_filename)//'/work/T1','r')
    do a = 1, nV
      do i = 1, nO
         read(iunit,'(F20.12)') t1(i,a)
      enddo
    enddo
    close(iunit)
  endif
  
end

! T2

subroutine read_t2(nO,nV,t2)

  implicit none

  BEGIN_DOC
  ! Read the T2 amplitudes for CC
  END_DOC

  ! in
  integer, intent(in)           :: nO, nV
  double precision, intent(out) :: t2(nO, nO, nV, nV)

  ! internal
  integer                       :: i,j,a,b, iunit
  logical                       :: ok

  integer, external :: getunitandopen

  if (read_amplitudes) then
    iunit = getUnitAndOpen(trim(ezfio_filename)//'/work/T2','r')
    do b = 1, nV
      do a = 1, nV
        do j = 1, nO
          do i = 1, nO
             read(iunit,'(F20.12)') t2(i,j,a,b)
          enddo
        enddo
      enddo
    enddo
    close(iunit)
  endif
  
end
