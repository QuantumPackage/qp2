subroutine print_summary(e_,pt2_data,pt2_data_err,n_det_,n_occ_pattern_,n_st,s2_)
  use selection_types
  implicit none
  BEGIN_DOC
! Print the extrapolated energy in the output
  END_DOC

  integer, intent(in)            :: n_det_, n_occ_pattern_, n_st
  double precision, intent(in)   :: e_(n_st), s2_(n_st)
  type(pt2_type)  , intent(in)   :: pt2_data, pt2_data_err
  integer                        :: i, k
  integer                        :: N_states_p
  character*(9)                  :: pt2_string
  character*(512)                :: fmt
  double precision               :: f(n_st)

  if (do_pt2) then
    pt2_string = '        '
  else
    pt2_string = '(approx)'
  endif

  N_states_p = min(N_det_,n_st)

  do i=1,N_states_p
    f(i) = 1.d0/(1.d0+pt2_data % norm2(i))
  enddo

  print *, ''
  print '(A,I12)',  'Summary at N_det = ', N_det_
  print '(A)',      '-----------------------------------'
  print *, ''

  write(fmt,*) '(''# ============'',', N_states_p, '(1X,''=============================''))'
  write(*,fmt)
  write(fmt,*) '(13X,', N_states_p, '(6X,A7,1X,I6,10X))'
  write(*,fmt) ('State',k, k=1,N_states_p)
  write(fmt,*) '(''# ============'',', N_states_p, '(1X,''=============================''))'
  write(*,fmt)
  write(fmt,*) '(A13,', N_states_p, '(1X,F14.8,15X))'
  write(*,fmt) '# E          ', e_(1:N_states_p)
  if (N_states_p > 1) then
    write(*,fmt) '# Excit. (au)', e_(1:N_states_p)-e_(1)
    write(*,fmt) '# Excit. (eV)', (e_(1:N_states_p)-e_(1))*27.211396641308d0
  endif
  write(fmt,*) '(A13,', 2*N_states_p, '(1X,F14.8))'
  write(*,fmt) '# PT2 '//pt2_string, (pt2_data % pt2(k), pt2_data_err % pt2(k), k=1,N_states_p)
  write(*,fmt) '# rPT2'//pt2_string, (pt2_data % pt2(k)*f(k), pt2_data_err % pt2(k)*f(k), k=1,N_states_p)
  write(*,'(A)') '#'
  write(*,fmt) '# E+PT2      ', (e_(k)+pt2_data % pt2(k),pt2_data_err % pt2(k), k=1,N_states_p)
  write(*,fmt) '# E+rPT2     ', (e_(k)+pt2_data % pt2(k)*f(k),pt2_data_err % pt2(k)*f(k), k=1,N_states_p)
  if (N_states_p > 1) then
    write(*,fmt) '# Excit. (au)', ( (e_(k)+pt2_data % pt2(k)-e_(1)-pt2_data % pt2(1)), &
      dsqrt(pt2_data_err % pt2(k)*pt2_data_err % pt2(k)+pt2_data_err % pt2(1)*pt2_data_err % pt2(1)), k=1,N_states_p)
    write(*,fmt) '# Excit. (eV)', ( (e_(k)+pt2_data % pt2(k)-e_(1)-pt2_data % pt2(1))*27.211396641308d0, &
      dsqrt(pt2_data_err % pt2(k)*pt2_data_err % pt2(k)+pt2_data_err % pt2(1)*pt2_data_err % pt2(1))*27.211396641308d0, k=1,N_states_p)
  endif
  write(fmt,*) '(''# ============'',', N_states_p, '(1X,''=============================''))'
  write(*,fmt)
  print *,  ''

  print *,  'N_det             = ', N_det_
  print *,  'N_states          = ', n_st
  if (s2_eig) then
    print *,  'N_sop             = ', N_occ_pattern_
  endif
  print *,  ''

  do k=1, N_states_p
    print*,'* State ',k
    print *,  '< S^2 >         = ', s2_(k)
    print *,  'E               = ', e_(k)
    print *,  'Variance        = ', pt2_data % variance(k), ' +/- ', pt2_data_err % variance(k)
    print *,  'PT norm         = ', dsqrt(pt2_data % norm2(k)), ' +/- ', 0.5d0*dsqrt(pt2_data % norm2(k)) * pt2_data_err % norm2(k) / pt2_data % norm2(k)
    print *,  'PT2             = ', pt2_data % pt2(k), ' +/- ', pt2_data_err % pt2(k)
    print *,  'rPT2            = ', pt2_data % pt2(k)*f(k), ' +/- ', pt2_data_err % rpt2(k)
    print *,  'E+PT2 '//pt2_string//' = ', e_(k)+pt2_data % pt2(k), ' +/- ', pt2_data_err % pt2(k)
    print *,  'E+rPT2'//pt2_string//' = ', e_(k)+pt2_data % pt2(k)*f(k), ' +/- ', pt2_data_err % pt2(k)*f(k)
    print *,  ''
  enddo

  print *,  '-----'
  if(n_st.gt.1)then
    print *, 'Variational Energy difference (au | eV)'
    do i=2, N_states_p
      print*,'Delta E = ', (e_(i) - e_(1)), &
        (e_(i) - e_(1)) * 27.211396641308d0
    enddo
    print *,  '-----'
    print*, 'Variational + perturbative Energy difference (au | eV)'
    do i=2, N_states_p
      print*,'Delta E = ', (e_(i)+ pt2_data % pt2(i) - (e_(1) + pt2_data % pt2(1))), &
        (e_(i)+ pt2_data % pt2(i) - (e_(1) + pt2_data % pt2(1))) * 27.211396641308d0
    enddo
    print *,  '-----'
    print*, 'Variational + renormalized perturbative Energy difference (au | eV)'
    do i=2, N_states_p
      print*,'Delta E = ', (e_(i)+ pt2_data % pt2(i)*f(i) - (e_(1) + pt2_data % pt2(1)*f(1))), &
        (e_(i)+ pt2_data % pt2(i)*f(i) - (e_(1) + pt2_data % pt2(1)*f(1))) * 27.211396641308d0
    enddo
  endif

  call print_energy_components()

end subroutine

