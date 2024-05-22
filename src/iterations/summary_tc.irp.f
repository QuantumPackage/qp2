subroutine print_summary_tc(e_,pt2_data,pt2_data_err,n_det_,n_configuration_,n_st,s2_)
  use selection_types
  implicit none
  BEGIN_DOC
! Print the extrapolated energy in the output
  END_DOC

  integer, intent(in)            :: n_det_, n_configuration_, n_st
  double precision, intent(in)   :: e_(n_st), s2_(n_st)
  type(pt2_type)  , intent(in)   :: pt2_data, pt2_data_err
  integer                        :: i, k
  integer                        :: N_states_p
  character*(9)                  :: pt2_string
  character*(512)                :: fmt
  double precision, allocatable :: pt2_minus(:),pt2_plus(:),pt2_tot(:), pt2_abs(:),pt1_norm(:),rpt2_tot(:)
  double precision, allocatable :: error_pt2_minus(:), error_pt2_plus(:), error_pt2_tot(:), error_pt2_abs(:)

  if (do_pt2) then
    pt2_string = '        '
  else
    pt2_string = '(approx)'
  endif

  N_states_p = min(N_det_,n_st)

  allocate(pt2_minus(N_states_p),pt2_plus(N_states_p),pt2_tot(N_states_p), pt2_abs(N_states_p),pt1_norm(N_states_p),rpt2_tot(N_states_p))
  allocate(error_pt2_minus(N_states_p), error_pt2_plus(N_states_p), error_pt2_tot(N_states_p), error_pt2_abs(N_states_p))
  do k = 1, N_states_p
    pt2_plus(k)  = pt2_data % variance(k)
    pt2_minus(k) = pt2_data % pt2(k)
    pt2_abs(k)   = pt2_plus(k) - pt2_minus(k)
    pt2_tot(k)   = pt2_plus(k) + pt2_minus(k)
    pt1_norm(k)  = pt2_data % overlap(k,k)
    rpt2_tot(k)  = pt2_tot(k) / (1.d0 + pt1_norm(k))
    error_pt2_minus(k) = pt2_data_err % pt2(k)
    error_pt2_plus(k)  = pt2_data_err % variance(k)
    error_pt2_tot(k)   = dsqrt(error_pt2_minus(k)**2+error_pt2_plus(k)**2)
    error_pt2_abs(k)   = error_pt2_tot(k) ! same variance because independent variables 
  enddo
  k=1
  write(*,'(A40,X,I10,X,100(F16.8,X))')'Ndet,E,E+PT2,pt2_minus,pt2_plus,pt2_abs=',n_det_,e_(k),e_(k) + pt2_tot(k),e_(k) + rpt2_tot(k),pt2_minus(k), pt2_plus(k),pt2_abs(k)

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
  write(*,fmt) '# PT2 '//pt2_string, (pt2_tot(k), error_pt2_tot(k), k=1,N_states_p)
  write(*,fmt) '# rPT2'//pt2_string, (rpt2_tot(k), error_pt2_tot(k), k=1,N_states_p)
  write(*,'(A)') '#'
  write(*,fmt) '# E+PT2      ', (e_(k)+pt2_tot(k) ,error_pt2_tot(k), k=1,N_states_p)
  write(*,fmt) '# E+rPT2     ', (e_(k)+rpt2_tot(k),error_pt2_tot(k), k=1,N_states_p)
  if (N_states_p > 1) then
    write(*,fmt) '# Excit. (au)', ( (e_(k)+pt2_tot(k)-e_(1)-pt2_tot(1)), &
      dsqrt(error_pt2_tot(k)*error_pt2_tot(k)+error_pt2_tot(1)*error_pt2_tot(1)), k=1,N_states_p)
    write(*,fmt) '# Excit. (eV)', ( (e_(k)+pt2_tot(k)-e_(1)-pt2_tot(1))*27.211396641308d0, &
      dsqrt(error_pt2_tot(k)*error_pt2_tot(k)+error_pt2_tot(1)*error_pt2_tot(1))*27.211396641308d0, k=1,N_states_p)
  endif
  write(fmt,*) '(''# ============'',', N_states_p, '(1X,''=============================''))'
  write(*,fmt)
  print *,  ''

  print *,  'N_det             = ', N_det_
  print *,  'N_states          = ', n_st
  if (s2_eig) then
    print *,  'N_cfg             = ', N_configuration_
    if (only_expected_s2) then
      print *,  'N_csf             = ', N_csf
    endif
  endif
  print *,  ''

  do k=1, N_states_p
    print*,'* State ',k
    print *,  '< S^2 >         = ', s2_(k)
    print *,  'E               = ', e_(k)
    print *,  'PT norm         = ', pt1_norm(k)
    print *,  'PT2             = ', pt2_tot(k),  ' +/- ', error_pt2_tot(k)
    print *,  'rPT2            = ', rpt2_tot(k), ' +/- ', error_pt2_tot(k)
    print *,  'E+PT2 '//pt2_string//' = ', e_(k)+pt2_tot(k) , ' +/- ', error_pt2_tot(k)
    print *,  'E+rPT2'//pt2_string//' = ', e_(k)+rpt2_tot(k), ' +/- ', error_pt2_tot(k)
    print *,  'Positive PT2    = ',pt2_plus(k),' +/- ',error_pt2_plus(k)
    print *,  'Negative PT2    = ',pt2_minus(k),' +/- ',error_pt2_minus(k)
    print *,  'Abs PT2         = ',pt2_abs(k),  ' +/- ',error_pt2_abs(k) 
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
      print*,'Delta E = ', (e_(i)+ pt2_tot(i) - (e_(1) + pt2_tot(1))), &
        (e_(i)+ pt2_tot(i) - (e_(1) + pt2_tot(1))) * 27.211396641308d0
    enddo
    print *,  '-----'
    print*, 'Variational + renormalized perturbative Energy difference (au | eV)'
    do i=2, N_states_p
      print*,'Delta E = ', (e_(i)+ rpt2_tot(i) - (e_(1) + rpt2_tot(1))), &
        (e_(i)+ rpt2_tot(i) - (e_(1) + rpt2_tot(1))) * 27.211396641308d0
    enddo
  endif

!  call print_energy_components()

end subroutine

