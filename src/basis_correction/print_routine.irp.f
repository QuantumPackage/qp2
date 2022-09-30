subroutine print_basis_correction
 implicit none
 integer :: istate
 provide mu_average_prov
 if(mu_of_r_potential.EQ."hf")then
  provide ecmd_lda_mu_of_r ecmd_pbe_ueg_mu_of_r
 else if(mu_of_r_potential.EQ."cas_ful".or.mu_of_r_potential.EQ."cas_truncated")then
  provide ecmd_lda_mu_of_r ecmd_pbe_ueg_mu_of_r 
  provide ecmd_pbe_on_top_mu_of_r ecmd_pbe_on_top_su_mu_of_r
 endif

 print*, ''
 print*, ''
 print*, '****************************************'
 print*, '****************************************'
 print*, 'Basis set correction for WFT using DFT Ecmd functionals'
 print*, 'These functionals are accurate for short-range correlation'
 print*, ''
 print*, 'For more details look at Journal of Chemical Physics 149, 194301 1-15 (2018)        '
 print*, '                         Journal of Physical Chemistry Letters 10, 2931-2937 (2019) '
 print*, '                         ???REF SC?'
 print*, '****************************************'
 print*, '****************************************'
 print*, 'mu_of_r_potential = ',mu_of_r_potential
 if(mu_of_r_potential.EQ."hf")then
   print*, ''
   print*,'Using a HF-like two-body density to define mu(r)'
   print*,'This assumes that HF is a qualitative representation of the wave function ' 
   print*,'********************************************'
   print*,'Functionals more suited for weak correlation'
   print*,'********************************************'
   print*,'+) LDA Ecmd functional     : purely based on the UEG (JCP,149,194301,1-15 (2018)) '
   do istate = 1, N_states
    write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD LDA           , state ',istate,' = ',ecmd_lda_mu_of_r(istate)
   enddo
   print*,'+) PBE-UEG Ecmd functional : PBE at mu=0, UEG ontop pair density at large mu (JPCL, 10, 2931-2937 (2019))'
   do istate = 1, N_states
    write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD PBE-UEG       , state ',istate,' = ',ecmd_pbe_ueg_mu_of_r(istate)
   enddo

  else if(mu_of_r_potential.EQ."cas_ful")then
   print*, ''
   print*,'Using a CAS-like two-body density to define mu(r)'
   print*,'This assumes that the CAS is a qualitative representation of the wave function ' 
   print*,'********************************************'
   print*,'Functionals more suited for weak correlation'
   print*,'********************************************'
   print*,'+) LDA Ecmd functional     : purely based on the UEG (JCP,149,194301,1-15 (2018)) '
   do istate = 1, N_states
    write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD LDA           , state ',istate,' = ',ecmd_lda_mu_of_r(istate)
   enddo
   print*,'+) PBE-UEG Ecmd functional : PBE at mu=0, UEG ontop pair density at large mu (JPCL, 10, 2931-2937 (2019))'
   do istate = 1, N_states
    write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD PBE-UEG       , state ',istate,' = ',ecmd_pbe_ueg_mu_of_r(istate)
   enddo
   print*,''
   print*,'********************************************'
   print*,'********************************************'
   print*,'+) PBE-on-top Ecmd functional : (??????? REF-SCF ??????????)' 
   print*,'PBE at mu=0, extrapolated ontop pair density at large mu, usual spin-polarization'
   do istate = 1, N_states
    write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD PBE-OT        , state ',istate,' = ',ecmd_pbe_on_top_mu_of_r(istate)
   enddo
   print*,''
   print*,'********************************************'
   print*,'+) PBE-on-top no spin polarization Ecmd functional : (??????? REF-SCF ??????????)' 
   print*,'PBE at mu=0, extrapolated ontop pair density at large mu, and ZERO SPIN POLARIZATION'
   do istate = 1, N_states
    write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD SU-PBE-OT     , state ',istate,' = ',ecmd_pbe_on_top_su_mu_of_r(istate)
   enddo
   print*,''

  endif
  print*,''
  print*,'**************'
  do istate = 1, N_states
    write(*, '(A29,X,I3,X,A3,X,F16.10)') '  Average mu(r)      , state ',istate,' = ',mu_average_prov(istate)
  enddo

end



subroutine print_all_basis_correction
 implicit none
 integer :: istate
 provide mu_average_prov
 provide ecmd_lda_mu_of_r ecmd_pbe_ueg_mu_of_r 
 provide ecmd_pbe_on_top_mu_of_r ecmd_pbe_on_top_su_mu_of_r

 print*, ''
 print*, ''
 print*, '****************************************'
 print*, '****************************************'
 print*, 'Basis set correction for WFT using DFT Ecmd functionals'
 print*, 'These functionals are accurate for short-range correlation'
 print*, ''
 print*, 'For more details look at Journal of Chemical Physics 149, 194301 1-15 (2018)        '
 print*, '                         Journal of Physical Chemistry Letters 10, 2931-2937 (2019) '
 print*, '                         ???REF SC?'
 print*, '****************************************'
 print*, '****************************************'
 print*, 'mu_of_r_potential = ',mu_of_r_potential
 print*, ''
 print*,'Using a CAS-like two-body density to define mu(r)'
 print*,'This assumes that the CAS is a qualitative representation of the wave function ' 
 print*,'********************************************'
 print*,'Functionals more suited for weak correlation'
 print*,'********************************************'
 print*,'+) LDA Ecmd functional     : purely based on the UEG (JCP,149,194301,1-15 (2018)) '
 do istate = 1, N_states
  write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD LDA           , state ',istate,' = ',ecmd_lda_mu_of_r(istate)
 enddo
 print*,'+) PBE-UEG Ecmd functional : PBE at mu=0, UEG ontop pair density at large mu (JPCL, 10, 2931-2937 (2019))'
 do istate = 1, N_states
  write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD PBE-UEG       , state ',istate,' = ',ecmd_pbe_ueg_mu_of_r(istate)
 enddo
 print*,''
 print*,'********************************************'
 print*,'********************************************'
 print*,'+) PBE-on-top Ecmd functional : (??????? REF-SCF ??????????)' 
 print*,'PBE at mu=0, extrapolated ontop pair density at large mu, usual spin-polarization'
 do istate = 1, N_states
  write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD PBE-OT        , state ',istate,' = ',ecmd_pbe_on_top_mu_of_r(istate)
 enddo
 print*,''
 print*,'********************************************'
 print*,'+) PBE-on-top no spin polarization Ecmd functional : (??????? REF-SCF ??????????)' 
 print*,'PBE at mu=0, extrapolated ontop pair density at large mu, and ZERO SPIN POLARIZATION'
 do istate = 1, N_states
  write(*, '(A29,X,I3,X,A3,X,F16.10)') '  ECMD SU-PBE-OT     , state ',istate,' = ',ecmd_pbe_on_top_su_mu_of_r(istate)
 enddo
 print*,''

  print*,''
  print*,'**************'
  do istate = 1, N_states
    write(*, '(A29,X,I3,X,A3,X,F16.10)') '  Average mu(r)      , state ',istate,' = ',mu_average_prov(istate)
  enddo

end


