integer, parameter :: max_dim = 511
integer, parameter :: SIMD_vector = 32
integer, parameter :: N_int_max = 128

double precision, parameter :: pi =  dacos(-1.d0)
double precision, parameter :: inv_pi =  1.d0/dacos(-1.d0)
double precision, parameter :: sqpi =  dsqrt(dacos(-1.d0))
double precision, parameter :: pi_5_2 =  34.9868366552d0
double precision, parameter :: dfour_pi =  4.d0*dacos(-1.d0)
double precision, parameter :: dtwo_pi =  2.d0*dacos(-1.d0)
double precision, parameter :: inv_sq_pi =  1.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: c_mu_gauss =  27.d0/(8.d0*dsqrt(dacos(-1.d0)))
double precision, parameter :: c_mu_gauss_tot =  1.5d0*27.d0/(8.d0*dsqrt(dacos(-1.d0)))+3.d0/dsqrt(dacos(-1.d0))
double precision, parameter :: alpha_mu_gauss =  1.5d0
double precision, parameter :: inv_sq_pi_2 = 0.5d0/dsqrt(dacos(-1.d0))
double precision, parameter :: thresh = 1.d-15
double precision, parameter :: cx_lda = -0.73855876638202234d0
double precision, parameter :: c_2_4_3 = 2.5198420997897464d0
double precision, parameter :: cst_lda = -0.93052573634909996d0
double precision, parameter :: c_4_3 = 4.d0/3.d0
double precision, parameter :: c_1_3 = 1.d0/3.d0
double precision, parameter :: sq_op5 = dsqrt(0.5d0)
double precision, parameter :: dlog_2pi =  dlog(2.d0*dacos(-1.d0))
double precision, parameter :: const_rs=  (3.d0/(4.D0*pi))**(1.d0/3.d0)

! physical constants and units conversion factors
double precision, parameter :: k_boltzman_si    = 1.38066d-23       ! K k^-1 
double precision, parameter :: k_boltzman_au    = 3.1667d-6         ! Hartree k^-1 
double precision, parameter :: k_boltzman_m1_au = 315795.26d0       ! Hartree^-1 k 
double precision, parameter :: bohr_radius_si   = 0.529177d-10      ! m
double precision, parameter :: bohr_radius_cm   = 0.529177d-8       ! cm
double precision, parameter :: bohr_radius_angs = 0.529177d0        ! Angstrom
double precision, parameter :: electronmass_si  = 9.10953d-31       ! Kg
double precision, parameter :: electronmass_uma = 5.4858d-4         ! uma
double precision, parameter :: electronvolt_si  = 1.6021892d-19     ! J  
double precision, parameter :: uma_si           = 1.66057d-27       ! Kg
double precision, parameter :: debye_si         = 3.33564d-30       ! coulomb meter
double precision, parameter :: debye_au         = 0.393427228d0     ! e * Bohr
double precision, parameter :: angstrom_to_au   = 1.889727d0        ! au
double precision, parameter :: au_to_ohmcmm1    = 46000.0d0         ! (ohm cm)^-1
double precision, parameter :: au_to_kb         = 294210.0d0        ! kbar
double precision, parameter :: au_to_eV         = 27.211652d0
double precision, parameter :: uma_to_au        = 1822.89d0
double precision, parameter :: au_to_terahertz  = 2.4189d-5
double precision, parameter :: au_to_sec        = 2.4189d-17
double precision, parameter :: au_to_fsec       = 2.4189d-2
double precision, parameter :: Wcm2             = 3.5d16
double precision, parameter :: amconv           = 1.66042d-24/9.1095d-28*0.5d0  !  mass conversion: a.m.u to a.u. (ry)
double precision, parameter :: uakbar           = 147105.d0         !  pressure conversion from ry/(a.u)^3 to k


