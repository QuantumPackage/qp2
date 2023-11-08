 BEGIN_PROVIDER [integer, n_max_fit_slat]
 implicit none
 BEGIN_DOC
! number of gaussian to fit exp(-x)
!
! I took 20 gaussians from the program bassto.f
 END_DOC
 n_max_fit_slat = 20
 END_PROVIDER

 BEGIN_PROVIDER [double precision, coef_fit_slat_gauss, (n_max_fit_slat)]
&BEGIN_PROVIDER [double precision, expo_fit_slat_gauss, (n_max_fit_slat)]
 implicit none
  include 'constants.include.F'
 BEGIN_DOC
 ! fit the exp(-x) as 
 !
 ! \sum_{i = 1, n_max_fit_slat} coef_fit_slat_gauss(i) * exp(-expo_fit_slat_gauss(i) * x**2)
 !
 ! The coefficient are taken from the program bassto.f
 END_DOC


      expo_fit_slat_gauss(01)=30573.77073000000
      coef_fit_slat_gauss(01)=0.00338925525
      expo_fit_slat_gauss(02)=5608.45238100000
      coef_fit_slat_gauss(02)=0.00536433869
      expo_fit_slat_gauss(03)=1570.95673400000
      coef_fit_slat_gauss(03)=0.00818702846
      expo_fit_slat_gauss(04)=541.39785110000
      coef_fit_slat_gauss(04)=0.01202047655
      expo_fit_slat_gauss(05)=212.43469630000
      coef_fit_slat_gauss(05)=0.01711289568
      expo_fit_slat_gauss(06)=91.31444574000
      coef_fit_slat_gauss(06)=0.02376001022
      expo_fit_slat_gauss(07)=42.04087246000
      coef_fit_slat_gauss(07)=0.03229121736
      expo_fit_slat_gauss(08)=20.43200443000
      coef_fit_slat_gauss(08)=0.04303646818
      expo_fit_slat_gauss(09)=10.37775161000
      coef_fit_slat_gauss(09)=0.05624657578
      expo_fit_slat_gauss(10)=5.46880754500
      coef_fit_slat_gauss(10)=0.07192311571
      expo_fit_slat_gauss(11)=2.97373529200
      coef_fit_slat_gauss(11)=0.08949389001
      expo_fit_slat_gauss(12)=1.66144190200
      coef_fit_slat_gauss(12)=0.10727599240
      expo_fit_slat_gauss(13)=0.95052560820
      coef_fit_slat_gauss(13)=0.12178961750
      expo_fit_slat_gauss(14)=0.55528683970
      coef_fit_slat_gauss(14)=0.12740141870
      expo_fit_slat_gauss(15)=0.33043360020
      coef_fit_slat_gauss(15)=0.11759168160
      expo_fit_slat_gauss(16)=0.19982303230
      coef_fit_slat_gauss(16)=0.08953504394
      expo_fit_slat_gauss(17)=0.12246840760
      coef_fit_slat_gauss(17)=0.05066721317
      expo_fit_slat_gauss(18)=0.07575825322
      coef_fit_slat_gauss(18)=0.01806363869
      expo_fit_slat_gauss(19)=0.04690146243
      coef_fit_slat_gauss(19)=0.00305632563
      expo_fit_slat_gauss(20)=0.02834749861
      coef_fit_slat_gauss(20)=0.00013317513



END_PROVIDER 

double precision function slater_fit_gam(x,gam)
 implicit none
 double precision, intent(in) :: x,gam
 BEGIN_DOC
! fit of the function exp(-gam * x) with gaussian functions 
 END_DOC
 integer :: i
 slater_fit_gam = 0.d0
 do i = 1, n_max_fit_slat
  slater_fit_gam += coef_fit_slat_gauss(i) * dexp(-expo_fit_slat_gauss(i) * gam * gam * x * x)
 enddo
end

subroutine expo_fit_slater_gam(gam,expos)
 implicit none
 BEGIN_DOC
! returns the array of the exponents of the gaussians to fit exp(-gam*x)
 END_DOC
 double precision, intent(in)  :: gam
 double precision, intent(out) :: expos(n_max_fit_slat)
 integer :: i
 do i = 1, n_max_fit_slat
  expos(i) = expo_fit_slat_gauss(i) * gam * gam
 enddo
end

