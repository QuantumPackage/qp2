BEGIN_PROVIDER [double precision, ha_to_ev]
&BEGIN_PROVIDER [double precision, au_to_D]
&BEGIN_PROVIDER [double precision, planck_cte]
&BEGIN_PROVIDER [double precision, light_speed]
&BEGIN_PROVIDER [double precision, Ha_to_J]
&BEGIN_PROVIDER [double precision, Ha_to_nm]

  implicit none

  BEGIN_DOC
  ! Some conversion between different units
  END_DOC

  ! Hartree to eV
  Ha_to_eV = 27.211396641308d0

  ! au to Debye
  au_to_D = 2.5415802529d0

  ! Planck's constant in SI units
  planck_cte = 6.62606957d-34

  ! Light speed in SI units
  light_speed = 2.99792458d10

  ! Hartree to Joule
  Ha_to_J = 4.35974434d-18

  ! Hartree to nm
  Ha_to_nm = 1d9 * (planck_cte * light_speed) / Ha_to_J

END_PROVIDER
