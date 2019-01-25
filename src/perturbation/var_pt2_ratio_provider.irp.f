BEGIN_PROVIDER [ double precision, var_pt2_ratio  ]
  implicit none
  BEGIN_DOC
  ! The selection process stops when the energy ratio variational/(variational+PT2)
  ! is equal to var_pt2_ratio
  END_DOC

  var_pt2_ratio = correlation_energy_ratio_max
END_PROVIDER

