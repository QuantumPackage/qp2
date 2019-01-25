 BEGIN_PROVIDER [ integer, n_pt_max_integrals ]
&BEGIN_PROVIDER [ integer, n_pt_max_i_x]
 implicit none
 BEGIN_DOC
! Number of points used in the numerical integrations.
 END_DOC
 integer :: n_pt_sup
 integer :: prim_power_l_max
 include 'utils/constants.include.F'
 prim_power_l_max = maxval(ao_power)
 n_pt_max_integrals = 24 * prim_power_l_max + 4
 n_pt_max_i_x = 8 * prim_power_l_max
 ASSERT (n_pt_max_i_x-1 <= max_dim)
 if (n_pt_max_i_x-1 > max_dim) then
   print *, 'Increase max_dim in utils/constants.include.F to ', n_pt_max_i_x-1
   stop 1
 endif
END_PROVIDER

