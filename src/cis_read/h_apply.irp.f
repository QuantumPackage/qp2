! Generates subroutine H_apply_cis
! --------------------------------

BEGIN_SHELL [ /usr/bin/env python3 ]
from generate_h_apply import H_apply
H = H_apply("cis",do_double_exc=False)
print(H)

H = H_apply("cis_sym",do_double_exc=False)
H.filter_only_connected_to_hf()
print(H)

END_SHELL

