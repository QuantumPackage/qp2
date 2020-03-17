! Generates subroutine H_apply_cisd
! ----------------------------------

BEGIN_SHELL [ /usr/bin/env python3 ]
from generate_h_apply import H_apply
H = H_apply("cisd",do_double_exc=True)
print(H)

H = H_apply("cisd_sym",do_double_exc=True)
H.filter_only_connected_to_hf()
print(H)
END_SHELL

