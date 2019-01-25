! Generates subroutine H_apply_cis
! --------------------------------

BEGIN_SHELL [ /usr/bin/env python2 ]
from generate_h_apply import H_apply
H = H_apply("cis",do_double_exc=False)
print H
END_SHELL

