=========
electrons
=========

Number of electrons.

The |EZFIO| parameters :option:`electrons elec_alpha_num` and
:option:`electrons elec_beta_num` give the number of alpha ("up") and
beta ("down") electrons. From them the providers :c:data:`elec_num`
and :c:data:`elec_num_tab` give the total number of electrons and the
(alpha, beta) pair.

Assumptions
===========

* `elec_num` >= 0
* `elec_alpha_num` >= 0
* `elec_beta_num` >= 0
* `elec_alpha_num` >= `elec_beta_num`
