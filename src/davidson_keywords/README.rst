=================
davidson_keywords
=================

Keywords of the Davidson diagonalization, read from the |EZFIO|
database:

* `n_states_diag`: the number of states diagonalized, a fraction of
  `n_states`
* `threshold_davidson` and `threshold_nonsym_davidson`: the
  convergence thresholds of the symmetric and of the non symmetric
  Davidson algorithm, when `threshold_davidson_from_pt2` is |false|
* `threshold_davidson_from_pt2`: when |true|, the threshold is set to
  :math:`E_{\text{rPT2}} \times` `threshold_davidson_from_pt2`
* `disk_based_davidson`: when |true|, a memory-mapped file may be
  used to store the Davidson vectors if not enough RAM is available
* `distributed_davidson`: when |true|, the distributed algorithm is
  used; set it before running multi-node calculations
* `davidson_sze_max`: the number of micro-iterations before
  re-contracting
* `state_following`, `only_expected_s2`, `csf_based`,
  `n_det_max_full`: behavior options of the diagonalization
