.. _module_davidson_keywords: 
 
.. program:: davidson_keywords 
 
.. default-role:: option 
 
=================
davidson_keywords
=================

Keywords used for Davidson algorithms.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: threshold_davidson
 
    Thresholds of Davidson's algorithm if threshold_davidson_from_pt2 is false.
 
    Default: 1.e-10
 
.. option:: threshold_nonsym_davidson
 
    Thresholds of non-symetric Davidson's algorithm
 
    Default: 1.e-10
 
.. option:: davidson_sze_max
 
    Number of micro-iterations before re-contracting
 
    Default: 15
 
.. option:: state_following
 
    If |true|, the states are re-ordered to match the input states
 
    Default: False
 
.. option:: disk_based_davidson
 
    If |true|, a memory-mapped file may be used to store the W and S2 vectors if not enough RAM is availabl
 
    Default: True
 
.. option:: n_states_diag
 
    Controls the number of states to consider during the Davdison diagonalization. The number of states is n_states * n_states_diag
 
    Default: 4
 
.. option:: n_det_max_full
 
    Maximum number of determinants where |H| is fully diagonalized
 
    Default: 1000
 
.. option:: threshold_davidson_from_pt2
 
    Thresholds of Davidson's algorithm is set to E(rPT2)*threshold_davidson_from_pt2
 
    Default: false
 
.. option:: distributed_davidson
 
    If |true|, use the distributed algorithm. If you plan to run multi-node calculations, set this to true before running.
 
    Default: False
 
.. option:: only_expected_s2
 
    If |true|, use filter out all vectors with bad |S^2| values
 
    Default: True
 
.. option:: csf_based
 
    If |true|, use the CSF-based algorithm. Keeps only vectors with expected_s2
 
    Default: False
 
 
Providers 
--------- 
 
.. c:var:: n_states_diag


    File : :file:`davidson_keywords/input.irp.f`

    .. code:: fortran

        integer	:: n_states_diag	


    Number of states to consider during the Davdison diagonalization

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`n_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
       * :c:data:`psi_energy`
       * :c:data:`sxeigenvec`

 
.. c:var:: threshold_davidson_pt2


    File : :file:`davidson_keywords/usef.irp.f`

    .. code:: fortran

        double precision	:: threshold_davidson_pt2	


    Threshold of Davidson's algorithm, using PT2 as a guide

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`threshold_davidson`


