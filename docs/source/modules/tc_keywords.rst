.. _module_tc_keywords: 
 
.. program:: tc_keywords 
 
.. default-role:: option 
 
tc_keywords
===========

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
EZFIO parameters 
---------------- 
 
.. option:: read_rl_eigv
 
    If |true|, read the right/left eigenvectors from ezfio
 
    Default: False
 
.. option:: comp_left_eigv
 
    If |true|, computes also the left-eigenvector
 
    Default: False
 
.. option:: three_body_h_tc
 
    If |true|, three-body terms are included
 
    Default: False
 
.. option:: three_e_3_idx_term
 
    If |true|, the diagonal 3-idx terms of the 3-e interaction are taken
 
    Default: True
 
.. option:: three_e_4_idx_term
 
    If |true|, the off-diagonal 4-idx terms of the 3-e interaction are taken
 
    Default: True
 
.. option:: three_e_5_idx_term
 
    If |true|, the off-diagonal 5-idx terms of the 3-e interaction are taken
 
    Default: True
 
.. option:: pure_three_body_h_tc
 
    If |true|, pure triple excitation three-body terms are included
 
    Default: False
 
.. option:: double_normal_ord
 
    If |true|, contracted double excitation three-body terms are included
 
    Default: False
 
.. option:: noL_standard
 
    If |true|, standard normal-ordering for L (to be used with three_body_h_tc |false|)
 
    Default: True
 
.. option:: core_tc_op
 
    If |true|, takes the usual Hamiltonian for core orbitals (assumed to be doubly occupied)
 
    Default: False
 
.. option:: full_tc_h_solver
 
    If |true|, you diagonalize the full TC H matrix
 
    Default: False
 
.. option:: thresh_it_dav
 
    Thresholds on the energy for iterative Davidson used in TC
 
    Default: 1.e-5
 
.. option:: thrsh_cycle_tc
 
    Thresholds to cycle the integrals with the envelop
 
    Default: 1.e-10
 
.. option:: max_it_dav
 
    nb max of iteration in Davidson used in TC
 
    Default: 1000
 
.. option:: thresh_psi_r
 
    Thresholds on the coefficients of the right-eigenvector. Used for PT2 computation.
 
    Default: 0.000005
 
.. option:: thresh_psi_r_norm
 
    If |true|, you prune the WF to compute the PT1 coef based on the norm. If False, the pruning is done through the amplitude on the right-coefficient.
 
    Default: False
 
.. option:: state_following_tc
 
    If |true|, the states are re-ordered to match the input states
 
    Default: False
 
.. option:: symmetric_fock_tc
 
    If |true|, using F+F^t as Fock TC
 
    Default: False
 
.. option:: selection_tc
 
    if +1: only positive is selected, -1: only negative is selected, :0 both positive and negative
 
    Default: 0
 
.. option:: mu_r_ct
 
    a parameter used to define mu(r)
 
    Default: 1.5
 
.. option:: beta_rho_power
 
    a parameter used to define mu(r)
 
    Default: 0.33333
 
.. option:: zeta_erf_mu_of_r
 
    a parameter used to define mu(r)
 
    Default: 1.
 
.. option:: thr_degen_tc
 
    Threshold to determine if two orbitals are degenerate in TCSCF in order to avoid random quasi orthogonality between the right- and left-eigenvector for the same eigenvalue
 
    Default: 1.e-5
 
.. option:: maxovl_tc
 
    If |true|, maximize the overlap between orthogonalized left- and right eigenvectors
 
    Default: False
 
.. option:: test_cycle_tc
 
    If |true|, the integrals of the three-body jastrow are computed with cycles
 
    Default: False
 
.. option:: thresh_biorthog_diag
 
    Threshold to determine if diagonal elements of the bi-orthogonal condition L.T x R are close enouph to 1
 
    Default: 1.e-6
 
.. option:: thresh_lr_angle
 
    Maximum value of the angle between the couple of left and right orbital for the rotations
 
    Default: 5.0
 
.. option:: thresh_biorthog_nondiag
 
    Threshold to determine if non-diagonal elements of L.T x R are close enouph to 0
 
    Default: 1.e-6
 
.. option:: var_tc
 
    If |true|, use VAR-TC
 
    Default: False
 
.. option:: io_tc_integ
 
    Read/Write integrals int2_grad1_u12_ao, tc_grad_square_ao and tc_grad_and_lapl_ao from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: io_tc_norm_ord
 
    Read/Write normal_two_body_bi_orth from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: only_spin_tc_right
 
    If |true|, only the right part of WF is used to compute spin dens
 
    Default: False
 
.. option:: save_sorted_tc_wf
 
    If |true|, save the bi-ortho wave functions in a sorted way
 
    Default: True
 
.. option:: use_ipp
 
    If |true|, use Manu IPP
 
    Default: True
 
.. option:: tc_grid1_a
 
    size of angular grid over r1: [ 6 | 14 | 26 | 38 | 50 | 74 | 86 | 110 | 146 | 170 | 194 | 230 | 266 | 302 | 350 | 434 | 590 | 770 | 974 | 1202 | 1454 | 1730 | 2030 | 2354 | 2702 | 3074 | 3470 | 3890 | 4334 | 4802 | 5294 | 5810 ]
 
    Default: 50
 
.. option:: tc_grid1_r
 
    size of radial grid over r1
 
    Default: 30
 
.. option:: tc_grid2_a
 
    size of angular grid over r2: [ 6 | 14 | 26 | 38 | 50 | 74 | 86 | 110 | 146 | 170 | 194 | 230 | 266 | 302 | 350 | 434 | 590 | 770 | 974 | 1202 | 1454 | 1730 | 2030 | 2354 | 2702 | 3074 | 3470 | 3890 | 4334 | 4802 | 5294 | 5810 ]
 
    Default: 266
 
.. option:: tc_grid2_r
 
    size of radial grid over r2
 
    Default: 70
 
.. option:: tc_integ_type
 
    approach used to evaluate TC integrals [ analytic | numeric | semi-analytic ]
 
    Default: numeric
 
.. option:: minimize_lr_angles
 
    If |true|, you minimize the angle between the left and right vectors associated to degenerate orbitals
 
    Default: False
 
.. option:: thresh_de_tc_angles
 
    Thresholds on delta E for changing angles between orbitals
 
    Default: 1.e-6
 
.. option:: ao_to_mo_tc_n3
 
    If |true|, memory scale of TC ao -> mo: O(N3)
 
    Default: False
 
.. option:: tc_save_mem_loops
 
    If |true|, use loops to save memory TC
 
    Default: False
 
.. option:: tc_save_mem
 
    If |true|, more calc but less mem
 
    Default: False
 
.. option:: im_thresh_tc
 
    Thresholds on the Imag part of TC energy
 
    Default: 1.e-7
 
.. option:: transpose_two_e_int
 
    If |true|, you duplicate the two-electron TC integrals with the transpose matrix. Acceleates the PT2.
 
    Default: False
