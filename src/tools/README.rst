=====
tools
=====

A collection of small programs to inspect, transform and convert the
results stored in the |EZFIO| database.

* inspecting the wave function and the energies:
  :ref:`print_wf`, :ref:`print_energy`, :ref:`print_var_energy`,
  :ref:`print_dipole`, :ref:`print_ci_vectors`,
  :ref:`print_detweights`, :ref:`print_sorted_wf_coef`
* the Hamiltonian matrix: :ref:`print_hamiltonian` and
  :ref:`print_hmat` print it in the space of the determinants stored
  in the |EZFIO| database, and :ref:`diagonalize_h` extracts the
  :option:`determinants n_states` lowest states of it
* the |MOs|: :ref:`rotate_mos`, :ref:`swap_mos`,
  :ref:`sort_by_fock_energies`, :ref:`save_natorb` (with the
  `save_natorb_no_ov_rot` and `save_natorb_no_ref` variants),
  :ref:`save_ortho_mos`, :ref:`save_one_e_dm`, :ref:`attachement_orb`
* the wave function: :ref:`truncate_wf`, :ref:`sort_wf` and
  :ref:`cas_complete`, which diagonalizes the Hamiltonian in the
  complete active space
* the file formats: :ref:`fcidump` and :ref:`fcidump_pyscf` write an
  `FCIDUMP` file, :ref:`molden` a Molden file,
  :ref:`four_idx_transform` transforms the two-electron integrals
  from the |AO| to the |MO| basis, and :ref:`write_integrals_erf`
  writes the integrals of the
  :math:`\operatorname{erf}(\mu r_{12})/r_{12}` operator
* the analyses: :ref:`analyze_ccsd`, :ref:`check_orthonormality` and
  :ref:`benchmark`

The complete list, with a link to each program page, is in the
Programs section below.
