============================================
Developing new functionals for KS and RS-DFT
============================================

How the functionals are used
============================

A |DFT| (or |RSDFT|) calculation needs two kinds of quantities: the
*exchange/correlation* energy functionals, which give the energy, and the
*exchange/correlation* potentials of the alpha and beta electrons, which
are used to optimize the wave function.

All the |DFT| programs get these quantities from six general providers,
declared in the :ref:`module_dft_one_e` module:

* :c:data:`energy_x` and :c:data:`energy_c`: the *exchange* and
  *correlation* energies,

* :c:data:`potential_x_alpha_ao` and :c:data:`potential_x_beta_ao`: the
  exchange potentials on the |AO| basis,

* :c:data:`potential_c_alpha_ao` and :c:data:`potential_c_beta_ao`: the
  correlation potentials on the |AO| basis.

The potentials on the |MO| basis are then obtained automatically by the
|AO| to |MO| transformation, which is done downstream of these
providers.

Which functional fills the general providers is selected by two
keywords of the :ref:`module_dft_keywords` module, which the programs
such as :ref:`ks_scf` and :ref:`rs_ks_scf` read:

* :option:`dft_keywords exchange_functional`
* :option:`dft_keywords correlation_functional`

The exchange and the correlation functionals are chosen independently,
so that range-separated hybrids (|RSDFT|) are built for instance by
combining a short-range exchange functional with a long-range
correlation functional. The fraction of exact (Hartree-Fock) exchange
is controlled by the :option:`dft_keywords HF_exchange` keyword: the
contribution of the exchange functional is scaled by (1 - HF_exchange).

The selection is done at pre-processing time. The files
:file:`e_xc_general.irp.f` and :file:`pot_general.irp.f` of the
:ref:`module_dft_one_e` module contain a shell block which scans the
:file:`functionals` directory and generates an if/else chain on the
keyword values. As a consequence, adding a functional requires no manual
change of these files: it is enough to add the new providers and to
rebuild.

Adding a new functional
=======================

Create a file :file:`<name>.irp.f` in the
:file:`$QP_ROOT/src/functionals` directory. The current functionals
are :file:`lda.irp.f`, :file:`pbe.irp.f`, :file:`sr_lda.irp.f`,
:file:`sr_pbe.irp.f` and :file:`none.irp.f`. The name must be in
lowercase, and it must be the value given to
:option:`dft_keywords exchange_functional` or
:option:`dft_keywords correlation_functional` when the functional is
used.

Both keywords can take the name, so the file has to define six
providers:

* ``energy_x_<name>`` and ``energy_c_<name>``: the exchange and
  correlation energies (``double precision`` arrays of size
  ``(N_states)``),

* ``potential_x_alpha_ao_<name>``, ``potential_x_beta_ao_<name>``,
  ``potential_c_alpha_ao_<name>``, ``potential_c_beta_ao_<name>``: the
  exchange and correlation potentials of the alpha and beta electrons
  on the |AO| basis (``double precision`` arrays of size
  ``(ao_num,ao_num,N_states)``).

They can be coupled (several ``BEGIN_PROVIDER`` blocks chained with
``&``) as in the existing files. As an example, the file
:file:`sr_pbe.irp.f` defines the energies ``energy_x_sr_pbe`` and
``energy_c_sr_pbe`` together with the four potential arrays on the
|AO| basis, and the file
:file:`none.irp.f` sets the same six quantities to zero.

After adding the file, rebuild the repository from the root
(:command:`ninja`), which regenerates the if/else chain with the new
functional included. The functional is then used by setting
:option:`dft_keywords exchange_functional` and/or
:option:`dft_keywords correlation_functional` to its name.

For a multi-configurational |RSDFT| use, external functional plugins
exist, see `the plugin repository of E. Giner
<https://gitlab.com/eginer/qp_plugins_eginer>`_.


The density used in DFT
=======================

Different ways of building the density
--------------------------------------

The keyword which controls the density is
:option:`density_for_dft density_for_dft`, declared in the
:ref:`module_density_for_dft` module. The accepted values are:

* ``WFT`` (the default): the density is the one of the wave function
  stored in the |EZFIO| directory (the determinants and coefficients),

* ``KS``: the density is the one of a single Slater determinant built
  with the first :c:data:`elec_alpha_num` and
  :c:data:`elec_beta_num` molecular orbitals,

* ``input_density``: a one-body density matrix on the |MO| basis is
  read from the |EZFIO| directory (it is stored in the
  :file:`aux_quantities` group),

* ``input_density_ao``: a one-body density matrix on the |AO| basis is
  read from the |EZFIO| directory and transformed to the |MO| basis,

* ``damping_rs_dft``: a damped mixture of the ``WFT`` and
  ``input_density`` densities, with the damping factor
  :option:`density_for_dft damping_for_rs_dft`,

* ``state_average_dens``: the state-average density of all states,
  stored in the slot of the first state.

The programs :ref:`ks_scf` and :ref:`rs_ks_scf` run with the ``WFT``
choice. If an |MO| basis is already stored in the |EZFIO| directory,
the one-body density matrices are built on this |MO| basis; for the
``KS`` choice, the determinant is built with the current |MO| basis
stored in the |EZFIO| directory.

The resulting density matrices are the providers
:c:data:`one_e_dm_mo_alpha_for_dft`, :c:data:`one_e_dm_mo_beta_for_dft`
and :c:data:`one_e_dm_mo_for_dft` (|MO| basis) and
:c:data:`one_e_dm_alpha_ao_for_dft`, :c:data:`one_e_dm_beta_ao_for_dft`
(|AO| basis) of the :ref:`module_density_for_dft` module.

The density and its gradients on the grid
-----------------------------------------

The values of the density and of its gradients on the |DFT| grid points
are provided by the :ref:`module_dft_utils_in_r` module:

* :c:data:`one_e_dm_and_grad_alpha_in_r` and
  :c:data:`one_e_dm_and_grad_beta_in_r` hold, for each grid point and
  state, the density (line 4) and its three spatial gradients (lines
  1 to 3),

* the :c:func:`density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r`
  subroutine evaluates the density, its gradients, the |AOs| and the
  gradients of the |AOs| at a single point of space.

If the density and its gradients are read through these providers,
their values are coherent with the choice of
:option:`density_for_dft density_for_dft`, because they are built from
the ``for_dft`` density matrices. Symmetrically, the energy and the
potentials are read by every |DFT| program from the general providers
of the :ref:`module_dft_one_e` module, so that the choice of the
:option:`dft_keywords exchange_functional` and
:option:`dft_keywords correlation_functional` keywords is coherent
through the whole calculation.
