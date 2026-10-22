.. _cipsi:

=====
CIPSI
=====

The :ref:`fci` program runs selected configuration interaction in the full
|FCI| space with the |CIPSI| algorithm. The theory is described in
:doc:`../intro/selected_ci`; this page describes the practical workflow.
There is no program called ``cipsi``: the entry point is ``fci``.


Preparing the calculation
-------------------------

|CIPSI| needs a set of |MOs| and the two-electron integrals on this basis in
the |EZFIO| database. A Hartree-Fock calculation provides both:

.. code:: bash

    qp run scf

The :ref:`qp_set_frozen_core` command freezes the core orbitals, and
:ref:`qp_set_mo_class` gives full control over the orbital classes:

.. code:: bash

    qp set_frozen_core
    qp set_mo_class --core="[1-2]" --act="[3-20]"

All the determinants that can be built in the space of the active orbitals are
reachable by the selection. Orbitals of the *core* class are always doubly
occupied, and orbitals of the *deleted* class are never occupied.

If the two-electron integrals have been transformed and stored on disk by a
previous run, they are read back instead of being recomputed. The
:ref:`four_idx_transform` program performs this transformation explicitly,
which is useful when a long calculation is restarted several times.


Running the calculation
-----------------------

.. code:: bash

    qp run fci | tee h2o.fci.out

The program starts from the Hartree-Fock determinant, or from the wave
function stored in the |EZFIO| database when
:option:`determinants read_wf` is |true|. Each iteration performs the following
steps:

#. The Hamiltonian is diagonalized in the current space, with the Davidson
   algorithm, or with LAPACK when the number of determinants is at or below
   ``n_det_max_full``.
#. The second-order perturbative contribution of the external determinants is
   estimated (by sampling the external space in the default path, or by
   evaluating the candidates generated from the generator space in the
   deterministic path) and used to select the determinants that enter the
   space at the next iteration.
#. The new wave function is saved in the |EZFIO| database, together with the
   variational energy and the |PT2| correction.

The program prints a summary block per iteration with the number of
determinants, the variational energy, the |PT2| correction, the variance, and
the statistical error on the |PT2|. The last eight iterations are kept in
memory to extrapolate the |FCI| energy with
:math:`E_\text{FCI} = E - \alpha E_\text{PT2}`; the extrapolated value is
printed at each iteration.

With :option:`determinants s2_eig` set to |true|, which is the default, the
wave function is constrained to be an eigenfunction of |S^2| with the value
:option:`determinants expected_s2`. Before the first diagonalization, and
after any truncation of a restarted wave function, the program completes each
configuration with all the determinants that can be built from it, so that
each configuration is invariant under |S^2| and the Hamiltonian eigenstates
found within it are spin eigenstates. In the full-matrix path of the
diagonalization, the program diagonalizes
:math:`H + 0.1\,(S^2 - S^2_\text{target})` instead of :math:`H`, which
penalizes the states whose |S^2| deviates from the target value. The
:doc:`excited_states` page describes the ``only_expected_s2`` filter and the
``state_following`` option.


Deterministic and stochastic perturbation theory
------------------------------------------------

The default path is the semi-stochastic algorithm of
:cite:`Garniron_2017b`, selected by :option:`perturbation do_pt2` set to
|true|. The selection and the |PT2| correction are computed together by
sampling the external space, so the |PT2| is unbiased and carries a statistical
error bar. The sampling stops when the relative error on the |PT2|,
:math:`\sigma_\text{PT2}/|E_\text{PT2}|`, falls below
:option:`perturbation pt2_relative_error` (default :math:`2\times10^{-3}`).

With :option:`perturbation do_pt2` set to |false|, the program calls the
``run_cipsi`` path instead: the candidates generated from the generator space
are scanned explicitly to select the determinants, and no stochastic |PT2|
sampling is performed. The |PT2|
printed in this mode comes from the selection step and is biased, because it
uses a truncated wave function. In the small-wave-function regime this path
is cheaper than the hybrid algorithm, since every generated candidate is
evaluated exactly once, without the overhead of the sampling machinery. It is
appropriate when the goal is to obtain a wave function, for instance as a
trial wave function for |QMC|, rather than an accurate energy.

The zeroth-order Hamiltonian of the perturbative correction is chosen with
:option:`perturbation h0_type`. The default, ``EN``, uses the variational
energy of the state in the Epstein-Nesbet denominator. ``HF`` replaces this
denominator by the diagonal energy of the determinant with the largest
coefficient. ``CFG`` uses configuration energies in the selection step, with
the same denominator as ``EN`` in the |PT2|.


Stopping criteria
-----------------

The main loop continues while all the following conditions hold: the number of
determinants is below :option:`determinants n_det_max`, the largest absolute
|PT2| over the states is above :option:`perturbation pt2_max`, the largest
variance is above :option:`perturbation variance_max`, and the correlation
energy ratio is at or below
:option:`perturbation correlation_energy_ratio_max`. The calculation stops as
soon as one of these conditions fails.

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Parameter
     - Default
     - Stops the calculation when
   * - ``n_det_max``
     - 1000000
     - the space reaches this size
   * - ``pt2_max``
     - :math:`10^{-4}`
     - max :math:`|E_\text{PT2}|` is at or below this value
   * - ``variance_max``
     - 0.0
     - max variance is at or below this value (0 disables)
   * - ``correlation_energy_ratio_max``
     - 1.0
     - the ratio exceeds this value (1 disables)

The correlation energy ratio is
:math:`(E_\text{CI}-E_\text{HF})/(E_\text{CI}+E_\text{PT2}-E_\text{HF})`.
Its default value of 1 makes the criterion inactive; setting it to a value
below 1, for example 0.99, stops the calculation at a fixed fraction of the
correlation energy and is convenient to compare molecules of different sizes
at a comparable accuracy.

A limit on the size of the space is not a convergence criterion. When the
calculation stops because of ``n_det_max``, one extra iteration computes the
|PT2| correction for the final wave function, and the printed extrapolated
energy is the best estimate of the |FCI| energy for this space.


Selection parameters
--------------------

:option:`determinants selection_factor` sets the number of determinants passed
to the selection: the target is
:math:`\sqrt{N_\text{states}}\, N_\text{det}\times\text{factor}`, with a
minimum of ``n_states_diag``. With the default factor
of 1, the size of the space roughly doubles at each iteration.

:option:`determinants threshold_generators` defines the generator space. The
determinants of the current wave function are sorted by decreasing weight,
and only those carrying this fraction of the squared norm generate the
external determinants. The default of 0.999 keeps almost the whole wave
function as generator; lowering it restricts the selection to the largest
contributions and reduces the cost of an iteration.

The selectors are the determinants on which the matrix elements
:math:`\langle i | H | \psi \rangle` of the perturbative estimate are
computed. Their threshold, ``threshold_selectors``, is the square root of
``threshold_generators``. Because both thresholds are applied to the squared
norm, the selector threshold is larger, and the selector space is never
smaller than the generator space.


Restarting a calculation
------------------------

The wave function is saved in the |EZFIO| database after each iteration, so an
interrupted calculation is restarted from the last saved iteration by setting

.. code:: bash

    qp set determinants read_wf true
    qp run fci | tee h2o.fci.restart.out

The :ref:`pt2` program recomputes the |PT2| correction of a stored wave
function without performing any selection. It reads the wave function,
samples the |PT2| with the stochastic scheme until the relative error
reaches :option:`perturbation pt2_relative_error`, and writes the energy
and its |PT2| correction back to the database:

.. code:: bash

    qp run pt2

The convergence history of a run is in its standard output: every
iteration prints a summary block with the number of determinants, the
variational energy, the |PT2| correction and the corrected energy, and
the end of the output adds the extrapolated values.


Multi-state calculations
------------------------

For excited states, :option:`determinants n_states` is set to the number of
requested states and the selection uses the largest contribution over the
states. :doc:`excited_states` describes the multi-state workflow, the use of
:option:`determinants s2_eig` and the state-following algorithm.


.. seealso::

    The documentation of the :ref:`module_cipsi` and
    :ref:`module_perturbation` modules, and of the :ref:`fci` and :ref:`pt2`
    programs.
