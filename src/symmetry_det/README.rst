============
symmetry_det
============

Determinant and CI-state symmetry: the third layer of the symmetry support,
after the point group (:ref:`module_symmetry`) and the MO labels
(:ref:`module_symmetry_mo`). It labels the determinants of the wave function,
projects the CI states on the irreps, prints the state table, and implements
the :option:`symmetry irrep_states` filter and the irrep-aware state following used by
Davidson.

Needed Modules
==============
* :ref:`module_symmetry_mo` -- point group, character tables, MO irreps
* :ref:`module_determinants` -- ``psi_det``, ``psi_coef``, the sorted-bit
  index
* :ref:`module_utils`

The module deliberately does **not** need :ref:`module_davidson`. The
dependency goes the other way: ``davidson`` needs ``symmetry_det``. This is
what lets the CI programs use the labels without a cycle.

Determinants: families and projectors
=====================================

Two determinants belong to the same *family* when they have the same number
of spin-up and spin-down electrons in every symmetry block. A family is the
set of determinants that a group element can mix, so the projector is built
once per family, not once per determinant.

For an element :math:`g` and two members :math:`D, D'` of a family,

.. math::

   \langle D' | U(g) | D \rangle =
   \prod_b \det\left( X_g[A'_b, A_b] \right)
   \prod_b \det\left( X_g[B'_b, B_b] \right) ,

where :math:`A_b` and :math:`B_b` are the alpha and beta occupied MOs of the
determinant in block :math:`b`, and :math:`X_g` is the MO representation
matrix of :ref:`module_symmetry_mo`. Each block contributes a square
submatrix because the family has the same count in every block.

The projector on irrep :math:`a` is the same formula as everywhere else:

.. math::

   T^a = \frac{d_a}{m_a |G|} \sum_c \chi_a(c) U_c ,

restricted to the family. ``pg_det_irrep(D)`` is the argmax of
:math:`T^a_{DD}` over :math:`a`; it is 0 when
:math:`\max_a T^a_{DD} < 1 -` :option:`symmetry pg_thresh`. ``pg_det_purity(D)`` is
that maximum, the weight of the largest isotypic component of the
determinant.

A determinant that fills whole blocks is an eigenvector of every
:math:`U(g)`, with eigenvalue :math:`\prod_b \det X_g[\mathrm{occ}_b,
\mathrm{occ}_b]`: its family has one member and its purity is 1. A family
of size greater than one exists because some block holds between 1 and
:math:`d_b - 1` electrons of one spin, and every member then shares that
partial occupation. A purity below 1 means that the determinant is a
superposition of irreps: the closed-shell occupation of a degenerate shell
gives 0.5 for one :math:`e^2` determinant of :math:`C_{3v}`, for instance.
The converse deserves a sentence: a one-electron determinant in a
degenerate shell has a family of size :math:`d_b` but purity 1, because it
lies wholly inside that isotypic component; it is a component of the
multiplet, not an eigenvector, and any state built from it carries the
irrep. ``pg_det_block`` is the family size, ``pg_det_family`` the family id,
and ``pg_family_first`` the ``psi_det`` index of the representative the
module enumerates.

A truncated CI space is handled exactly: a family member that is not in
``psi_det`` has coefficient 0, so the missing part of the projector
contributes nothing. The label of a determinant never depends on which
other members happen to be selected.

``pg_det_weights(D, w)`` returns the whole distribution ``w`` of one
determinant, and the report prints it for every determinant without a
label, up to four weights and at most 50 determinants::

    determinants with no definite irrep (probability of each irrep):
        det  family  |F|   purity  composition
          8       8    4  0.793412   A1 0.793412  B2 0.206588

It evaluates only the diagonal of the family projector, :math:`|G| \times
|F|` matrix elements instead of :math:`|G| \times |F|^2`, which is what
makes one call per printed row affordable. The report cross-checks
``max_a w_a`` against ``pg_det_purity``, which comes from the full matrix,
and prints a bookkeeping line when the two differ by more than ``1e-8``.

When the MO basis is not symmetry adapted, the module prints a warning
instead of stopping, and everything above still holds. The blocks are
connected components of ``X_g``, so ``X_g`` is block diagonal on them by
construction, the families stay invariant, and :math:`T^a_{DD}` is the exact
probability that determinant :math:`D` carries irrep :math:`a`.

CI states
=========

The weight of irrep :math:`a` in state :math:`s` is the exact diagonal
matrix element of the projector in the state,

.. math::

   w_a(s) = \sum_{D,D'} c_D(s)\, T^a_{DD'}\, c_{D'}(s) ,

off-diagonal family terms included. ``pg_state_irrep(s)`` is the argmax, or
0 when no weight reaches ``1 - pg_thresh``;
``pg_state_purity(s)`` is the maximum. ``pg_det_mixed_frac(s)`` is the
weight of the state carried by determinants with no definite irrep, and
``pg_state_norm(s) = sum_D c_D(s)^2``.

``sum_a w_a(s)`` equals ``pg_state_norm(s)`` for every state. A deficit is
only a bookkeeping statement about ``psi_coef`` (it is not normalized in the
determinant basis); it never means that the CI space is not group
invariant. The sum rule is reported, not gated: the exact mathematical
content is already in the projector construction.

The state table
---------------

``pg_print_state_irreps(energy, s2)`` is called by :ref:`cis`, :ref:`cisd`
and :ref:`fci` after diagonalization. It prints one line per state::

    symmetry: state irreps: group C2v, 4 irreps, N_det = 679
    symmetry:    state      energy (Ha)       S^2  irrep   purity
    symmetry:         1    -76.1140581951    0.0000  A1   1.000

and nothing when the detection failed or was not requested. The
:math:`S^2` column is the caller's value, so a symmetry label and a spin
label are always read together. There is no printing switch: the table is
part of a CI run whose labels are usable.

The table uses the exact weights above, with the projectors for the families
of size greater than one. For a pure determinant the diagonal is a Kronecker
delta and the contribution is taken from the label, which is why the cost is
one pass over the determinants per state plus one projector per multi-member
family.

When the MO basis is not symmetry adapted the table gains two notice lines
before the header::

    symmetry: the MO basis is not symmetry adapted: a state weight is its
    symmetry: exact irrep probability for the stored CI space (not a guess).

and a state without a definite irrep ends its line with up to four
``name weight`` pairs, the irreps that carry its weight.

The weights stay exact for the stored determinant space. ``pg_mo_block`` is
built from ``X_g`` itself, so ``X_g`` is block diagonal on the blocks
whether or not the MOs are eigenvectors of the projectors: a family keeps its
per-block occupations under every group element, and the projectors close on
the stored space. What is lost is the right to call a weight a label, which
is why state following and the :option:`irrep_states` filter are off on such
a basis (``pg_labels_usable`` tests ``pg_mo_labels_ok``).

Fast labels for Davidson
------------------------

Davidson cannot afford the exact family projectors at every iteration. It
uses ``pg_det_state_labels``, which sums ``coef(i,k)**2`` over the
determinants of each irrep and over the unlabelled determinants, and returns
a relative purity and a label. Everything is a ratio because Davidson hands
in vectors that are not normalized. A state whose weight sits on
multi-member families gets label 0 and is not constrained; the exact
``pg_state_irrep_weight`` remains the value the report prints.

The labels are used to follow states: the label of the current vector is
compared with the label of the previous wave function, so state following
does not swap two states of different irreps.

The ``irrep_states`` filter
===========================

:option:`symmetry irrep_states` is a list of irrep names. When it is non-empty and
the symmetry machinery is usable, the diagonalization keeps only the states
whose known irrep is in the list. The semantics are deliberately narrow:

* it is applied **after** the :math:`S^2` selection and, in Davidson, after
  the state-following attribution. It never widens a selection; it is a
  filter on the states that already survived;
* a state with no definite irrep (label 0) is always kept: a filter must not
  drop a state it cannot judge;
* the kept states are moved to the front of the CI vector set, preserving
  their energy order;
* when fewer states than requested are allowed, the run completes and says
  so: ``symmetry: irrep filter: only N state(s) in irrep_states for n_states
  = M``. This message is printed once per diagonalization call, so a CI run
  with several calls shows it several times;
* every state whose known irrep is rejected prints
  ``symmetry: irrep filter: state <k> is <name>, not in irrep_states``;
* CSF-based Davidson (``csf_based = true``) cannot be filtered. It prints
  ``symmetry: irrep_states is ignored when csf_based is true`` and runs
  unchanged;
* an unknown irrep name prints the available names and stops: malformed
  input is the one place where stopping is allowed.

On the Lapack path (small determinant spaces, ``N_det`` at most
``n_det_max_full``), the filter is a *stable permutation* of the list that
the :math:`S^2` selection already produced. That list is "states of the
expected :math:`S^2`, then fill states". The filter promotes the allowed
fill states too, so asking for ``A1`` with ``n_states = 3`` can return two
singlets and a triplet: it is the third ``A1`` state of the computed slots,
not necessarily the third singlet. The energies are unchanged; only the
order and the count of kept states are affected.

State following and the filter compose as "keep the states that follow,
then move the allowed ones forward". The filter therefore acts on the states
that passed the :math:`S^2` test, not on the whole subspace; when the
requested irreps are not among them, raise
:option:`davidson_keywords n_states_diag` (the multiplier of
:option:`determinants n_states` that sets the number of computed roots) or
:option:`determinants n_states`, and re-run.

Report program
==============

``qp_run symmetry_det DIR`` prints the families, the determinant table (the
individual lines are hidden above 200 determinants), the per-irrep weights
of every state and the gates D1-D6. It exits 1 if a gate fails, after every
number has been printed, so it can be chained in a test script.

The program reads the wave function stored in the EZFIO. It is meaningful
after a CI calculation; after an SCF only, it labels the Hartree-Fock
determinant.

Gates
=====

.. list-table::
   :header-rows: 1

   * - Gate
     - Test
   * - D1
     - :math:`\sum_a T^a = 1` on every family (validates the characters,
       :math:`m_a` and the family structure at once)
   * - D2
     - :math:`(T^a)^2 = T^a` and :math:`T^a` Hermitian, every family,
       every irrep
   * - D3
     - :math:`-10^{-8} \le` purity :math:`\le 1 + 10^{-8}` for every
       determinant
   * - D4
     - :math:`|\sum_a w_a(s) - \sum_D c_D(s)^2| \le 10^{-8}` (the
       bookkeeping sum rule)
   * - D5
     - every labelled state has :math:`w_{\max} \ge 1 - 10^{-6}`; the mixed
       determinant fraction is reported
   * - D6
     - the dominant determinant of state 1 has irrep 1 for a closed-shell
       Hartree-Fock wave function, when that determinant is pure

Families of size one are skipped by D1 and D2: there :math:`T^a` is a
Kronecker delta by construction, and the label already comes from it.

Physics note
============

An open degenerate shell has no one-determinant irrep, and neither has a
closed-shell occupation of a partially filled shell: with one electron in a
two-dimensional ``e`` block, the determinant is a 50/50 mixture of two
irreps. That is the SCF solution breaking the point group, not a labelling
error. The determinant is reported with a purity below 1 and no label; the
state weights carry the symmetry instead. On an Abelian group every
determinant is pure, whatever its open shells.

Limits
======

* a determinant family is enumerated only up to ``pg_family_cap`` (1000)
  members;
* one symmetry block may hold at most ``pg_block_cap`` (30) MOs, because the
  occupation of a block has to fit in one integer bit mask;
* the projector refuses a block with more than 12 occupied MOs in a family;
* when a limit is hit, ``pg_label_ierr`` and ``pg_label_msg`` say which
  determinant was refused, and the CI programs run without symmetry labels.
  Nothing hangs and nothing is silently truncated.

The first two limits are unreachable on a symmetry-adapted basis: a block is
then one irrep, with at most five MOs. They bound the mixed-basis path, where
a block is a connected component of ``X_g`` and can be as wide as the number
of MOs that mixing touches.
