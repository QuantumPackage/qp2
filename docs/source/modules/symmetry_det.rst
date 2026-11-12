.. _module_symmetry_det: 
 
.. program:: symmetry_det 
 
.. default-role:: option 
 
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

Three consumers read three different objects of that one formula, and that
is what keeps the cost down:

* the **determinant labels** need only the diagonal :math:`T^a_{DD}` of a
  family: :math:`|G| \times |F|` matrix elements instead of
  :math:`|G| \times |F|^2`, and no :math:`|F|^2` projector storage;
* the **state weights** need the full :math:`U_g` matrix, but only the
  scalars :math:`q(g,s) = c^T U_g c`, one number per element and per state.
  The matrix of an element is built once and contracted with every column,
  and only then multiplied by the projector weight, so the irrep loop sits
  outside the pair loop: :math:`n_a n_{st}` operations per element, not
  :math:`n_a n_{st} |F|^2`;
* the **gates** are the only place that builds the full :math:`T^a`.

Families are enumerated once, at their representative
(``pg_family_first(f)``, the smallest ``psi_det`` index of the family);
every other member is copied from that pass. Inside a family, the
classification of the sectors (one sector is one spin and one block) does
not depend on the group element nor on the pair: an empty sector is
nothing, a full sector has the same block determinant for every member, and
only a partly occupied sector costs per-pair work. The occupied MO lists of
the partial sectors are computed once per family, and the block determinant
of a full sector is read from a table built once per (block, element).
``pg_xg(g,:,:)`` stores the matrix of one element as a contiguous plane,
which is the layout those reads want, and it is the only place the MO
matrices are kept.

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
makes one call per printed row affordable. This routine is the reference
oracle: it is the only remaining caller of ``pg_elem_me``, the literal
per-pair kernel, and the label and state-weight passes never call it. The
report cross-checks ``max_a w_a`` against ``pg_det_purity``, which comes
from the full matrix, and prints a bookkeeping line when the two differ by
more than ``1e-8``.

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
``pg_state_dom(s)`` is the same argmax without that condition -- the most
probable irrep, printed when no definite one exists;
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

    State irreps: group C2v, 4 irreps, N_det = 679

     ========= ================= ========= ============ ======
         State       Energy (Ha)       S^2 Irrep        Purity
     ========= ================= ========= ============ ======
             1    -76.1140581951    0.0000 A1            1.000
     ========= ================= ========= ============ ======

and nothing when the detection failed or was not requested. The
:math:`S^2` column is the caller's value, so a symmetry label and a spin
label are always read together. There is no printing switch: the table is
part of a CI run whose labels are usable. The irrep is left-justified in a
fixed 12-character column, so short labels line up with long ones. When a
state has no *definite* irrep, the column shows its most probable one
(``pg_state_dom``), and the ``Purity`` column next to it says how far the
state is from a pure symmetry: a truncated CIPSI state prints
``Sigma_g+`` with purity ``0.998`` instead of ``(none)``. ``(none)`` remains
for a state whose irrep weights are all negligible.

The table uses the exact weights above, with the projectors for the families
of size greater than one. For a pure determinant the diagonal is a Kronecker
delta and the contribution is taken from the label, so those determinants
cost one pass over the coefficients per state. For a multi-member family the
matrix of one group element is built once and contracted with every state
column before the irrep loop, so the extra cost is one :math:`|F|^2` matrix
per element and family plus one scalar product per element and state, not
one projector per state.

When the MO basis is not symmetry adapted the table gains two notice lines
before the header::

    The MO basis is not symmetry adapted: a state weight is its
    exact irrep probability for the stored CI space (not a guess).

and a state without a definite irrep shows its most probable one in the
irrep column and gets one indented ``composition:`` line with up to four
``name weight`` pairs, the irreps that carry its weight. The composition is
not appended to the table row: it would run past the closing rule.

The weights stay exact for the stored determinant space. ``pg_mo_block`` is
built from ``X_g`` itself, so ``X_g`` is block diagonal on the blocks
whether or not the MOs are eigenvectors of the projectors: a family keeps its
per-block occupations under every group element, and the projectors close on
the stored space. What is lost is the right to call a weight a label, which
is why state following and the :option:`irrep_states` filter are off on such
a basis (``pg_labels_usable`` tests ``pg_mo_labels_ok``).

Exact zeros in the saved vectors
--------------------------------

The CI diagonalization leaves a coefficient that the MO point group forbids
at roundoff (:math:`10^{-12}` to :math:`10^{-15}`), not at zero, and every
later quantity that reads the wave function inherits that noise: the 1-RDM,
the natural orbitals, a perturbative correction. :ref:`cis`, :ref:`cisd` and
:ref:`fci` call ``pg_zero_state_irreps()`` before their last save:
``psi_coef(i,s)`` is set to zero when the determinant ``i`` and the state
``s`` both carry a definite irrep and the two differ. A label 0 means
"unknown" -- a determinant of a degenerate family, a state that mixes
irreps -- and such a row or a column is never touched, because a best guess
would destroy amplitude the labels do not justify.

The rule is the one that already produced the labels: same projectors, same
``pg_thresh``, so no new |EZFIO| option and no new tolerance. Nothing
printed changes -- what is removed is below the printed precision -- but the
stored vector carries exact zeros, and the empty blocks of a natural-orbital
run then have occupations that are exactly ``0.0``. A plain run
(``auto_symmetry = false``) never calls it.

State irreps in a CIPSI run
---------------------------

A CIPSI run prints the state labels at three points beyond the final table:
the ``# Irrep`` row of the "Summary at N_det" table, the per-state headings
of the extrapolation block, and one ``Irrep of state`` line after every
diagonalization. All three read the **exact** weights, through
``pg_state_irreps_psi`` (labels of :c:data:`psi_coef`) and
``pg_state_irreps_vec`` (labels of the vectors that were just diagonalized).
A state can therefore never be selected on one definition of its symmetry
and printed on another, and the end-of-run table is the reference the other
three agree with.

Each entry point is a thin dispatcher that tests ``pg_labels_usable()`` and
nothing else, then calls a worker that names the providers. The split is
required by IRPF90, which hoists the build of every provider named by a
procedure to that procedure's entry: a guard and a ``PROVIDE`` cannot share
a body. With ``auto_symmetry = false``, or when the detection fails or the
MO basis is not symmetry adapted, the dispatcher returns before a single
provider of the chain is built and the run prints exactly what it printed
before the feature existed.

``# Irrep`` is a row of the summary table, between the second rule and the
``# E`` row, with one 30-character field per state (the width of the numeric
columns above and below it), left-justified like the table's other labels::

    # ============ ============================= =============================
                         State      1                  State      2
    # ============ ============================= =============================
    # Irrep       A1                            B2
    # E             -76.11405819                  -75.47913080

It prints the label, and a purity row directly below it when at least one
state has no definite irrep: a truncated CIPSI state whose purity is below
:math:`1 - \text{pg\_thresh}` still shows its most probable irrep, with the
weight fraction that supports it::

    # ============ ============================= =============================
                         State      1                  State      2
    # ============ ============================= =============================
    # Irrep       A1                            B2
    # Purity              0.9977                        0.8109

The purity row is omitted when every state is definite, so a pure run keeps
the historical summary. The extrapolation block appends the label, and its
purity when the label is not definite, to the heading:
``State 2   Irrep: B2 (purity 0.8109)``.

After each diagonalization the :c:data:`CI_energy` provider prints
``Irrep of state 2     Sigma_g+ (purity 0.9921)``: the most probable irrep
with the purity, because that is where a state may still be contaminated by
the selection. The labels of that line come from ``CI_eigenvectors``, not
from :c:data:`psi_coef`: inside the provider :c:data:`psi_coef` still holds
the previous iteration's vectors and is overwritten only once the provider
has been built. A CIPSI state with a few-percent contamination therefore
shows the same dominant irrep and the same purity at every stage of the run
and in the final table.

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
 
 
 
Programs 
-------- 
 
 * :ref:`symmetry_det` 
 
Providers 
--------- 
 
.. c:var:: int pg_block_cap


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer	:: pg_block_cap	


    Widest symmetry block the family enumeration can represent. The
    occupation of a block is a bit mask, so a block wider than the bits of
    a default integer cannot be unranked. Blocks are irrep spaces
    (``d_a`` <= 5) unless the MO basis is not symmetry adapted, which this
    module now accepts: the cap is what keeps that acceptance bounded.

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`

 
.. c:var:: int pg_block_dim


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_block_start	(pg_n_mo_block)
        integer, allocatable	:: pg_block_mo	(mo_num)
        integer, allocatable	:: pg_block_dim	(pg_n_mo_block)


    ``pg_block_mo(pg_block_start(b) : pg_block_start(b)+pg_block_dim(b)-1)``
    holds the MOs of symmetry block ``b`` in ascending MO order.
    
    The MOs of one block are NOT contiguous in the MO list (the same
    irrep can appear at several energies: the A1g MOs of benzene are
    spread over the whole spectrum), which is why the family enumeration
    of the family goes through this indirection and never through
    ``pg_mo_block(i)`` arithmetic.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_block_mo


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_block_start	(pg_n_mo_block)
        integer, allocatable	:: pg_block_mo	(mo_num)
        integer, allocatable	:: pg_block_dim	(pg_n_mo_block)


    ``pg_block_mo(pg_block_start(b) : pg_block_start(b)+pg_block_dim(b)-1)``
    holds the MOs of symmetry block ``b`` in ascending MO order.
    
    The MOs of one block are NOT contiguous in the MO list (the same
    irrep can appear at several energies: the A1g MOs of benzene are
    spread over the whole spectrum), which is why the family enumeration
    of the family goes through this indirection and never through
    ``pg_mo_block(i)`` arithmetic.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_block_start


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_block_start	(pg_n_mo_block)
        integer, allocatable	:: pg_block_mo	(mo_num)
        integer, allocatable	:: pg_block_dim	(pg_n_mo_block)


    ``pg_block_mo(pg_block_start(b) : pg_block_start(b)+pg_block_dim(b)-1)``
    holds the MOs of symmetry block ``b`` in ascending MO order.
    
    The MOs of one block are NOT contiguous in the MO list (the same
    irrep can appear at several energies: the A1g MOs of benzene are
    spread over the whole spectrum), which is why the family enumeration
    of the family goes through this indirection and never through
    ``pg_mo_block(i)`` arithmetic.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_col_diag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_det_block


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label.
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights use
    this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A family is enumerated once, at
    ``pg_family_first(f)``, the smallest psi_det index of the family: the
    loop skips every determinant that this representative already
    labelled (``pg_det_family /= 0``), so a member costs nothing beyond
    the copy of its label. The labelling itself is
    :c:func:`pg_family_labels`, whose working set is ``np*pg_n_irrep``
    instead of ``np^2*pg_n_irrep``.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_cap`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_mo_mixed`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_det_family


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label.
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights use
    this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A family is enumerated once, at
    ``pg_family_first(f)``, the smallest psi_det index of the family: the
    loop skips every determinant that this representative already
    labelled (``pg_det_family /= 0``), so a member costs nothing beyond
    the copy of its label. The labelling itself is
    :c:func:`pg_family_labels`, whose working set is ``np*pg_n_irrep``
    instead of ``np^2*pg_n_irrep``.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_cap`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_mo_mixed`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_det_irrep


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label.
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights use
    this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A family is enumerated once, at
    ``pg_family_first(f)``, the smallest psi_det index of the family: the
    loop skips every determinant that this representative already
    labelled (``pg_det_family /= 0``), so a member costs nothing beyond
    the copy of its label. The labelling itself is
    :c:func:`pg_family_labels`, whose working set is ``np*pg_n_irrep``
    instead of ``np^2*pg_n_irrep``.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_cap`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_mo_mixed`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: double pg_det_mixed_frac


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        integer, allocatable	:: pg_state_dom	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the MO projector
    in the state, not an occupation-number guess. Computed by
    :c:func:`pg_state_weight_matrix`, which the irrep filter also calls.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess must not drive the filter or the zeroing.
    ``pg_state_dom(s)`` = the same argmax without that condition: the most
    probable irrep, 0 only when every weight is below :c:data:`pg_thresh`.
    It is for the printed report, which names the dominant irrep together
    with its purity instead of saying "none"; nothing may constrain a state
    on it.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_n_col`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_proj_w`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: double pg_det_purity


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label.
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights use
    this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A family is enumerated once, at
    ``pg_family_first(f)``, the smallest psi_det index of the family: the
    loop skips every determinant that this representative already
    labelled (``pg_det_family /= 0``), so a member costs nothing beyond
    the copy of its label. The labelling itself is
    :c:func:`pg_family_labels`, whose working set is ``np*pg_n_irrep``
    instead of ``np^2*pg_n_irrep``.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_cap`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_mo_mixed`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: double pg_detfull


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_detfull	(pg_n_mo_block,pg_n_elem)


    ``pg_detfull(b,g)`` = det( X_g[occ_b, occ_b] ) for a block ``b`` that is
    entirely occupied: the block determinant of element ``g``, the factor
    every determinant of the family shares.
    
    A sector that is full has a block determinant that does not depend on
    the determinant: it is read once per (block, element) here instead of
    being recomputed for every pair. The family descriptor of
    :c:func:`pg_fam_setup` sends a block with more than
    :c:data:`pg_det_nmax` occupied MOs to the overflow branch before any
    of this is read, so the 0.d0 written for such a block is never used:
    it keeps the array defined, nothing more.
    
    Cost ``nblock*|G|`` small determinants, built once.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_xg`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_elem_col


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_elem_col	(pg_n_elem)


    Inverse of :c:data:`pg_col_member`: ``pg_elem_col(g)`` is the file column that
    element ``g`` belongs to. :c:data:`pg_char` is indexed by column and both
    projectors sum over elements, so this is the map between the two.
    
    A column can stand for several geometric classes (T, Td, Oh, Ih have
    classes of 4, 6, 8, ... elements); the character is constant on the
    column, which is why the MO projectors sum over columns, not over
    classes.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_proj_w`

 
.. c:var:: int pg_family_cap


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer	:: pg_family_cap	


    Largest family the module agrees to enumerate. A family has
    ``prod_b C(d_b,n_b)`` members per spin, so a determinant with several
    partially occupied degenerate shells grows combinatorially. Above
    this size the module stops and says which determinant it was,
    instead of hanging or silently truncating the projector.

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`

 
.. c:var:: int pg_family_first


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label.
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights use
    this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A family is enumerated once, at
    ``pg_family_first(f)``, the smallest psi_det index of the family: the
    loop skips every determinant that this representative already
    labelled (``pg_det_family /= 0``), so a member costs nothing beyond
    the copy of its label. The labelling itself is
    :c:func:`pg_family_labels`, whose working set is ``np*pg_n_irrep``
    instead of ``np^2*pg_n_irrep``.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_cap`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_mo_mixed`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_label_ierr


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer	:: pg_label_ierr	
        character*(512)	:: pg_label_msg	


    0 when the last labelling attempt succeeded; otherwise
    
    ``1``
        there is no determinant to label
    
    ``2``
        there is no MO or no symmetry block
    
    ``3``
        one symmetry block is wider than :c:data:`pg_block_cap` MOs
    
    ``4``
        a determinant family is larger than pg_family_cap
    
    ``5``
        the family enumeration is inconsistent
    
    ``6``
        one symmetry block holds more than 12 occupied MOs
    
    ``7``
        the per-sector radices do not decompose the family (internal
        inconsistency between :c:func:`pg_fam_setup` and the odometer of
        :c:func:`pg_family_fill`)
    
    :c:data:`pg_label_msg` is the one-sentence reason. Reporting is left to the
    caller: a provider that stops cannot be reported on.

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: char pg_label_msg


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer	:: pg_label_ierr	
        character*(512)	:: pg_label_msg	


    0 when the last labelling attempt succeeded; otherwise
    
    ``1``
        there is no determinant to label
    
    ``2``
        there is no MO or no symmetry block
    
    ``3``
        one symmetry block is wider than :c:data:`pg_block_cap` MOs
    
    ``4``
        a determinant family is larger than pg_family_cap
    
    ``5``
        the family enumeration is inconsistent
    
    ``6``
        one symmetry block holds more than 12 occupied MOs
    
    ``7``
        the per-sector radices do not decompose the family (internal
        inconsistency between :c:func:`pg_fam_setup` and the odometer of
        :c:func:`pg_family_fill`)
    
    :c:data:`pg_label_msg` is the one-sentence reason. Reporting is left to the
    caller: a provider that stops cannot be reported on.

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_n_family


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_det_irrep	(N_det)
        double precision, allocatable	:: pg_det_purity	(N_det)
        integer, allocatable	:: pg_det_block	(N_det)
        integer, allocatable	:: pg_det_family	(N_det)
        integer	:: pg_n_family	
        integer, allocatable	:: pg_family_first	(N_det)


    ``pg_det_irrep(D)`` = argmax_a T^a_{DD}, or 0 when max_a T^a_{DD} <
    1 - pg_thresh: a determinant with one electron in a degenerate block
    (:math:`|e_x>`, :math:`|e_y>`) is not an eigenfunction of the projectors and must not
    get a best-guess label.
    
    ``pg_det_purity(D)`` = max_a T^a_{DD} = sum_{D'} :math:`(T^a_{D'D})^2` at the
    argmax, 1 for a pure determinant and 1/d_a for one member of an open
    degenerate shell.
    
    ``pg_det_block(D)`` = :math:`|F|`, the dimension of the family of D;
    ``pg_det_family(D)`` = its family id, in order of first appearance;
    ``pg_family_first(f)`` = the psi_det index of the first member of
    family f, i.e. the representative the other modules enumerate.
    
    A determinant is pure if and only if its family has size 1: a family
    of size > 1 exists only because a block has 1 <= n <= d_b-1, and then
    every member shares that partial occupation. The state weights use
    this to skip the projectors for pure determinants.
    
    :c:data:`pg_n_family` <= N_det. A family is enumerated once, at
    ``pg_family_first(f)``, the smallest psi_det index of the family: the
    loop skips every determinant that this representative already
    labelled (``pg_det_family /= 0``), so a member costs nothing beyond
    the copy of its label. The labelling itself is
    :c:func:`pg_family_labels`, whose working set is ``np*pg_n_irrep``
    instead of ``np^2*pg_n_irrep``.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`pg_block_cap`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_mo_mixed`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: double pg_proj_w


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_proj_w	(pg_n_irrep,pg_n_elem)


    ``pg_proj_w(a,g)`` = :math:`(d_a/(m_a |G|)) \chi_a(col(g))`, the weight
    of element ``g`` in the projector on irrep ``a``.
    
    This is the ONLY character lookup of the determinant and state
    kernels: the pair loop of :c:func:`pg_family_labels`,
    :c:func:`pg_family_vg` and :c:func:`pg_family_wadd` reads it once per
    (element, irrep) instead of going back to :c:data:`pg_char` through
    :c:data:`pg_elem_col` for every determinant pair.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_elem_col`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_proj_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_state_dom


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        integer, allocatable	:: pg_state_dom	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the MO projector
    in the state, not an occupation-number guess. Computed by
    :c:func:`pg_state_weight_matrix`, which the irrep filter also calls.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess must not drive the filter or the zeroing.
    ``pg_state_dom(s)`` = the same argmax without that condition: the most
    probable irrep, 0 only when every weight is below :c:data:`pg_thresh`.
    It is for the printed report, which names the dominant irrep together
    with its purity instead of saying "none"; nothing may constrain a state
    on it.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_n_col`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_proj_w`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: int pg_state_irrep


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        integer, allocatable	:: pg_state_dom	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the MO projector
    in the state, not an occupation-number guess. Computed by
    :c:func:`pg_state_weight_matrix`, which the irrep filter also calls.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess must not drive the filter or the zeroing.
    ``pg_state_dom(s)`` = the same argmax without that condition: the most
    probable irrep, 0 only when every weight is below :c:data:`pg_thresh`.
    It is for the printed report, which names the dominant irrep together
    with its purity instead of saying "none"; nothing may constrain a state
    on it.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_n_col`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_proj_w`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: double pg_state_irrep_weight


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        integer, allocatable	:: pg_state_dom	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the MO projector
    in the state, not an occupation-number guess. Computed by
    :c:func:`pg_state_weight_matrix`, which the irrep filter also calls.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess must not drive the filter or the zeroing.
    ``pg_state_dom(s)`` = the same argmax without that condition: the most
    probable irrep, 0 only when every weight is below :c:data:`pg_thresh`.
    It is for the printed report, which names the dominant irrep together
    with its purity instead of saying "none"; nothing may constrain a state
    on it.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_n_col`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_proj_w`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: double pg_state_norm


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        integer, allocatable	:: pg_state_dom	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the MO projector
    in the state, not an occupation-number guess. Computed by
    :c:func:`pg_state_weight_matrix`, which the irrep filter also calls.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess must not drive the filter or the zeroing.
    ``pg_state_dom(s)`` = the same argmax without that condition: the most
    probable irrep, 0 only when every weight is below :c:data:`pg_thresh`.
    It is for the printed report, which names the dominant irrep together
    with its purity instead of saying "none"; nothing may constrain a state
    on it.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_n_col`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_proj_w`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: double pg_state_purity


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_state_irrep_weight	(pg_n_irrep,N_states)
        integer, allocatable	:: pg_state_irrep	(N_states)
        integer, allocatable	:: pg_state_dom	(N_states)
        double precision, allocatable	:: pg_state_purity	(N_states)
        double precision, allocatable	:: pg_det_mixed_frac	(N_states)
        double precision, allocatable	:: pg_state_norm	(N_states)


    ``pg_state_irrep_weight(a,s)`` = w_a(s), the weight of irrep ``a`` in
    state ``s``: the diagonal matrix element of the MO projector
    in the state, not an occupation-number guess. Computed by
    :c:func:`pg_state_weight_matrix`, which the irrep filter also calls.
    
    ``pg_state_irrep(s)`` = argmax_a w_a(s), or 0 when the state has no
    definite irrep (spin contamination, an incomplete family, a broken
    MO basis): a best guess must not drive the filter or the zeroing.
    ``pg_state_dom(s)`` = the same argmax without that condition: the most
    probable irrep, 0 only when every weight is below :c:data:`pg_thresh`.
    It is for the printed report, which names the dominant irrep together
    with its purity instead of saying "none"; nothing may constrain a state
    on it.
    ``pg_state_purity(s)`` = max_a w_a(s).
    
    ``pg_det_mixed_frac(s)`` = sum of c_D(s)^2 over the determinants whose
    :c:data:`pg_det_irrep` is 0: the weight of the state carried by determinants
    that belong to no symmetry group on their own.
    
    ``pg_state_norm(s)`` = sum_D c_D(s)^2, the exact value of
    sum_a w_a(s). A deficit means psi_coef is not
    normalized in the determinant basis; it never means that the CI space
    is not group invariant.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_n_col`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_proj_w`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`


 
.. c:var:: double pg_xg


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_xg	(pg_n_elem,mo_num,mo_num)


    ``pg_xg(g,:,:)`` = C^T S M_g C, the matrix of group element ``g`` in the
    MO basis, built once. The element index is the leading one, so the
    ``mo_num x mo_num`` matrix of one element is STRIDED: element
    ``(g,i,j)`` sits at 1-based offset ``g + (i-1)*pg_n_elem +
    (j-1)*pg_n_elem*mo_num`` of the flat array. The determinant kernels
    therefore address it flat, through the ``xrow``/``xcol`` tables of
    :c:func:`pg_fam_setup`; :c:func:`pg_elem_me` keeps the readable
    ``pg_xg(g,i,j)`` form it always had.
    
    One scratch matrix at a time was the first plan, but the
    determinant projectors need every element for every family, and a
    rebuild is a full AO transform (three matmuls over ao_num). :math:`|G|` is at
    most 120 and the MO basis of a plugin-sized system a few hundred
    wide, so the set costs :math:`8|G|N_{MO}^2` bytes: 0.8 MB for benzene in
    6-31G. The report prints the number. A production version can keep
    only the rows and columns of the partially occupied blocks; nothing
    below cares, because the only readers are :c:func:`pg_elem_me` (the
    reference oracle), :c:data:`pg_detfull` and :c:func:`pg_family_vg`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`pg_n_elem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_state_irrep_weight`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void pg_binom(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer function pg_binom(dim, n)


    C(dim,n), 0 outside 0<=n<=dim. The product is built so that every
    division is exact (r*(dim-n+k)/k is an integer at each step), which
    is what keeps C(5,2) = 10 and C(12,6) = 924 worth of exactness in a
    four-byte integer.

 
.. c:function:: void pg_det_counts(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_det_counts(idet, nsel)


    ``nsel(s,b)`` = number of spin-``s`` electrons of determinant ``idet`` in
    symmetry block ``b``. Two determinants are in the same family if and
    only if their ``nsel`` agree.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`pg_col_diag`
       * :c:data:`psi_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_gates`
       * :c:data:`pg_det_irrep`
       * :c:func:`pg_det_weights`
       * :c:func:`pg_state_weight_matrix`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`

 
.. c:function:: void pg_det_gates(void)


    File : :file:`symmetry_det.irp.f`

    .. code:: fortran

        logical function pg_det_gates()


    ``D1``
        sum_a T^a = 1 on every family (validates chi, mult_a and the
        family structure at once).
    
    ``D2``
        (T^a)^2 = T^a and T^a Hermitian, every family, every a.
    
    D2b tests the production element kernel (:c:func:`pg_family_vg`, which
    computes the ordered pairs), while the diagonal used for labels and
    purity is independently cross-checked against :c:func:`pg_det_weights`
    -- which still runs on :c:func:`pg_elem_me` -- by the report's
    self-check.
    
    ``D3``
        -1e-8 <= purity <= 1+1e-8 for every determinant.
    
    ``D4``
        ``sum_a w_a(s) - sum_D c_D(s)^2`` <= 1e-8: the accumulation loop
        reproduces the sum rule of the state weights.
    
    ``D5``
        every labelled state has w_max >= 1-1e-6, and the mixed
        determinant fraction is reported.
    
    ``D6``
        the dominant determinant of state 1 has irrep 1 for a closed-shell
        Hartree-Fock wave function (the cross-check that the determinant
        and MO labellings agree). Gated only when the
        determinant's family has one member: a block that is partially
        filled makes even a closed-shell occupation a mixture, so such a
        configuration has no one-determinant irrep to compare with. That
        case is printed, and the state weights carry the answer.
    
    Families of size 1 are skipped by D1/D2: there T^a = delta is exact
    by construction, and the label already comes from that diagonal.
    
    ``D7``
        the tensor contraction of :c:func:`pg_family_wadd` against the
        reference double sum ``sum_{ip,jp} c(ip,s) V_g(ip,jp) c(jp,s)``
        built from :c:func:`pg_family_vg`, on a deterministic sample of
        families. This is the gate that makes the W3 rewrite honest: the
        contraction is a reassociation of that sum, so it may differ in
        the last bits but nowhere else. At most ``d7_nfam`` families with
        ``np > 1`` are sampled (the first ones, so the sample is a fixed
        function of the wave function), and each family is checked on
        ``max(2, min(|G|, d7_pairs/np**2))`` evenly strided elements --
        g = 1 and at least one non-identity element are always included.
        The relative measure is ``max|q_contract - q_ref| / max|q_ref|``
        over the sample, gated at 1e-13.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`pg_xg`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_counts`
       * :c:func:`pg_fam_setup`
       * :c:func:`pg_family_fill`
       * :c:func:`pg_family_index`
       * :c:func:`pg_family_projector`
       * :c:func:`pg_family_vg`
       * :c:func:`pg_family_wadd`

 
.. c:function:: void pg_det_report(void)


    File : :file:`symmetry_det.irp.f`

    .. code:: fortran

        subroutine pg_det_report


    Prints the labels and the gates, and stops with status 1 if any gate
    fails. Everything measured is printed before the verdict: one run
    shows the whole story.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_group_name`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_mo_mixed`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`read_wf`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`symmetry_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list_ab`
       * :c:func:`pg_det_weights`
       * :c:func:`pg_mo_top2`

 
.. c:function:: void pg_det_state_labels(void)


    File : :file:`state_labels.irp.f`

    .. code:: fortran

        subroutine pg_det_state_labels(ndet, coef, lda, nst, label, purity, mixed)


    Fast irrep labels of ``nst`` CI vectors stored in ``coef``.
    
    For each column k: ``purity(k)`` = wmax/tot, ``mixed(k)`` = wmix/tot and
    ``label(k)`` = argmax when wmax >= (1 - pg_thresh)*tot, else 0, where
    w(a) sums coef(i,k)**2 over the determinants of irrep a, and wmix sums
    the determinants whose :c:data:`pg_det_irrep` is 0 (mixed families).
    
    Everything is a ratio: the callers hand in vectors that are not
    normalized (Davidson's in-loop ``u_in`` is W*y), so an absolute test
    would label nothing. A zero column gives label 0, purity 0, mixed 0.
    
    The output arrays are intent(inout), not intent(out): they are
    assumed-size dummies, and an intent(out) dummy would leave the whole
    actual argument undefined whatever bounds the caller used. The per-k
    zeroing below defines them.
    
    Fast path, not the exact formula: a state whose weight sits on
    multi-member determinant families gets label 0 and is not constrained.
    :c:data:`pg_state_irrep_weight` (stage 3) stays the exact value, off-diagonal
    family terms included, and is what the report prints.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_thresh`

 
.. c:function:: void pg_det_weights(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_det_weights(idet, w)


    ``w(a)`` = T^a_{DD} = <D|P^a|D> for determinant ``idet``: the probability
    that D carries irrep ``a``. Sums to 1 over ``a`` (gate D1's identity on
    the diagonal), so a determinant with no definite irrep is described by
    a distribution, not by a label.
    
    This is the diagonal of :c:func:`pg_family_projector` at one member, not
    the whole matrix: |G|*|F| element matrix elements instead of |G|*|F|^2,
    which is what makes it affordable to call once per printed row of the
    report. The coefficients are the ones of the projector; the caller can
    check ``max_a w(a)`` against :c:data:`pg_det_purity`, computed from the full
    matrix, and the driver does exactly that.
    
    ``w`` is left at 0 when the determinant is outside 1:N_det, when the
    labelling failed, or when its family exceeds :c:data:`pg_family_cap`; the
    caller decides what to print. Call it only after :c:data:`pg_det_irrep` is
    built: it re-enumerates the family and does not relabel anything.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_family_cap`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_report`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_counts`
       * :c:func:`pg_elem_me`
       * :c:func:`pg_family_fill`
       * :c:func:`pg_family_index`
       * :c:func:`pg_family_size`

 
.. c:function:: void pg_elem_me(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_elem_me(ig, np, mask, ip, jp, val)


    ``val`` = :math:`<D_ip|U(g_ig)|D_jp>` = prod_b det( X_g[A'_b,A_b] ) * prod_b det( X_g[B'_b,B_b] ).
    
    Every block contributes a square submatrix because the family has
    the same count in every block. The counts are computed first and the
    submatrix filled afterwards: the two must not share a running index,
    or a block with a zero count inherits the previous block's value.
    
    An over-full or inconsistent pair returns exactly 0 rather than a
    number, and an exact zero is propagated exactly, so the caller may
    skip the accumulation of that pair.
    
    Kept as the independent reference used by :c:func:`pg_det_weights`;
    no hot path calls it (the labels use :c:func:`pg_family_labels`, the
    state weights :c:func:`pg_family_wadd`).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_xg`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_weights`

 
.. c:function:: void pg_fam_amat_size(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer function pg_fam_amat_size(npart, mrad)


    Number of entries of the flat sector-matrix buffer of
    :c:func:`pg_family_vg`: ``sum_p mrad(p)**2``. The radices multiply to
    the family size ``np``, so this is at most ``np**2``: the buffer never
    needs more memory than the ``v(np,np)`` matrix it feeds. The caller
    allocates it once per family, never per group element.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_col_diag`

 
.. c:function:: void pg_fam_cfac(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        double precision function pg_fam_cfac(g, nfull, fsec)


    Product of the block determinants of the full sectors of element ``g``:
    the factor common to every pair of the family. 1.d0 when no sector is
    full, exactly 0.d0 as soon as one factor is 0 (the whole family is
    annihilated by that element, and an exact zero is propagated exactly).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_col_diag`
       * :c:data:`pg_detfull`

 
.. c:function:: void pg_fam_setup(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_fam_setup(np, nsel, mask, npart, nfull, fsec, psec, pn, pos, xrow, xcol, &
                                mrad, strd, rep, aoff, si)


    Classifies the sectors of one family once, independently of the group
    element and of the determinant pair.
    
    A sector is one (spin, block) pair. The family has the same count in
    every block, so a sector is either
    
    * empty: nothing to multiply;
    * full: all its MOs are occupied, its block determinant is the same
      for every member and every pair, and it is folded into
      :c:data:`pg_detfull` (``fsec``);
    * partial: the only sectors the pair kernel has to look at (``psec``,
      with ``pn(p)`` the number of occupied MOs).
    
    ``pos(k,ip,p)`` is the MO of the ``k``-th occupied orbital of the
    partial sector ``p`` for family member ``ip``, in ascending MO order:
    the precomputed occupied-MO table that replaces the per-(g,pair,sector)
    bit scans of :c:func:`pg_elem_me`. Ascending order is what makes the
    submatrix ordering identical to the reference kernel.
    
    ``xrow(k,ip,p)`` and ``xcol(k,ip,p)`` are that same MO written as the
    two flat offsets of :c:data:`pg_xg`, which is
    ``(pg_n_elem,mo_num,mo_num)``: element ``(g,i,j)`` is at 1-based
    offset ``g + (i-1)*pg_n_elem + (j-1)*pg_n_elem*mo_num``. Both products
    depend only on the MO, so hoisting them here turns every access of
    :c:func:`pg_family_vg` and :c:func:`pg_fam_diag` into two integer adds
    and two loads instead of two multiplications.
    
    The family enumeration (:c:func:`pg_family_fill`) is an odometer over
    the partial sectors: member ``ip`` carries one digit per sector, and
    the sector submatrix depends only on that digit. This routine also
    returns that mixed-radix description:
    
    * ``mrad(p)`` = C(d_b, pn(p)), the number of occupations of the
      partial sector ``p``;
    * ``strd(p)`` = the stride of the sector in the odometer of
      :c:func:`pg_family_fill`, i.e. the product of the radices of the
      (spin, block) slots that precede it in that enumeration;
    * ``si(ip,p)`` = ``mod((ip-1)/strd(p), mrad(p)) + 1``, the digit of
      member ``ip`` in sector ``p``, computed once here and never in a
      loop over pairs;
    * ``rep(p,i)`` = the first member whose digit in sector ``p`` is
      ``i``. Every member with that digit occupies the same MOs of the
      sector (the digit is the rank of the subset, and the unranking is
      injective), so ``xrow``/``xcol`` read at ``rep(p,i)`` are the
      subset of the digit: the sector matrix of :c:func:`pg_family_vg`
      reuses the occupied-MO tables without a second unranking;
    * ``aoff(p)`` = the 1-based offset of entry ``(1,1)`` of sector ``p``
      in the flat sector-matrix buffer of :c:func:`pg_family_vg`, whose
      length is :c:func:`pg_fam_amat_size`.
    
    ``psec`` is ordered ``s`` outer / ``b`` inner while the odometer is
    ``b`` outer / ``s`` inner. ``psec`` is NOT reordered (the product over
    sectors in :c:func:`pg_family_vg` follows it, and the multiplication
    order is part of the result); the strides are looked up in the
    odometer order. ``prod_p mrad(p) == np`` is asserted here, the same
    invariant as the ``mtot /= np`` test of :c:func:`pg_family_fill`.
    
    The tables are transient, one family at a time: ``12*np*2*nblock``
    integers each at the deepest call, plus ``np*2*nblock`` integers for
    each of ``si`` and ``rep``. ``npart = -1`` means the family was
    refused because a block holds more than :c:data:`pg_det_nmax`
    occupied MOs; :c:data:`pg_label_ierr` then carries the reason and no
    caller reads ``pos`` or the offsets. The flat offsets require
    ``pg_n_elem*mo_num**2`` to fit in a default integer, which is implied
    by the size of :c:data:`pg_xg` itself.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_elem`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_gates`
       * :c:func:`pg_family_labels`
       * :c:func:`pg_family_projector`
       * :c:func:`pg_family_wadd`

 
.. c:function:: void pg_family_diag_acc(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_diag_acc(np, nfull, fsec, npart, psec, pn, xrow, xcol, xgf, acc)


    Adds ``pg_proj_w(a,g) * <D_ip|U(g)|D_ip>`` to ``acc(ip,a)`` for every
    element ``g`` and every member ``ip`` of the family.
    
    This is the diagonal-only kernel of :c:func:`pg_family_labels`. It
    computes exactly the per-pair product that :c:func:`pg_family_vg`
    obtained before the sector matrices were factored out, with
    ``jp = ip``: the same 1-term / 2-term / ``pg_small_det`` expressions
    (which now also build the entries of every ``A_p``), the same partial
    product starting at 1.d0, the same exact-zero short circuit and the
    same ``acc`` accumulation. The EXPRESSIONS and the multiplication
    order MUST stay identical; that identity is what makes the labels
    bit-identical to the pre-W2 reference. What changes is only the cost:
    one call per family instead of ``|G|*np`` calls, and the offsets stay
    in registers.
    
    ``g`` is the outer loop and ``ip`` inner, i.e. the order in which the
    labels have always been accumulated, so each ``acc(ip,a)`` receives
    the same summands in the same order as before.
    
    No ``pos``: the offsets carry all the addressing, and the occupied-MO
    table is not read here.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_col_diag`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_proj_w`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_family_labels`

 
.. c:function:: void pg_family_fill(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_fill(nsel, np, mask)


    ``mask(s,b,ip)`` is the bit mask of the MOs of block ``b`` occupied with
    spin ``s`` by family member ``ip``, bit ``k`` standing for the
    MO ``pg_block_mo(pg_block_start(b)+k)``.
    
    An odometer over the per-(spin,block) combination ranks: ``scalar``
    index q in [0,np-1] decomposed in mixed radix, each digit unranked to
    a subset. Cost np * sum_b d_b, paid once per family.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_label_ierr`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_gates`
       * :c:data:`pg_det_irrep`
       * :c:func:`pg_det_weights`
       * :c:func:`pg_state_weight_matrix`

 
.. c:function:: void pg_family_index(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_index(np, mask, idx, imin)


    ``idx(ip)`` is the index in :c:data:`psi_det` of family member ``ip``,
    or 0 when that member is not in the wave function (a truncated CI
    space: the missing member simply has coefficient 0).
    ``imin`` is the smallest nonzero index.
    
    The lookup is QP2's own binary search on :c:data:`psi_det_sorted_bit`, which
    is why the members are built as bit strings through :c:func:`list_to_bitstring`
    and not compared determinantal by determinantal.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`psi_det_sorted_bit`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_gates`
       * :c:data:`pg_det_irrep`
       * :c:func:`pg_det_weights`
       * :c:func:`pg_state_weight_matrix`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`list_to_bitstring`

 
.. c:function:: void pg_family_labels(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_labels(np, nsel, mask, lab, pur)


    ``pur(ip)`` = max_a T^a_{ip,ip} and ``lab(ip)`` = the argmax when it
    reaches 1 - :option:`symmetry pg_thresh`, else 0.
    
    :math:`T^a_{DD}` needs only the DIAGONAL of the family: ``|G|*np``
    matrix elements instead of ``|G|*np^2``. ``pur`` is what
    :c:data:`pg_det_purity` was, and ``lab`` follows the same
    ``1 - pg_thresh`` rule as :c:data:`pg_det_irrep`.
    
    The accumulation is per (member, irrep): the irrep loop is inside the
    element loop but outside the pair loop, and the (irrep, element)
    weight is read from :c:data:`pg_proj_w`. It is done by
    :c:func:`pg_family_diag_acc`, the diagonal kernel: one call per
    family instead of one :c:func:`pg_family_vg` call per (element,
    member), with the same sector order and the same accumulation order.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_col_diag`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_thresh`
       * :c:data:`pg_xg`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_fam_setup`
       * :c:func:`pg_family_diag_acc`

 
.. c:function:: void pg_family_projector(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_projector(np, nsel, mask, t)


    ``t(ip,jp,a)`` = T^a_{ip,jp} = (d_a/(mult_a*:math:`|G|`)) sum_g chi_a(col(g)) :math:`<D_ip|U(g)|D_jp>`.
    
    This is the determinant version of the MO projector: the same
    coefficients, the same characters, the same ``mult_a`` as
    :c:data:`pg_mo_weight`. T^a is a Hermitian
    projector on the family (gates D1/D2), and its diagonal gives the
    determinant labels.
    
    The REFERENCE full matrix, kept for gates D1/D2 only: the labels read
    the diagonal through :c:func:`pg_family_labels` and the state weights
    the scalar :math:`q(g,s)` through :c:func:`pg_family_wadd`. Ordered
    pairs are summed, so gate D2b (symmetry of :math:`T^a`) stays
    meaningful. Skipping an exactly zero character is exact: it would add
    a product whose first factor is 0.
    
    The mixed-radix tables and the sector-matrix buffer of
    :c:func:`pg_family_vg` are built once here, before the ``do g`` loop,
    exactly as in :c:func:`pg_family_wadd`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_col_diag`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_proj_w`
       * :c:data:`pg_xg`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_gates`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_fam_setup`
       * :c:func:`pg_family_vg`

 
.. c:function:: void pg_family_size(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_size(nsel, np)


    ``np`` = prod_{s,b} C(d_b, nsel(s,b)), the number of determinants with
    the same per-block spin counts. A block that is empty or full
    contributes one choice. The product is cut before any multiplication
    that could overflow, so the cap test of the caller is reliable, and
    ``np`` is left negative when a block is wider than
    :c:data:`pg_block_cap` -- the caller reports that case with its own
    message.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_block_cap`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_family_cap`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:func:`pg_det_weights`

 
.. c:function:: void pg_family_vg(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_vg(g, np, nfull, fsec, npart, psec, pn, pos, xrow, xcol, &
                                mrad, rep, aoff, si, amat, xgf, v)


    ``v(ip,jp)`` = :math:`<D_{ip}|U(g)|D_{jp}>` on the family, i.e. the
    N-electron matrix of ONE group element: the factor of the full sectors
    times the partial-sector product.
    
    This is the reference element kernel, no longer on the production
    path: :c:func:`pg_family_projector` sums it over the elements for
    gates D1/D2, and gate D7 compares it with the tensor contraction of
    :c:func:`pg_family_wadd`. All ordered pairs are computed by both, and
    no caller assumes ``v`` symmetric, which is what keeps gate D2b
    meaningful.
    
    FACTORISATION. The family enumeration is an odometer, so member ``ip``
    carries one digit ``si(ip,p)`` per partial sector and the sector
    submatrix depends only on the two digits:
    
      :math:`V(ip,jp) = cfac * prod_p A_p(si(ip,p), si(jp,p))`,
      :math:`A_p(i,j) = det X_g[subset_i, subset_j]`.
    
    ``A_p`` is built once per (family, g) into the caller-owned ``amat``:
    ``sum_p mrad(p)**2`` sector determinants instead of ``np**2`` for
    every element, and the pair loop is then one load and one multiply per
    (pair, sector). The four tables come from :c:func:`pg_fam_setup`, so
    the digit of a member is a load, not a division.
    
    BIT-IDENTITY with the per-pair kernel of W1. Every ``A_p(i,j)`` is
    computed by the same expression, from the same ``xgf`` entries, as the
    pair ``(ip,jp)`` carrying those digits: ``rep(p,i)`` is a member that
    occupies the same MOs of the sector, so ``xrow``/``xcol`` read at it
    are the same integers. ``val`` starts at 1.d0, the sectors multiply it
    in ascending ``p`` (``psec`` order, which this routine must not
    reorder), an exact zero exits the sector loop, and ``cfac`` is applied
    LAST, so ``v = cfac*(((1*d1)*d2)*...)`` exactly as before. All sector
    entries are built even when another sector is exactly zero: the pair
    product still exits at the same sector and never reads the unused
    entries, so no value moves. ``xgf`` holds finite representation
    matrices, so building an entry W1 would have skipped reads a finite
    number and computes nothing else.
    
    ``xgf`` is :c:data:`pg_xg` seen as one flat array (a sequence
    association on ``pg_xg(1,1,1)``), and ``xrow``/``xcol`` are the flat
    offsets of the two MOs, precomputed per family by
    :c:func:`pg_fam_setup`. Each access is then ``xgf(g + row + col)``:
    the two multiplications of the address computation are gone.
    
    ``amat`` is scratch: its first :c:func:`pg_fam_amat_size` entries are
    overwritten for the current ``g``. It is caller-owned and allocated
    once per family, next to ``v``, and never inside the ``do g`` loop:
    ``np**2`` doubles is 1.28 MB at ``np = 400``, which must not be an
    automatic array inside an OpenMP region. ``pos`` is kept in the
    interface (and is still the documented occupied-MO table) although the
    offsets carry all the addressing information.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_col_diag`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_gates`
       * :c:func:`pg_family_projector`

 
.. c:function:: void pg_family_wadd(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        subroutine pg_family_wadd(np, nsel, mask, c, ldc, nst, xgf, gsel, qtot)


    ``qtot(g,s) =! ``qtot(g,s) +( c^T V_g c`` for the ``nst`` CI columns ``c`` on this)
    family. The irrep weight is NOT applied here: the caller sums
    ``qtot`` over the families and applies :c:data:`pg_proj_w` once per
    element,
    
      sum_F sum_g pmw(a,g) q_F(g,s) = sum_g pmw(a,g) (sum_F q_F(g,s)),
    
    so the number of additions to the weight vector drops from
    ``nst*n_irrep*n_family*|G|`` to ``nst*n_irrep*|G|``. The epilogue
    lives in :c:func:`pg_state_weight_matrix`.
    
    CONTRACTION, not the element matrix. With
    :math:`V_g(ip,jp) = cfac * prod_p A_p(si(ip,p),si(jp,p))`, the ket
    digit of each partial sector is contracted in turn:
    
      u <- c ;  for p:  u'(i,rest,s) = sum_q A_p(i,q) u(q,rest,s)
    
    and :math:`q(s) = cfac * sum_ip c(ip,s) u(ip,s)`. The odometer is
    mixed radix, so the members that share the digits of every sector
    other than ``p`` form the block
    ``{j0 + l + (k-1)*strd(p), k=1..mrad(p)}``, enumerated by ``h``
    (higher digits) and ``l`` (lower digits). Cost per (family, g):
    ``nst*sum_p np*mrad(p)`` multiplies plus ``sum_p mrad(p)**2`` sector
    determinants and ``nst*np`` for the final dot, instead of
    ``nst*np**2*npart`` pair products on the full ``np**2`` matrix.
    
    REASSOCIATIONS relative to the ``c^T V c`` double sum, all named here
    because this routine is allowed to differ from W1/W2 in the last
    bits and nowhere else:
    
    * the ``jp`` sum is performed mode by mode instead of per pair. The
      exact sum is the same; rounding differs by O(npart*max_p mrad(p)*eps)
      relative, i.e. at most ~1e-12 for the largest sector the family
      enumeration supports (C(12,6) = 924), six orders of magnitude below
      the ``1 - pg_thresh`` = 1-1e-6 label rule that consumes these
      weights;
    * ``cfac`` is applied once to ``q(s)`` instead of once per pair
      (``c*(cfac*v)*c`` -> ``cfac*(c*u)``): a rescaling by a common
      factor, same argument, one rounding per state instead of np**2;
    * an exactly zero ``A_p`` zeroes ``u`` exactly, so an element that
      annihilates the family still contributes exactly 0; the old
      ``val == 0`` short circuit is not needed and would not change a
      value.
    
    ``c`` is indexed by family member, not by determinant: the caller
    gathers it, and a missing member is a zero column contribution.
    ``c(ldc,nst)`` holds the family in its first ``np`` rows, ``ldc >= np``:
    the caller keeps one buffer per thread and reuses it for families of
    different size, which an exact ``c(np,nst)`` shape would forbid (the
    leading extent is the one that enters the address).
    
    ``xgf`` is :c:data:`pg_xg` seen as one flat array (a sequence
    association on ``pg_xg(1,1,1)``, as in :c:func:`pg_family_vg`), so the
    sector build reads ``xgf(g + row + col)`` and pays no multiplication
    for the address.
    
    The sector-matrix buffer is built ONCE per family, before the
    ``do g`` loop, and overwritten for every element; it is
    ``sum_p mrad(p)**2`` doubles (:c:func:`pg_fam_amat_size`), i.e. never
    more than ``np**2``. ``u``, ``un`` and ``amat`` are allocated in one
    statement per family, never inside the ``do g`` loop. ``si`` and
    ``pos`` are still returned by :c:func:`pg_fam_setup` but the
    contraction addresses members by the odometer block structure and
    does not read them.
    
    :c:func:`pg_family_vg` is no longer called here: it stays as the
    independent reference for gates D2b and D7.
    
    ``gsel(g)`` selects the group elements to contract. It is a
    gate/instrumentation hook, not a physics switch: the production caller
    passes an all-true mask, so the result is unchanged, and the only other
    caller is D7, which selects the elements it can afford to check. An
    unselected element is skipped by the first statement of the ``do g``
    loop, before :c:func:`pg_fam_cfac`, so it costs one logical test and
    nothing else; ``qtot`` keeps its incoming value for that element, which
    is exactly what a caller that never reads it expects.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_col_diag`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_elem`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_gates`
       * :c:func:`pg_state_weight_matrix`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`move_alloc`
       * :c:func:`pg_fam_setup`

 
.. c:function:: void pg_labels_usable(void)


    File : :file:`state_labels.irp.f`

    .. code:: fortran

        logical function pg_labels_usable()


    Cheap guard for the stage-4 symmetry work: true when there is a
    determinant, the point group was detected and the MOs carry definite
    irreps. Never stops and never builds the determinant labelling: it
    names ``N_det`` and the two stage-2 functions only, so a caller can test
    it before paying for anything.
    
    It has no printing switch to test: the state table is always printed
    when the labels are usable, and a printing decision must not disable
    the irrep filter. :c:data:`auto_symmetry` stays the
    master switch, through :c:func:`pg_symmetry_ready`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`

 
.. c:function:: void pg_print_state_irreps(void)


    File : :file:`state_report.irp.f`

    .. code:: fortran

        subroutine pg_print_state_irreps(energy, s2)


    Prints the irrep of each CI state. Never stops. It is silent when the
    detection was not requested or failed: :c:func:`pg_symmetry_ready` returns
    .false. and the failure has already printed its reason. A state with no
    definite irrep is described by its irrep probabilities, which are exact
    for the stored CI space even when the MO basis is not symmetry adapted.
    ``energy`` and ``s2`` are the caller's arrays, at
    least ``N_states`` long.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_irreps_label`

 
.. c:function:: void pg_print_state_irreps_work(void)


    File : :file:`state_report.irp.f`

    .. code:: fortran

        subroutine pg_print_state_irreps_work(energy, s2)


    The table. :c:data:`pg_state_irrep_weight` delivers the weights, the labels,
    the purities and the norm in one provider block.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`pg_group_name`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_irrep_filter`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_mo_mixed`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`pg_thresh`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_irreps_label`

 
.. c:function:: void pg_small_det(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        double precision function pg_small_det(a, n)


    det(a(1:n,1:n)), ``a`` untouched. ``n`` <= :c:data:`pg_det_nmax` (a
    block never has more than d_a <= 5 occupied MOs; the guard is in
    :c:func:`pg_fam_setup`).
    
    The function returns exactly 0 whenever the computed determinant is
    exactly 0. The ``n = 1`` and ``n = 2`` branches are explicit products
    (a 2x2 product detects the cancellation of ``ad - bc`` that the
    pivoted elimination does not see), the general case is Gaussian
    elimination with partial pivoting. The two forms can differ by one
    ulp, never more.

 
.. c:function:: void pg_state_apply(void)


    File : :file:`state_labels.irp.f`

    .. code:: fortran

        subroutine pg_state_apply(label, n, ok)


    Narrow ``ok(1:n)`` to the states whose known irrep is in :c:data:`irrep_states`.
    
    ``label(k)`` = 0 means "unknown": such a state is kept, because a filter
    must not drop a state it may not be subject to.
    
    Called once, after the diagonalization, on the states that were actually
    computed: the filter reorders them, it does not steer the Davidson
    iteration (a swap inside the loop made the iteration wander and, for some
    fixtures, fail to converge).
    
    :c:data:`pg_irrep_allowed` is provided here, not in the callers: IRPF90 hoists
    the explicit PROVIDE to the entry of this routine, so the mask and the
    detection chain are built at the first filtered call and by no
    caller's entry block.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_irrep_allowed`

 
.. c:function:: void pg_state_irrep_head(void)


    File : :file:`cipsi_report.irp.f`

    .. code:: fortran

        subroutine pg_state_irrep_head(label, dom, purity, head)


    ``"   Irrep: <label>"`` for one state, used by the extrapolation
    headings. A state with no definite label shows its most probable irrep
    and the weight fraction behind it, instead of ``none``; it falls back to
    ``none`` only when even the dominant weight is negligible. Called only
    after :c:func:`pg_state_irreps_psi` returned true.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_col`

 
.. c:function:: void pg_state_irreps_label(void)


    File : :file:`state_report.irp.f`

    .. code:: fortran

        subroutine pg_state_irreps_label(energy, s2)


    Forces the determinant labelling and gates the report on
    :c:data:`pg_label_ierr`. Separate from the dispatcher because IRPF90 hoists
    both the implicit needs and an explicit PROVIDE to the procedure
    entry: a body holding the PROVIDE would label before the guards.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_print_state_irreps`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_print_state_irreps_work`

 
.. c:function:: void pg_state_irreps_psi(void)


    File : :file:`cipsi_report.irp.f`

    .. code:: fortran

        logical function pg_state_irreps_psi(nst, label, dom, purity)


    Irrep labels of the states stored in :c:data:`psi_coef`. ``label(k)``
    is the definite label of state ``k`` (0 when the state has no definite
    irrep), ``dom(k)`` the most probable one and ``purity(k)`` its weight
    fraction -- what the report prints in place of ``none``. ``nst`` entries.
    Returns .false. -- and writes nothing -- when the symmetry is not usable,
    so a caller can skip the whole row instead of printing blanks.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_irreps_psi_work`

 
.. c:function:: void pg_state_irreps_psi_work(void)


    File : :file:`cipsi_report.irp.f`

    .. code:: fortran

        subroutine pg_state_irreps_psi_work(nst, label, dom, purity)


    Entered only after :c:func:`pg_state_irreps_psi` returned true.
    :c:data:`pg_state_irrep` is the exact weight-matrix label of
    :c:data:`psi_coef`, the value the end-of-run table prints;
    :c:data:`pg_state_dom` is the same argmax without the purity test.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`pg_state_irrep_weight`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_irreps_psi`

 
.. c:function:: void pg_state_irreps_vec(void)


    File : :file:`cipsi_report.irp.f`

    .. code:: fortran

        logical function pg_state_irreps_vec(nst, coef, lda, label, dom, purity)


    Irrep labels of ``nst`` CI vectors stored in ``coef(lda,*)``, for the
    vectors that were just diagonalized. ``label(k)`` is 0 when the state
    has no definite irrep, ``dom(k)`` the most probable irrep (0 when there
    is none) and ``purity(k)`` its dominant weight fraction. Exact weight
    matrix, the same call the irrep filter makes: a state cannot be selected
    on one definition of its symmetry and printed on another.
    
    Returns .false. -- and writes nothing -- when the symmetry is not
    usable.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_irreps_vec_work`

 
.. c:function:: void pg_state_irreps_vec_work(void)


    File : :file:`cipsi_report.irp.f`

    .. code:: fortran

        subroutine pg_state_irreps_vec_work(nst, coef, lda, label, dom, purity)


    Entered only after :c:func:`pg_state_irreps_vec` returned true.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_irreps_vec`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_labels_exact`

 
.. c:function:: void pg_state_labels_exact(void)


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        subroutine pg_state_labels_exact(ndet, coef, lda, nst, label, dom, purity)


    Exact irrep labels of ``nst`` CI vectors, for the irrep filter.
    
    ``label(k)`` = argmax_a w_a(k) when that maximum is at least
    1 - :c:data:`pg_thresh`, else 0 ("no definite irrep": the filter must
    not constrain such a state). Same rule as :c:data:`pg_state_irrep`,
    because it is the same weight matrix: the state that is selected and
    the state that is printed carry the same label.
    
    ``dom(k)`` = the same argmax without the purity condition: the most
    probable irrep, for the report only (0 when every weight is below
    :c:data:`pg_thresh`). ``purity(k)`` = max_a w_a(k)/norm(k).
    
    Not for the Davidson loop: this enumerates the families. The loop uses
    :c:func:`pg_det_state_labels`, which only sums coef^2 by determinant
    label and is deliberately allowed to answer 0.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_col`
       * :c:data:`pg_thresh`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_irreps_vec_work`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_state_weight_matrix`

 
.. c:function:: void pg_state_weight_matrix(void)


    File : :file:`state_weights.irp.f`

    .. code:: fortran

        subroutine pg_state_weight_matrix(ndet, coef, lda, nst, w, norm, dmix)


    Exact irrep weights of ``nst`` CI vectors stored in ``coef(lda,*)``.
    
    ``w(a,k)`` = weight of irrep ``a`` in column ``k``, ``norm(k)`` = sum_D
    coef_D^2 and ``dmix(k)`` = the part of ``norm(k)`` carried by
    determinants with :c:data:`pg_det_irrep` = 0.
    
    A determinant with :c:data:`pg_det_block` = 1 is pure and its family is
    itself, so T^a_{DD} = delta_{a,pg_det_irrep(D)} exactly: those
    determinants are accumulated from the label, and the projectors are
    built only for the families of size > 1 (which are exactly the
    determinants with no definite irrep). The cost is
    
      O( sum_F [ |G| * (nst*sum_p np*mrad(p) + sum_p mrad(p)^2 + nst*np) ]
         + nst*n_irrep*|G| )
    
    instead of O(n_irrep nst sum_F np^2) or O(|G| nst sum_F np^2 npart):
    the matrix of one group element is never formed, the coefficient
    tensor is contracted mode by mode (:c:func:`pg_family_wadd`), one
    scalar per (element, state) is kept, and the irrep projector is
    applied once over the families, because
    
      sum_F sum_g pmw(a,g) q_F(g,s) = sum_g pmw(a,g) (sum_F q_F(g,s)).
    
    This is not the O(n_irrep nst ndet) diagonal-only approximation of
    :c:func:`pg_det_state_labels`.
    
    The weights may differ from the previous implementation at the 1e-15
    relative level. Three reassociations, all below the ``1 - pg_thresh``
    rule (1-1e-6) by at least six orders of magnitude:
    
    * the inner ``jp`` sum of ``c^T V_g c`` is done mode by mode inside
      :c:func:`pg_family_wadd` (see there);
    * ``q_F(g,s)`` is summed over the families BEFORE the projector
      weight is applied, instead of applying ``pmw`` to every family and
      accumulating the products. The summands are the same, grouped by
      ``g``; the relative change is O(n_family*eps) of the group sum;
    * the pure-determinant weight stays the seed of ``w`` (below, before
      the family loop) and the whole projector contribution is added to
      it as ONE grouped sum at the end, where it used to be interleaved
      with the pure term family by family. Same addends, one grouping
      change.
    
    The printed precision and the label rule are unaffected.
    
    ``coef`` is assumed-size: the caller passes ``psi_coef`` or
    ``CI_eigenvectors``, both stored as (ndet, nst).
    
    The family loop is the part that parallelizes: each thread owns one
    contiguous block of :c:data:`pg_n_family` (``SCHEDULE(STATIC)``, so the
    partition is fixed by the family count and the thread count) and sums
    its families into its own ``qtot_c(:,:,ic)`` slice, in the serial
    family order. The slices are added afterwards in chunk order. With one
    chunk the group-element sums are accumulated in exactly the serial
    family order; with more chunks the same sum is reassociated by chunk,
    ``((q_1 + q_2) + q_3)`` instead of one running sum: last-bit
    differences only, reproducible for a fixed thread count, and eleven
    orders of magnitude below the ``1 - pg_thresh`` label rule.
    ``reduction(+:)`` is not used: its association is not specified, so
    two runs would not agree.
    
    ``qtot_c`` replaces the former weight slice ``wc``: it holds
    ``sum_F q_F(g,s)`` per element and state (``pg_n_elem*nst`` doubles
    per chunk, ~31 kB at |G| = 124, nst = 4, 8 chunks), not a weight. The
    pure-determinant seed stays in ``w`` and the irrep epilogue is applied
    to it once, after the merge.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_int`
       * :c:data:`pg_block_start`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_detfull`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_n_col`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_proj_w`
       * :c:data:`pg_xg`
       * :c:data:`psi_det`
       * :c:data:`psi_det_sorted_bit`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_state_irrep_weight`
       * :c:func:`pg_state_labels_exact`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_det_counts`
       * :c:func:`pg_family_fill`
       * :c:func:`pg_family_index`
       * :c:func:`pg_family_wadd`

 
.. c:function:: void pg_unrank(void)


    File : :file:`det_irrep.irp.f`

    .. code:: fortran

        integer function pg_unrank(dim, n, rank)


    Bit mask of the ``rank``-th ``n``-subset of {0,...,dim-1}. ``dim`` <= 30
    (a symmetry block is an irrep space and has d_a <= 5 elements; the
    guard is in :c:func:`pg_fam_setup`).

 
.. c:function:: void pg_write_state_irrep(void)


    File : :file:`state_report.irp.f`

    .. code:: fortran

        subroutine pg_write_state_irrep(iunit, istate, label, dom, purity)


    One line giving the irrep of state ``istate``, for the block printed
    by the :c:data:`CI_energy` provider after a diagonalization. Aligned on
    the ``A50`` label column of :c:func:`write_double`, so it sits under
    "Energy of state" and "S^2 of state".
    
    ``label`` = 0 means no definite irrep: the line then names the most
    probable irrep ``dom`` and gives ``purity``, the dominant weight
    fraction, rather than inventing a definite label. It says ``none`` only
    when every irrep weight is negligible. Never stops. Called only after
    the labels were computed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`pg_n_col`

 
.. c:function:: void pg_zero_state_irreps(void)


    File : :file:`zero_forbidden.irp.f`

    .. code:: fortran

        subroutine pg_zero_state_irreps()


    Zero the determinant coefficients that symmetry forbids.
    
    ``psi_coef(i,s)`` is set to zero when the determinant ``i`` and the
    state ``s`` both carry a definite irrep and the two differ. A label 0
    means "unknown" (a determinant in a degenerate family, a state mixing
    irreps): such a row or column is never touched, because a best guess
    would destroy amplitude the labels do not justify.
    
    Call it before the last save of the CI vector: the stored vector
    carries exact zeros instead of roundoff. Same projector and same
    threshold as :c:data:`pg_det_irrep` and :c:data:`pg_state_irrep`;
    no new |EZFIO| option.
    
    The guard lives in this dispatcher because IRPF90 hoists the providers
    named in a body to the procedure entry: a guard in the body would run
    after the labelling.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_zero_state_irreps_work`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`psi_coef`

 
.. c:function:: void pg_zero_state_irreps_work(void)


    File : :file:`zero_forbidden.irp.f`

    .. code:: fortran

        subroutine pg_zero_state_irreps_work()


    Worker of :c:func:`pg_zero_state_irreps`: names the providers, so it
    may only be called once the guard returned true.
    
    ``TOUCH psi_coef`` runs only when something was actually zeroed, so a
    wave function that is already clean does not force a re-sort.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_label_ierr`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`pg_state_irrep_weight`
       * :c:data:`psi_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_zero_state_irreps`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`psi_coef`

