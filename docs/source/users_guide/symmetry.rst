.. _symmetry:

==================
Molecular symmetry
==================

The :ref:`module_symmetry` module finds the point group of the molecule from
the nuclear geometry, labels the |MOs| with irreducible representations
(irreps), labels the determinants, and prints the irrep of every |CI| state.
This page describes the user workflow. The algorithms and the conventions are
on the module pages: :ref:`module_symmetry` (how the group is found),
:ref:`module_symmetry_mo` (|MO| labels), :ref:`module_symmetry_det`
(determinants, states, and the irrep filter).

Symmetry is a label, not a constraint
-------------------------------------

In this version the point group is used to *label* what the calculation
produces. It does not reduce the determinant space, it does not impose a block
structure on the Hamiltonian, and it does not symmetrize the geometry. The
:ref:`fci` run explores the same determinants as before; it merely says which
irrep each state belongs to, and it can keep a chosen subset of irreps in the
final diagonalization.

Two consequences are worth stating before looking at any output:

* the labels are as good as the geometry. A molecule that is symmetric only
  within :option:`symmetry pg_tol` is treated as symmetric, and every label
  inherits that geometric error;
* a molecular orbital that is not symmetry adapted gets **no** label. The
  programs then say so and continue without symmetry. They never stop and
  never guess.

What runs by itself
-------------------

At the end of :ref:`scf` and of the Kohn-Sham |SCF|, the point group is
detected and the orbitals are labelled. For water in 6-31G::

    symmetry: point group C2v, |G| = 4, geometry residual 1.8E-16 bohr
    symmetry:   A1           7 MOs
    symmetry:   A2           0 MOs
    symmetry:   B1           2 MOs
    symmetry:   B2           4 MOs
    symmetry:   4 irreps, 13 MOs, mo_basis/mo_symmetry written

The interval is the largest distance between a nucleus and its image under the
symmetry operations that were accepted, in bohr: it is the honest measure of
how symmetric the input geometry is. The labels are also stored in the |EZFIO|
database, in ``mo_basis/mo_symmetry``.

Nothing is printed and nothing is stored when :option:`symmetry auto_symmetry`
is |false|. If the detection fails, the reason is printed once and the rest of
the run is unaffected.

Checking the detection
----------------------

``qp_run symmetry_mo`` prints the complete report of the first two layers: the
canonical frame, the group, its classes, the character table with the
geometric meaning of every column, the |AO| basis check, every |MO| label and
six gates. It exits with status 1 when a gate fails, *after* printing every
number, so it can be chained in a test script.

Two numbers decide whether the labels can be trusted:

``delta_MO``
    how far the |MO| representation matrices are from being orthogonal,
    ``max_g |X_g^T X_g - 1|``. It is the accuracy floor of every purity and
    therefore of every label. It is a few ``1e-14`` on methane and ``1.6e-08``
    on benzene in 6-31G: the floor is set by the basis-set assembly, not by
    the point-group code.

the mixed |MO| count
    |MOs| with no definite irrep. A non-zero count means the |SCF| converged
    to a solution that breaks the point group, or that the geometry is not
    symmetric within :option:`symmetry pg_tol`, or that the basis does not
    support the group, or that the database held no |MOs| at all. The
    determinant and state layers then print probabilities instead of labels;
    see :ref:`symmetry-mixed-mos` below.

``qp_run symmetry_det`` does the same for the wave function stored in the
|EZFIO| database: families of determinants, every determinant label, the
per-irrep weights of every state and the gates D1-D6. Run it after a |CI|
calculation to inspect a state in detail.

Linear molecules
----------------

A collinear molecule takes a separate branch. Its group is not searched: every
plane through the axis maps each nucleus to itself, so the element key used for
the other groups cannot tell those operations apart, and the group is built
analytically instead. ``Cinfv`` is modelled as ``C31v`` (62 elements, 17
classes) and ``Dinfh`` as ``C31v x {E, sigma_h}`` (124, 34), renamed to the
chemical group before the table is opened. The labels are the usual
``Sigma_g+``, ``Sigma_u-``, ``Pi_g``, ``Delta_u``, ... because the sampled
rotation order, 31, is the largest the detector scans: the continuous
rotations act as the identity on a ``Sigma`` state and mix the two components
of a degenerate pair. A state of angular momentum :math:`|\Lambda| \ge 16`
would be indistinguishable from :math:`|\Lambda| = 15`; no basis set reaches
it.

The state table
---------------

:ref:`cis`, :ref:`cisd` and :ref:`fci` print one line per state after the
diagonalization::

    symmetry: state irreps: group C2v, 4 irreps, N_det = 81
    symmetry:    state      energy (Ha)       S^2  irrep   purity
    symmetry:         1    -75.9839974693    0.0000  A1   1.000
    symmetry:         2    -75.6375290051    0.0000  B1   1.000
    symmetry:         3    -75.5662848724    0.0000  A2   1.000

``purity`` is the largest irrep weight of the state, :math:`w_a(s)` of
:ref:`module_symmetry_det`, and is 1 for a state that belongs to a single
irrep. A purity below 1 means the state mixes irreps (a truncated or
non-converged vector, or a |SCF| solution that broke the symmetry); ``(none)``
in the irrep column is that same statement written out. The |S^2| column is
the value computed by the |CI| program, so a symmetry label and a spin label
are always read together.

The table is printed after every |CI| run whose labels are usable. There is
no printing switch; in particular the presence of the table never changes a
number.

.. _symmetry-mixed-mos:

Mixed orbitals: probabilities instead of labels
-----------------------------------------------

An |MO| basis can fail to be symmetry adapted. Three common causes: the
|SCF| converged to a solution that breaks the group, the geometry is
symmetric only within :option:`symmetry pg_tol`, or the |EZFIO| database
carries no ``mo_basis/mo_coef``, in which case QP2 falls back to the
orthonormalized |AOs|. Those |MOs| get index 0 in the |MO| report.

The determinant and state layers keep working, with a different reading:
where a label is not defined, they print the probability of each irrep. The
probabilities are exact for the determinant space that is stored. ``X_g`` is
block diagonal on the blocks found from the |MOs| themselves, so a family of
determinants keeps its per-block occupations under every group element; only
the one-line label is lost.

For example, a database without stored |MOs| (water, 6-31G*)::

    MO basis        : 18 MOs in 17 symmetry blocks, 2 without a definite irrep
      A mixed MO is a superposition: its two largest irrep weights are
        MO 8: A2  0.985564  B1  0.014436
        MO 9: B1  0.985564  A2  0.014436
      Determinant and state irreps below are probabilities wherever they are not 1.

A determinant that leaves a block partially occupied is printed with its
composition (this line comes from a fixture with a deliberately rotated
|MO| pair)::

    determinants with no definite irrep (probability of each irrep):
        det  family  |F|   purity  composition
          8       8    4  0.793412   A1 0.793412  B2 0.206588

Up to four weights above ``1e-4`` are shown per line, for at most 50
determinants. An unlabelled state gets the same treatment in the state
table, after the purity column.

On such a basis:

* state following in Davidson and :option:`symmetry irrep_states` are off:
  both need a label, and a distribution is not a label;
* ``qp_run symmetry_mo`` prints its full report and exits 1, because gate G5
  is "the |MO| set is symmetry adapted". That exit code is a verdict on the
  |MO| set, not a failure of the calculation;
* ``qp_run symmetry_det`` prints the probability tables and exits on its own
  gates.

The extra cost is bounded: the determinant table calls the diagonal of the
family projector, :math:`|G| \times |F|` matrix elements, for at most 50
rows. The state weights are the same exact projectors that are computed for
any other basis.

Keeping only some irreps
------------------------

:option:`symmetry irrep_states` restricts the states that the diagonalization
keeps to those whose irrep is listed. The value is a list of irrep names,
separated by commas, semicolons or blanks, matched without regard to case::

    "B1, A1"      "b1;a1"      "B1 A1"

Empty keeps every irrep. The default is the token ``None``, because |EZFIO|
cannot store a value-less default: it means the same empty list, and the
symmetry module turns it back into one. The keyword is part of the input
that ``qp_edit`` shows and writes, under the ``Symmetry`` section::

    irrep_states = B1, A1

The ``ezfio`` tool reaches the same field::

    source $QP_ROOT/quantum_package.rc
    ezfio set_file water.ezfio
    ezfio set symmetry irrep_states "B1, A1"
    qp_run cisd water.ezfio

The filter is a *narrowing* of what the calculation already selected. It is
applied after the |S^2| selection and, in Davidson, after state following:
it can only remove states, never add one, and a state with no definite irrep
is always kept, because a filter must not drop a state it cannot judge. The
kept states are moved to the front, in energy order.

When something is left out, the run says so::

    symmetry: irrep filter: state 3 is A2, not in irrep_states
    symmetry: irrep filter: only 1 state(s) in irrep_states for n_states = 2

The second line is the one that matters. The filter sees only the states that
survived the |S^2| selection, which are the lowest ``n_states`` of the
expected spin multiplicity. If the irreps you asked for are not among them,
the answer is not a different keyword but a larger window: raise
:option:`davidson_keywords n_states_diag`, the multiplier that sets how many
roots Davidson computes (:option:`determinants n_states` times
:option:`davidson_keywords n_states_diag`, at least 2), or
:option:`determinants n_states`, and diagonalize enough states for the wanted
irrep to appear. The second line is printed once per diagonalization call, so
a |CI| run that diagonalizes several times shows it several times.

Two situations are worth knowing:

* on the Lapack path (determinant spaces up to
  :option:`davidson_keywords n_det_max_full`), the filter is a stable
  permutation of the list produced by the |S^2| selection, which is
  "states of the expected |S^2|, then fill states". The filter promotes the
  allowed fill states too. Asking for ``A1`` with ``n_states = 3`` can
  therefore return two singlets and a triplet: it is the third ``A1`` state of
  the computed window, not the third singlet. The energies are unchanged;
* CSF-based Davidson cannot be filtered. With
  :option:`davidson_keywords csf_based` set, the program prints
  ``symmetry: irrep_states is ignored when csf_based is true`` and runs
  unchanged.

An unknown name is the one case where the program stops, because it is
malformed input and not a property of the calculation: it prints the names
that the detected group actually has, then stops.

State following
---------------

Davidson builds its guess vectors from the previous wave function, and
compares the irrep of each new vector with the irrep of the state it came
from. Two states of different irreps are never swapped, even when they are
nearly degenerate or cross along the iterations: the state of irrep ``A2``
stays the state of irrep ``A2``. This is automatic; there is no keyword for
it. When a state has no definite irrep the comparison is skipped for that
state and it follows the usual overlap criterion.

:option:`symmetry irrep_states` composes with state following as "follow the
states first, then keep the allowed ones". A filtered run therefore follows
the same states as the unfiltered one and only reorders them at the end,
which is what makes the reported energies directly comparable.

When the labels are unavailable
-------------------------------

The run continues in every one of these cases, and says why:

* **the basis does not support the group** -- the |AOs| of an atom, or a
  shell, do not map onto themselves under the operations. The report names
  the atom, the shell and the element;
* **an |MO| is not symmetry adapted** -- the reports switch to irrep
  probabilities, see :ref:`symmetry-mixed-mos`;
* **no character table** -- the message gives the path that was tried. The
  tables live in ``$QP_ROOT/data/Symmetry``, one file per group. A damaged
  file is reported and rejected, never partially used;
* **``C1``** -- a legitimate group with a one-row table. One irrep, ``A1``,
  one line per report, and the filter accepts ``A1`` only.

Keywords
--------

.. list-table::
   :header-rows: 1

   * - Keyword
     - Default
     - Meaning
   * - :option:`symmetry auto_symmetry`
     - |true|
     - Master switch of the whole detection.
   * - :option:`symmetry irrep_states`
     - ``None``
     - Irreps to keep in the |CI| diagonalization; ``None`` or an empty
       value keeps all.
   * - :option:`symmetry pg_tol`
     - ``1e-4``
     - Distance in bohr within which an atom maps onto another atom of the
       same element. The geometric tolerance of the detection.
   * - :option:`symmetry pg_thresh`
     - ``1e-6``
     - Weight threshold for the |MO| blocks, the determinant labels and the
       state purities.

Cost
----

The detection is a few dozen candidate directions and at most 128 group
elements: it is negligible next to one |SCF| cycle.

The remaining costs grow with the group, with the number of |MOs| and with
the number of determinants:

* |MO| labels: one matrix product ``C^T S M_g C`` per element, of order
  :math:`2 n_{MO} n_{AO}^2 + n_{MO}^2 n_{AO}` floating-point operations, plus
  a cache of :math:`|G| n_{MO}^2` doubles (0.80 MB for benzene in 6-31G, 24
  elements and 66 |MOs|). The |MO| report adds a few
  :math:`n_{irrep} n_{MO}` sums;
* determinant labels: one scan of the determinant list, then one projector
  per family of size greater than one, of order
  :math:`|G| |F|^2 \sum_b n_b^3/3` operations. For a symmetry-adapted basis
  the families are small and this is one pass over the list;
* state weights: those family projectors once, plus
  :math:`N_{states} \sum_F |F|^2 n_{irrep}` contractions. This is the exact
  table printed after a |CI| calculation. It never runs inside the |CI| loop:
  Davidson uses the cheap ``pg_det_state_labels``, which is one pass over the
  determinants per call;
* mixed basis: the determinant report adds one diagonal projector per printed
  row, :math:`|G| |F| \sum_b n_b^3/3` operations, for at most 50 rows. The
  state weights are unchanged.

Nothing runs for ``C1`` beyond the one-row report; everything can be
switched off with :option:`symmetry auto_symmetry`.
