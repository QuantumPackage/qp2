===========
symmetry_mo
===========

Molecular-orbital irreps. This module turns the point group found by
:ref:`module_symmetry` into labels for the molecular orbitals, writes them
to the EZFIO at the end of SCF/DFT, and provides the labels used by the
determinant and state layer (:ref:`module_symmetry_det`).

Needed Modules
==============
* :ref:`module_symmetry` -- point group, character tables
* :ref:`module_mo_basis`, :ref:`module_ao_basis`,
  :ref:`module_ao_one_e_ints` -- the MO coefficients and the AO overlap
* :ref:`module_utils`

What it adds
============

``pg_xg`` and the AO transformation
-----------------------------------

For a group element :math:`g`, the module builds its matrix
:math:`M_g` in the AO basis. Contracted shells are permuted by the atom
permutation of the element (``pg_shell_target``); inside a shell, the
block is the symmetric power :math:`\mathrm{Sym}^l(R_g)` of the rotation,
rescaled by the ratio of the radial factors of the primitives. Nothing is
hardcoded about the angular convention: the monomial representation is
expanded from :math:`R_g` directly.

The AO basis must support the group: it must be invariant under it, and the
coefficients of a shell must factorize as :math:`A_j B_p` (a radial factor
times an angular factor). ``pg_basis_ierr`` is 0 when this holds and 1, 2 or
3 otherwise, with ``pg_basis_msg`` naming the atom, shell, angular momentum
and element. The diagnostic is reported, never thrown from a provider.

The matrix of the element in the MO basis is

.. math::

   X_g = C^T S M_g C ,

where :math:`C` is the MO coefficient matrix and :math:`S` the AO overlap.
It is built once per element by ``pg_get_mo_rep``; the determinant layer
keeps the set in ``pg_xg``.

Symmetry blocks and MO labels
-----------------------------

A *symmetry block* is a connected component of the graph on the MOs with an
edge :math:`(i,j)` when :math:`|X_g(j,i)|^2 >` :option:`symmetry pg_thresh` for some
element :math:`g`. The test is on a weight, the same number and the same
units as the label cut, because the assembly defect is a weight defect: the
spurious amplitude between two MOs of different irreps is the square root of
the defect. One block is one irrep space; two copies of the same irrep at
different energies stay separate blocks.

The projector on irrep :math:`a` is the same formula everywhere in the
module:

.. math::

   P^a = \frac{d_a}{m_a |G|} \sum_c \chi_a(c) U_c ,

with :math:`d_a` the dimension (``pg_irrep_dim``), :math:`m_a` 2 when the
table row is the sum of a complex-conjugate pair and 1 otherwise
(``pg_irrep_mult``), and :math:`\chi_a(c)` the characters of the file column
:math:`c`. The normalisation :math:`d_a/(m_a |G|)` is stored once, in
``pg_proj_coef(a)``: the MO weights, the determinant and state projector,
and the determinant diagonal kernel all read that array instead of
recomputing the ratio, so the convention cannot drift between them.
The weight of irrep :math:`a` in MO :math:`i` is
``pg_mo_weight(i,a) = p_a(i) = <i|P^a|i>``. Then

* ``pg_mo_irrep(i)`` is the argmax of :math:`p_a(i)`;
* it is 0 when :math:`\max_a p_a(i) < 1 -` :option:`symmetry pg_thresh`: a mixed MO is
  not given a best guess;
* ``pg_mo_purity(i)`` is :math:`\max_a p_a(i)`, equal to 1 for a symmetry
  adapted MO and also 1 for every partner of a degenerate block, where
  :math:`P^a` is the identity;
* ``pg_irrep_weight(a)`` is the row sum of the weights, i.e. the number of
  MOs carrying irrep :math:`a`. It is never recomputed from a second trace
  formula: the gates test the very numbers the report prints.

The report at the end of SCF/DFT
--------------------------------

``pg_report_mo_symmetry`` is called at the end of :ref:`module_hartree_fock`
and :ref:`module_kohn_sham`. It prints the point group and the number of MOs
per irrep, and writes ``mo_symmetry`` when every MO has a definite irrep.
It never stops and prints nothing when :option:`symmetry auto_symmetry` is false.

``pg_write_mo_symmetry`` writes ``mo_basis/mo_symmetry``: the integer irrep
index of each MO, with the convention "same integer = same irrep", and
``mo_basis/mo_symmetry_name``: the irrep name of each MO, so that
:ref:`qp_edit` can label the MOs without rebuilding the character table. A
mixed MO has index 0 and nothing is written in that case, because writing 0
would declare all mixed MOs to belong to one irrep. The in-memory provider
is invalidated with ``soft_touch`` after the write, so the current process
and the EZFIO agree.

``pg_print_mo_energies(e, n)`` is called by the SCF, KS and range-separated
KS programs, at the program level because ``scf_utils`` must not depend on
this module. It prints the orbital energies in the order of the MOs, with the
irrep label of each MO in a third column. When :option:`symmetry
auto_symmetry` is false, when the detection gave up, or when some MO has no
definite irrep, it prints the historical ``Eigenvalues`` table instead, with
no irrep column and no extra message: ``pg_report_mo_symmetry`` prints the
reason right after. ``scf_utils`` stores the same numbers in
``mo_basis/mo_energy`` at the end of the SCF, in the order in which
``save_mos`` writes the MOs; that array is a snapshot and becomes stale if
the MOs are rotated afterwards.

Guards
------

``pg_symmetry_ready()`` is true when the group was detected, its table is
usable and the AO basis supports it. ``pg_ao_ok`` is the same statement for
the basis alone. ``pg_mo_labels_ok()`` is true when every MO has a definite
irrep. These functions never stop: they are what keeps a failed detection
from aborting a CI run. A provider of this module is only entered when its
guard has already passed (IRPF90 hoists a provider build to the entry of the
procedure that names it, so the guard is always written in the caller, never
in the worker).

Report program
==============

``qp_run symmetry_mo DIR`` (source: ``symmetry_mo.irp.f``, ``program
symmetry``) prints:

* the canonical frame and the measured geometry residual;
* the group, its order, its classes and the character table;
* the AO-basis gates and ``delta_MO``;
* every MO label and the gates G1-G6.

With ``PG_WRITE_MO_SYMMETRY=1`` in the environment it also writes
``mo_basis/mo_symmetry``. Nothing is written when a gate failed or when an
MO has no definite irrep, so a written field is always a checked partition.

Gates
=====

Every gate prints its measured value and its bound; a failing gate is
reported and the run continues, so one run shows every number, and the
program exits 1 at the very end.

.. list-table::
   :header-rows: 1

   * - Gate
     - Test
     - Bound
   * - G-iso
     - :math:`M_g^T S M_g = S` for every element, with the raw AO overlap
     - ``1e-7``
   * - G1
     - :math:`C^T S C = 1` (orthonormality of the MO set)
     - ``1e-7``
   * - G2
     - :math:`X_g X_h = X_{gh}`, products taken from the permutation table
     - ``1e-7``
   * - G3
     - every :math:`\mathrm{tr}\,P^a` is an integer multiple of :math:`d_a`
     - ``1e-6``
   * - G4
     - :math:`\sum_a \mathrm{tr}\,P^a =` ``mo_num``
     - ``1e-6``
   * - G5
     - :math:`0 \le p_a(i) \le 1`, :math:`\sum_a p_a(i) = 1`,
       :math:`\max_a p_a(i) = 1`, with the window ``max(1e-8, 10*delta_MO)``
     - window
   * - G6
     - :math:`P^a P^a v = P^a v` and :math:`\sum_a P^a v = v` on a fixed
       pseudo-random vector
     - ``1e-6``

``delta_MO = max_g |X_g^T X_g - 1|`` is printed but not gated: it is G-iso
transported into the MO metric, and it is the accuracy floor of every purity
of this module. G3 is the authoritative test of the column assignment of the
character table: a swap of ``sv`` and ``sd`` in an Abelian group keeps every
orthogonality sum intact but changes the number of MOs per irrep.

The ``1e-6`` bounds are not machine precision and must not be lowered:
the representation matrices are built from *representatives* of the group
elements, which are only as accurate as the geometry tolerance. Measured on
NH3, ``X_g X_gh`` is off by ``2.3e-8`` at ``pg_tol = 1e-4`` bohr, and
:math:`P^a` is a linear combination of the :math:`X_g`, so
:math:`||P^2 - P||` grows with the same floor. A wrong column assignment, a
wrong :math:`m_a` or a missing irrep is O(1) in these gates.

Mixed MOs
=========

A molecular orbital that is not symmetry adapted gets index 0 and purity
below 1. This is not a bug in the labels: it means the SCF converged to a
solution that breaks the point group, that the geometry is not symmetric
within :option:`symmetry pg_tol`, or that the EZFIO database holds no
``mo_basis/mo_coef`` and QP2 falls back to the orthonormalized atomic
orbitals. The report names the mixed MOs with the two largest weights of
each, the number of blocks and ``delta_MO``, so that the cause can be
identified.

The determinant and state layers keep working on such a basis. ``X_g`` is
block diagonal on the blocks of this module by construction, so the
per-block occupations of a determinant never change under the group: the
family projector stays exact, and the irrep of a determinant or of a state
is an exact probability for the stored determinant space.
:ref:`module_symmetry_det` prints those probabilities. Only the label is
lost, never the symmetry information.

What a mixed basis does disable is everything that needs a label instead of
a distribution: state following in Davidson and the :option:`irrep_states`
filter are off, because they call ``pg_labels_usable``, which tests
``pg_mo_labels_ok``; and ``mo_basis/mo_symmetry`` is not written. No program
stops.

``qp_run symmetry_mo`` is a diagnostic: on such a basis it prints the whole
report and then exits 1, because gate G5 is precisely "the MO set is
symmetry adapted". That exit code is the verdict on the MO set, and is
expected on a fixture that is meant to exercise the mixed path. The CI
programs and ``qp_run symmetry_det`` are unaffected; ``symmetry_det`` exits
on its own gates.

The cost of the mixed path is one diagonal projector per printed
determinant, and the usual family projectors for the states; both are
estimated in :ref:`module_symmetry_det`.
