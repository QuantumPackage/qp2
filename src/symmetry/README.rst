========
symmetry
========

Point-group detection and irrep labels for Quantum Package. The molecular
point group is *measured* from the nuclear geometry; the molecular orbitals
carry irreps; determinants and CI states are labelled with the same
character-table machinery.

The module is part of the core and always compiled in. It runs by itself at
the end of :ref:`module_hartree_fock` and :ref:`module_kohn_sham` when
:option:`auto_symmetry` is true (the default): it prints the group and the
number of MOs per irrep, and stores the MO labels in the EZFIO. Everything
else is on demand:

* :ref:`symmetry_mo` prints the complete group report and checks it;
* :ref:`symmetry_det` prints the determinant and state symmetry report and
  checks it;
* :ref:`cis`, :ref:`cisd` and :ref:`fci` print the irrep of every CI state
  after diagonalization;
* :option:`irrep_states` restricts a CI diagonalization to a list of
  irreps (see :ref:`module_symmetry_det`).

How the point group is found
============================

Nothing is read from the input and no table is selected before the geometry
has been measured. The geometry (:c:data:`nucl_coord`, :c:data:`nucl_charge`)
is the only input of the detection. This section describes the algorithm in
the order the code runs it.

The canonical frame
-------------------

The origin is the centre of mass, computed with ``element_mass`` exactly as
the rest of QP2 computes it. The inertia tensor about that origin is
diagonalized with ``dsyev``. Its eigenvectors are the principal axes; they
are sorted by ascending eigenvalue and each sign is fixed deterministically
(largest component positive), so two runs on the same geometry give
byte-identical directions.

``pg_is_planar`` is true when every nucleus lies within ``pg_tol`` of a
plane. The candidate normal is the inertia eigenvector with the largest
eigenvalue, which for a planar molecule is the normal to the plane.

The canonical frame is stored in ``pg_rot``: its rows are the canonical
``x``, ``y``, ``z`` unit vectors expressed in the input frame, so the
standard coordinates of a point are ``pg_rot . (r - pg_center)``.
Everything else in the module (candidate directions, rotation matrices,
element axes) stays in the input frame.

``z`` is the principal axis: the axis of the highest-order proper rotation.
If several distinct lines carry that order, the one keeping the most nuclei
on the axis is chosen, then the one closest to the most distinct inertia
axis, then the first candidate in the list. The tie-breaks are deterministic;
they do not change the group, only the labels of its axes.

``x`` is chosen among the branches of the standard Mulliken/Cotton
convention:

#. for a planar molecule whose principal axis lies in the plane, ``x`` is
   the plane normal and the molecular plane is ``yz`` (water, formaldehyde);
#. otherwise, if a two-fold proper axis perpendicular to ``z`` exists, ``x``
   is such an axis (the one with the most nuclei on it, then the one closest
   to the most distinct inertia axis);
#. otherwise, ``x`` is the direction of the lowest-index atom that is not on
   the ``z`` axis, projected on the plane perpendicular to ``z``;
#. if that fails too, a direction deduced from the inertia axes.

``y = z x x``. Cubic and icosahedral labels are invariant under relabelling
of their equivalent axes, so no "standard cube orientation" is imposed.

Candidate directions
--------------------

Symmetry axes and plane normals are searched in an over-complete list of
directions:

* every atom position relative to the centre;
* every midpoint of two atoms of the same element;
* the three inertia axes;
* the pairwise cross products of all the above.

Directions are normalized and deduplicated up to sign (tolerance ``1e-8`` on
the cosine), so the list holds distinct *lines*. The list is deliberately
generous: for the geometries of the data set, a symmetry axis is always
carried by an atom or a midpoint, or is an inertia axis or a cross product
of two of these.

For every candidate direction, the code tests proper rotations of order
``n = 32, 31, ..., 2`` by rotating the molecule and checking that every atom
lands on an atom of the same element within ``pg_tol``. The highest order
found is ``pg_n_max`` (up to 32, because ``D32h``/``D32d`` occur in the data
set); ``pg_n_axes_max`` counts the distinct lines carrying that order.

Enumeration of the group
------------------------

The group is *enumerated*, not looked up:

#. start with the identity;
#. add the proper rotations about the candidate axes (all powers of every
   order found), the reflections in the candidate planes, and the inversion
   when ``-1`` maps the atoms onto themselves;
#. close under composition until no new element appears.

Each element is accepted only if every atom maps onto an atom of the same
element within ``pg_tol``. The accepted map then defines an atom permutation
``pg_perm`` and an exact sign of the determinant ``pg_elem_det``.

Two elements are the same when they have **the same atom permutation and the
same determinant**. Matrices are never compared element by element. A
candidate direction read off a geometry that is symmetric only to
``pg_tol`` is tilted by ``pg_tol`` divided by the lever arm; a product of
two accepted elements then differs from the stored representative by far
more than any matrix threshold. Measured on NH3: the C3 axis is tilted by
``4.9e-9`` and 18 of the 36 products miss a ``1e-8`` matrix comparison,
while the permutation key is exact. The same key is used for deduplication,
for the multiplication table ``pg_elem_mul`` and for the classes.

``pg_geom_residual`` is the largest distance between an atom and its image
under an accepted element, in bohr. It is printed by the report and is the
honest measure of how symmetric the input geometry really is.

The largest group of the data set has 128 elements (``D32d``, ``D32h``);
``Ih`` has 120. The array bound is 160, so an unexpected larger group is
reported as a failure instead of being silently truncated.

Classes and invariants
----------------------

Classes are computed by direct conjugation on the atom permutations (with
the determinant to separate the identity from the reflection in the
molecular plane), never on the matrices, for the same reason. Class 1 is
the identity class.

The naming step uses only measured invariants:

* the group order ``|G|``;
* the highest proper-rotation order ``pg_n_max`` and the number of lines
  carrying it ``pg_n_axes_max``;
* the presence of the inversion ``pg_has_i``;
* the number of reflections whose normal is along ``z`` (``pg_n_sigma_h``,
  the horizontal plane) and the number of all other reflections
  (``pg_n_sigma_vtot``);
* the number of two-fold axes perpendicular to ``z`` (``pg_n_c2_perp``);
* whether the molecule is planar and whether ``z`` lies in the plane.

The "parallel / perpendicular to z" tests use a cosine threshold of ``1/2``
(45 degrees), not a tight tolerance: accepted elements carry the tilt of the
geometry, and in the groups of interest a line is either parallel or
perpendicular to ``z``, 90 degrees apart, so 45 degrees is the widest
possible margin.

Naming the group
----------------

``pg_group_name`` is produced by a decision tree on the invariants above.
It covers the 204 finite groups of the shipped data set; the tree is short
because two facts remove most of the cases:

* ``pg_n_axes_max`` counts the axes carrying a rotation of the highest
  order, so a ``D<n>`` group has ``pg_n_axes_max = 1`` for ``n >= 3`` (its
  perpendicular axes have order 2) and cannot be mistaken for a cubic group;
* in the ``n_max = 1`` branch the orientation of the single mirror is not
  usable, because ``z`` is the molecular-plane normal for a planar molecule
  but an inertia axis otherwise. ``C1``, ``Ci`` and ``Cs`` are therefore
  decided by ``|G|`` and by "there is exactly one reflection".

Five files of the data set (``C1h``, ``C1v``, ``D1``, ``D1d``, ``D1h``) are
degenerate aliases of ``Cs``, ``Cs``, ``C2``, ``C2h`` and ``C2v``: they are
the ``n = 1`` cases where the standard naming rules collapse two names onto
one group. The tree cannot produce them and their files are never opened.

A collinear molecule takes a separate branch. ``pg_is_linear`` flags it
before the enumeration: the moment about the axis vanishes and every
nucleus lies on the line through the centre. Its group is then built
analytically, as ``C31v`` (62 elements) or ``C31v x {E, sigma_h}`` (124)
when the horizontal plane maps the molecule onto itself, and reported as
``Cinfv`` or ``Dinfh``. The order of the sampled rotation is the largest
prime the detector supports, 31; the rotations that are not sampled act
inside a two-dimensional irrep, where they mix a degenerate pair only.
The shipped tables are the restrictions of the ``Cinfv``/``Dinfh`` irreps
to that group, so the labels are the usual ``Sigma_g+``, ``Pi_u``, ...;
``data/Symmetry/README.md`` gives their construction. The detector also
gives up when the invariants are inconsistent (a geometry that only looks
symmetric), rather than guessing a group.

The character table
-------------------

The group name is *the file name* of its character table. The file is read
from ``$QP_ROOT/data/Symmetry/<name>``. The first line of the file is
ignored; the name comes from the file name. The number of irreps is read
from the second line.

A file column is not a geometric class. In ``C3`` there are two geometric
classes (``C3`` and ``C3^2``) but one column; in ``T``, ``Td``, ``Oh`` and
``Ih`` a column is the union of several classes. The column matcher
therefore assigns to every file column the set of elements it stands for:
it parses the symbol (``E``, ``i``, ``sh``, ``sv``, ``sd``, ``s``, ``C2'``,
``C<n>``, ``C<n>^<k>``, ``S<n>``, ``S<n>^<k>``), selects the geometric
classes consistent with it, and closes the set under inversion (a class is
not closed: in ``C3`` the inverse of ``{C3}`` is ``{C3^2}``).

The axis predicates of the symbols use the canonical frame, never the raw
element vectors, which are in the input frame:

* ``C2'`` is the two-fold class containing ``x``;
* ``sv`` is the plane ``xz`` (normal along ``y``);
* ``sd`` and ``s`` are the plane ``yz`` (normal along ``x``);
* ``sh`` is the plane ``xy`` (normal along ``z``).

The parsed table is validated before use: the columns partition the group,
the first column is the identity, ``chi(E)`` is a positive integer, the
orthogonality sum ``sum_c |col(c)| chi_a(c) chi_b(c) = |G| mult_a delta_ab``
holds, and ``sum_a d_a^2/mult_a = |G|``. Row 1 must be the totally symmetric
irrep (every character 1), which is what makes irrep index 1 mean "totally
symmetric" everywhere. A damaged file prints the failing check and is
rejected; it is never silently used.

The gates of the MO report (G3 in particular) are the authoritative test of
the column assignment: a swap of ``sv`` and ``sd`` in an Abelian group keeps
every orthogonality sum intact but changes the number of MOs per irrep.

Why the labels match the literature
-----------------------------------

The axis convention above is the Mulliken/Cotton one. With it, a ``C2v``
molecule such as water has

* ``z`` = ``A1`` (the C2 axis),
* ``x`` = ``B1`` (the plane normal),
* ``y`` = ``B2`` (in the molecular plane),

and ``D2h`` molecules have ``x = B3u``, ``y = B2u``, ``z = B1u``. Getting
``sv`` and ``sd`` backwards swaps ``B1`` and ``B2`` in every Abelian group;
that is why the convention is enforced in code and why the MO gates test
the actual matrices rather than the table alone.

MO irreps
=========

The AO and MO representation matrices, the symmetry blocks and the MO labels
are described in :ref:`module_symmetry_mo`. In one paragraph: the module
builds, for every group element, the matrix ``X_g = C^T S M_g C`` of the
element in the MO basis, groups the MOs into symmetry blocks (connected
components of ``|X_g(j,i)|^2 > pg_thresh`` for some element), and evaluates
the projector

.. math::

   P^a = \frac{d_a}{m_a |G|} \sum_c \chi_a(c) U_c

on every MO. ``pg_mo_irrep(i)`` is the argmax of ``p_a(i) = <i|P^a|i>``, or
0 when ``max_a p_a(i) < 1 - pg_thresh``: a mixed MO gets no best-guess
label. ``m_a`` is 2 when the table row is the sum of a complex-conjugate
pair (``C3``, ``C10``, ...) and 1 otherwise; ``d_a`` is the dimension.

The accuracy floor of these labels is not machine precision. ``M_g`` is
built from the geometry and from the AO integrals, so both the geometric
tolerance ``pg_tol`` and the assembly defect of the basis enter. For
geometries symmetric to machine precision in 6-31G, the gate
``M_g^T S M_g = S`` measures ``5.7e-17`` on water, ``5.6e-14`` on methane and
``4.0e-09`` on benzene; the MO-metric defect ``delta_MO`` is ``7.0e-14`` on
methane and ``1.6e-08`` on benzene. The default ``pg_thresh = 1e-6`` is the
weight threshold, one to two orders of magnitude above those floors. The
full report (``qp_run symmetry_mo``) prints the measured floor next to the
number of mixed MOs.

Determinants and states
=======================

Determinant and state symmetry are described in :ref:`module_symmetry_det`.
The short version:

* determinants are grouped into *families* -- the determinants with the same
  per-block spin counts;
* ``<D'|U(g)|D>`` is a product of determinants of ``X_g`` submatrices, one per
  symmetry block;
* the same projector formula gives ``T^a`` on every family,
  ``pg_det_irrep`` and ``pg_det_purity``; a determinant is pure if and only
  if its family has one member, and it is unlabelled (irrep 0) when no
  projector weight reaches ``1 - pg_thresh``;
* the weight of irrep ``a`` in a CI state is the exact diagonal matrix
  element ``w_a(s) = sum_{D,D'} c_D(s) T^a_{DD'} c_D'(s)``, so degenerate
  shells and off-diagonal family terms are included. The report prints the
  energy, the :math:`S^2` and the label of every state after a CI
  calculation.

On a basis that is not symmetry adapted the projectors still close exactly
on the stored determinant space, because ``X_g`` is block diagonal on the
blocks found from the MOs themselves. The labels then become probabilities:
the determinant and state reports print the weight of every irrep, and state
following and the :option:`irrep_states` filter stay off. See
:ref:`module_symmetry_det`.

The CI programs use the same labels twice. In Davidson, the labels of the
current vectors are compared with the labels of the previous wave function
so that state following does not jump between irreps. The
:option:`irrep_states` keyword restricts the states that are kept: it is a
narrowing of the states that already passed the :math:`S^2` selection and
the state-following attribution, never a replacement for them. A state with
no definite irrep is always kept (a filter must not drop a state it cannot
judge).

Keywords
========

.. list-table::
   :header-rows: 1

   * - Keyword
     - Default
     - Meaning
   * - :option:`pg_tol`
     - ``1e-4``
     - Distance (bohr) within which an atom is considered to map onto another
       atom of the same element. It is the geometric tolerance of the whole
       detection.
   * - :option:`pg_thresh`
     - ``1e-6``
     - Weight threshold. An off-diagonal ``|X_g(j,i)|^2`` above it connects
       two MOs into one block; a projector weight above ``1 - pg_thresh``
       makes an MO or a determinant label definite.
   * - :option:`auto_symmetry`
     - ``true``
     - Master switch. When false, no provider of the detection chain is
       built, the SCF/DFT report is silent, and nothing is labelled.
   * - :option:`irrep_states`
     - ``All``
     - Comma-, semicolon- or blank-separated list of irrep names to keep in
       the CI diagonalization, for example ``B1, A1``. ``All`` (the
       default), the legacy token ``None`` and the empty string all keep
       every irrep; the provider turns those tokens back into the empty
       list.

.. note::

   :option:`irrep_states` appears in ``qp_edit``, in the ``Symmetry`` block,
   and is shown there as ``irrep_states = All`` by default. It can also be
   written with the ``ezfio`` tool::

       source $QP_ROOT/quantum_package.rc
       qp_run scf water.ezfio            # or any run that creates the EZFIO
       ezfio set_file water.ezfio
       ezfio set symmetry irrep_states "B1, A1"
       qp_run cisd water.ezfio

   An unknown irrep name prints the available names and stops: malformed
   user input is the one place where stopping is allowed.

What is stored
==============

The only EZFIO field written by the module is
``mo_basis/mo_symmetry`` (the provider :c:data:`mo_symmetry` of QP2). It
holds the irrep index of every MO, with the convention "same integer = same
irrep". It is written at the end of SCF/DFT only when *every* MO has a
definite irrep; a mixed MO has index 0, and writing 0 would silently declare
all mixed MOs to belong to one irrep. The labels are recomputed from the
geometry and the MOs at every run, so the stored field is a record, not an
input constraint.

Cost
====

* Detection: a few dozen candidate directions, ``|G| <= 128`` elements.
  Negligible compared with an SCF cycle.
* MO labels: one AO transformation per element (``|G|`` matrix products of
  the AO sizes) and ``|G|`` matrices of ``mo_num^2`` doubles when the
  determinant layer keeps ``pg_xg`` (0.8 MB for benzene in 6-31G, 66 MOs,
  ``|G| = 24``).
* State table: one pass over the determinants per printed state, plus one
  projector per family with more than one member. It is printed after every
  CI run whose labels are usable; for ``C1`` every state is ``A``, so the
  table only confirms that.

Limits and failure modes
========================

* **Linear molecules use a finite model.** ``Cinfv`` and ``Dinfh`` are
  sampled as ``C31v`` and ``C31v x {E, sigma_h}``. A rotation about the
  axis that is not sampled only mixes a degenerate pair, so it changes no
  label; the one thing the model cannot resolve is a state of angular
  momentum ``|Lambda| >= 16``, which aliases onto ``Lambda = 15`` (no
  basis set reaches it). The two tables carry 10 decimals because their
  two-dimensional characters are irrational; see
  ``data/Symmetry/README.md``.
* **The geometry is not symmetrized.** A molecule that is symmetric only
  within ``pg_tol`` is treated as symmetric, and every label inherits that
  geometric error. Use the reported ``pg_geom_residual`` and the
  ``delta_MO`` of the report to judge a run.
* **A mixed MO basis stops the labels, not the program.** Molecular orbitals
  that are not symmetry adapted (for example because the SCF mixed
  degenerate orbitals) are reported as mixed, their index is 0, and no
  state table can be trusted. Symmetrize the geometry, or converge the SCF
  more tightly.
* **Family cap.** A determinant family is enumerated only up to 1000
  members, and a symmetry block may hold at most 12 occupied MOs in the
  determinant projector. Beyond that the module reports which determinant
  it refused and labels nothing; it never hangs or silently truncates.
* **No table.** A group name without a readable table in the data directory
  is reported with the path that was tried, and the labels are dropped.
* **C1 is a legitimate group.** It has a one-row table and one irrep; the
  report prints one line, and the filter accepts ``A1`` only.

Data files
==========

``$QP_ROOT/data/Symmetry`` holds one character table per point group: 211
files, being the 204 finite groups reachable by the naming tree, the five
degenerate aliases listed above, and the two linear groups (``Cinfv``,
``Dinfh``), whose bodies are the ``C31v`` and ``D31h`` tables. The files
were generated by
the author of this module and are original data, not a copy of a printed
compilation; the provenance is recorded in the directory's ``README.md``.
The set is fixed data: it is extended by adding a file, never regenerated.
Every file is parsed and validated at use, so a damaged file is reported
and rejected instead of being used.

Conventions and design notes
============================

The comments in the sources carry the design rationale: the equations, the
gates and the measurements behind each choice. They are the reference for
anyone changing the module. The pages below restate what is needed to read
and use the code.

How to check the machinery
==========================

``qp_run symmetry_mo DIR`` prints the frame, the group, its classes, the
table, every MO label and gates G1-G6. ``qp_run symmetry_det DIR`` prints
the determinant table, the per-state irrep weights and gates D1-D6. Both
exit 1 after printing everything if a gate fails, so they can be chained in
a test script. The SCF run only prints the summary; the two programs are the
full diagnostic.
