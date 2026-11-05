.. _module_symmetry: 
 
.. program:: symmetry 
 
.. default-role:: option 
 
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
     - ``None``
     - Comma-, semicolon- or blank-separated list of irrep names to keep in
       the CI diagonalization, for example ``B1, A1``. ``None`` (the
       default) and the empty string both keep every irrep; the provider
       turns ``None`` back into the empty list.

.. note::

   :option:`irrep_states` is declared with ``interface: ezfio`` only, so it
   does **not** appear in ``qp_edit``. Write it with the ``ezfio`` tool::

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
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: pg_tol
 
    Geometries are considered symmetric when every atom maps onto an atom of the same element within this distance (bohr).
 
    Default: 1.e-4
 
.. option:: pg_thresh
 
    Threshold on representation-matrix elements used to decide whether an MO belongs to a symmetry block, and whether a determinant has a definite irrep.
 
    Default: 1.e-6
 
.. option:: auto_symmetry
 
    If true, detect the molecular point group at the end of SCF/DFT and label the MOs.
 
    Default: True
 
.. option:: irrep_states
 
    Comma-, semicolon- or blank-separated list of irrep names to keep in the CI diagonalization, for example "A1,B1". The token None, the default, keeps all irreps.
 
    Default: None
 
 
Providers 
--------- 
 
.. c:var:: bool auto_symmetry


    File : :file:`symmetry/settings.irp.f`

    .. code:: fortran

        logical	:: auto_symmetry	


    When true, the molecular point group is detected at the end of SCF/DFT,
    the MOs are labelled and ``mo_basis/mo_symmetry`` is written. Default
    .true.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_irrep_filter`

 
.. c:var:: char irrep_states


    File : :file:`symmetry/settings.irp.f`

    .. code:: fortran

        character*(256)	:: irrep_states	


    Comma-, semicolon- or blank-separated list of irrep names to keep in
    the CI diagonalization. Empty means all irreps; the
    qp_edit default is the token "None", which :c:func:`pg_normalize_irrep_states`
    turns into the empty string.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_n_irrep_allowed`

 
.. c:var:: double pg_center


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_center	(3)


    Centre of mass: origin of the canonical frame.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`center_of_mass`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_geom_residual`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_is_planar`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_max`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`

 
.. c:var:: double pg_char


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        integer	:: pg_n_col	
        integer	:: pg_n_irrep	
        character*(32), allocatable	:: pg_irrep_name	(pg_n_elem)
        double precision, allocatable	:: pg_char	(pg_n_elem,pg_n_elem)
        character*(8), allocatable	:: pg_data_sym	(pg_n_elem)
        integer, allocatable	:: pg_col_start	(pg_n_elem)
        integer, allocatable	:: pg_col_len	(pg_n_elem)
        integer, allocatable	:: pg_col_member	(pg_n_elem)
        integer, allocatable	:: pg_irrep_dim	(pg_n_elem)
        integer, allocatable	:: pg_irrep_mult	(pg_n_elem)


    Character table of :c:data:`pg_group_name`, with every file column
    matched to the set of group elements it stands for, and validated.
    
    Main outputs, for mo_irrep.irp.f:
    
    ``pg_irrep_name(a)``
        name of irrep a, as in the file (A1, B2u, T1g, E')
    
    ``pg_irrep_dim(a)``
        chi_a(E) = nint(pg_char(a,1))
    
    ``pg_irrep_mult(a)``
        1 for a genuine real irrep, 2 when the row is the sum of a
        complex-conjugate pair (C3, C10, ...). The MO projector is
        then (d_a/mult_a) sum_col chi_a(col) U_col.
    
    ``pg_char(a,c)``
        characters. Column c is file column c: the matching does not permute
        anything, so this is the parsed matrix verbatim.
    
    ``pg_col_len(c)``
        number of ELEMENTS of column c. This replaces the prototype's
        pg_class_size everywhere.
    
    ``pg_col_member(...)``
        the element indices, columns concatenated; column c is
        pg_col_member(pg_col_start(c) : ...).
    
    Validation: the columns partition the group, the first column is the
    identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
    orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
    sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
    :c:data:`pg_detection_failed`: the caller of this provider must give up.
    
    These identities do NOT detect a swap of two columns that have equal
    sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
    intact): the authoritative gate for the column *labelling* is the
    projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_class


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer	:: pg_n_class	
        integer, allocatable	:: pg_class	(pg_n_elem)


    Classes by conjugation. Element 1 is the identity, so class 1 is the
    identity class. Direct orbit enumeration (:math:`|G|` <= 120).
    
    The conjugation is done on the atom permutations, not on the
    matrices: the stored matrices carry the :c:data:`pg_tol`-sized tilt of the
    candidate directions, so a conjugated matrix is not within 1e-8 of any
    stored one. The permutation (with the determinant, to separate the
    identity from the reflection in the molecular plane) is the exact
    invariant used by :c:func:`pg_add_element` too. Direct orbit enumeration,
    :math:`|G|` <= 128 so :math:`|G|^2` = 16384 comparisons.
    
    Linear molecules are the exception: all their permutations are the
    identity, so the conjugation is done on the analytic matrices there
    (:c:func:`pg_find_by_mat`).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_col`

 
.. c:var:: int pg_class_rep


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_class_size	(pg_n_class)
        integer, allocatable	:: pg_class_rep	(pg_n_class)


    Size of every class and its lowest element index.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_group_name`
       * :c:data:`pg_n_col`

 
.. c:var:: int pg_class_size


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_class_size	(pg_n_class)
        integer, allocatable	:: pg_class_rep	(pg_n_class)


    Size of every class and its lowest element index.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_group_name`
       * :c:data:`pg_n_col`

 
.. c:var:: int pg_col_len


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        integer	:: pg_n_col	
        integer	:: pg_n_irrep	
        character*(32), allocatable	:: pg_irrep_name	(pg_n_elem)
        double precision, allocatable	:: pg_char	(pg_n_elem,pg_n_elem)
        character*(8), allocatable	:: pg_data_sym	(pg_n_elem)
        integer, allocatable	:: pg_col_start	(pg_n_elem)
        integer, allocatable	:: pg_col_len	(pg_n_elem)
        integer, allocatable	:: pg_col_member	(pg_n_elem)
        integer, allocatable	:: pg_irrep_dim	(pg_n_elem)
        integer, allocatable	:: pg_irrep_mult	(pg_n_elem)


    Character table of :c:data:`pg_group_name`, with every file column
    matched to the set of group elements it stands for, and validated.
    
    Main outputs, for mo_irrep.irp.f:
    
    ``pg_irrep_name(a)``
        name of irrep a, as in the file (A1, B2u, T1g, E')
    
    ``pg_irrep_dim(a)``
        chi_a(E) = nint(pg_char(a,1))
    
    ``pg_irrep_mult(a)``
        1 for a genuine real irrep, 2 when the row is the sum of a
        complex-conjugate pair (C3, C10, ...). The MO projector is
        then (d_a/mult_a) sum_col chi_a(col) U_col.
    
    ``pg_char(a,c)``
        characters. Column c is file column c: the matching does not permute
        anything, so this is the parsed matrix verbatim.
    
    ``pg_col_len(c)``
        number of ELEMENTS of column c. This replaces the prototype's
        pg_class_size everywhere.
    
    ``pg_col_member(...)``
        the element indices, columns concatenated; column c is
        pg_col_member(pg_col_start(c) : ...).
    
    Validation: the columns partition the group, the first column is the
    identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
    orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
    sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
    :c:data:`pg_detection_failed`: the caller of this provider must give up.
    
    These identities do NOT detect a swap of two columns that have equal
    sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
    intact): the authoritative gate for the column *labelling* is the
    projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_col_member


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        integer	:: pg_n_col	
        integer	:: pg_n_irrep	
        character*(32), allocatable	:: pg_irrep_name	(pg_n_elem)
        double precision, allocatable	:: pg_char	(pg_n_elem,pg_n_elem)
        character*(8), allocatable	:: pg_data_sym	(pg_n_elem)
        integer, allocatable	:: pg_col_start	(pg_n_elem)
        integer, allocatable	:: pg_col_len	(pg_n_elem)
        integer, allocatable	:: pg_col_member	(pg_n_elem)
        integer, allocatable	:: pg_irrep_dim	(pg_n_elem)
        integer, allocatable	:: pg_irrep_mult	(pg_n_elem)


    Character table of :c:data:`pg_group_name`, with every file column
    matched to the set of group elements it stands for, and validated.
    
    Main outputs, for mo_irrep.irp.f:
    
    ``pg_irrep_name(a)``
        name of irrep a, as in the file (A1, B2u, T1g, E')
    
    ``pg_irrep_dim(a)``
        chi_a(E) = nint(pg_char(a,1))
    
    ``pg_irrep_mult(a)``
        1 for a genuine real irrep, 2 when the row is the sum of a
        complex-conjugate pair (C3, C10, ...). The MO projector is
        then (d_a/mult_a) sum_col chi_a(col) U_col.
    
    ``pg_char(a,c)``
        characters. Column c is file column c: the matching does not permute
        anything, so this is the parsed matrix verbatim.
    
    ``pg_col_len(c)``
        number of ELEMENTS of column c. This replaces the prototype's
        pg_class_size everywhere.
    
    ``pg_col_member(...)``
        the element indices, columns concatenated; column c is
        pg_col_member(pg_col_start(c) : ...).
    
    Validation: the columns partition the group, the first column is the
    identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
    orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
    sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
    :c:data:`pg_detection_failed`: the caller of this provider must give up.
    
    These identities do NOT detect a swap of two columns that have equal
    sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
    intact): the authoritative gate for the column *labelling* is the
    projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_col_start


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        integer	:: pg_n_col	
        integer	:: pg_n_irrep	
        character*(32), allocatable	:: pg_irrep_name	(pg_n_elem)
        double precision, allocatable	:: pg_char	(pg_n_elem,pg_n_elem)
        character*(8), allocatable	:: pg_data_sym	(pg_n_elem)
        integer, allocatable	:: pg_col_start	(pg_n_elem)
        integer, allocatable	:: pg_col_len	(pg_n_elem)
        integer, allocatable	:: pg_col_member	(pg_n_elem)
        integer, allocatable	:: pg_irrep_dim	(pg_n_elem)
        integer, allocatable	:: pg_irrep_mult	(pg_n_elem)


    Character table of :c:data:`pg_group_name`, with every file column
    matched to the set of group elements it stands for, and validated.
    
    Main outputs, for mo_irrep.irp.f:
    
    ``pg_irrep_name(a)``
        name of irrep a, as in the file (A1, B2u, T1g, E')
    
    ``pg_irrep_dim(a)``
        chi_a(E) = nint(pg_char(a,1))
    
    ``pg_irrep_mult(a)``
        1 for a genuine real irrep, 2 when the row is the sum of a
        complex-conjugate pair (C3, C10, ...). The MO projector is
        then (d_a/mult_a) sum_col chi_a(col) U_col.
    
    ``pg_char(a,c)``
        characters. Column c is file column c: the matching does not permute
        anything, so this is the parsed matrix verbatim.
    
    ``pg_col_len(c)``
        number of ELEMENTS of column c. This replaces the prototype's
        pg_class_size everywhere.
    
    ``pg_col_member(...)``
        the element indices, columns concatenated; column c is
        pg_col_member(pg_col_start(c) : ...).
    
    Validation: the columns partition the group, the first column is the
    identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
    orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
    sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
    :c:data:`pg_detection_failed`: the caller of this provider must give up.
    
    These identities do NOT detect a swap of two columns that have equal
    sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
    intact): the authoritative gate for the column *labelling* is the
    projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: char pg_data_sym


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        integer	:: pg_n_col	
        integer	:: pg_n_irrep	
        character*(32), allocatable	:: pg_irrep_name	(pg_n_elem)
        double precision, allocatable	:: pg_char	(pg_n_elem,pg_n_elem)
        character*(8), allocatable	:: pg_data_sym	(pg_n_elem)
        integer, allocatable	:: pg_col_start	(pg_n_elem)
        integer, allocatable	:: pg_col_len	(pg_n_elem)
        integer, allocatable	:: pg_col_member	(pg_n_elem)
        integer, allocatable	:: pg_irrep_dim	(pg_n_elem)
        integer, allocatable	:: pg_irrep_mult	(pg_n_elem)


    Character table of :c:data:`pg_group_name`, with every file column
    matched to the set of group elements it stands for, and validated.
    
    Main outputs, for mo_irrep.irp.f:
    
    ``pg_irrep_name(a)``
        name of irrep a, as in the file (A1, B2u, T1g, E')
    
    ``pg_irrep_dim(a)``
        chi_a(E) = nint(pg_char(a,1))
    
    ``pg_irrep_mult(a)``
        1 for a genuine real irrep, 2 when the row is the sum of a
        complex-conjugate pair (C3, C10, ...). The MO projector is
        then (d_a/mult_a) sum_col chi_a(col) U_col.
    
    ``pg_char(a,c)``
        characters. Column c is file column c: the matching does not permute
        anything, so this is the parsed matrix verbatim.
    
    ``pg_col_len(c)``
        number of ELEMENTS of column c. This replaces the prototype's
        pg_class_size everywhere.
    
    ``pg_col_member(...)``
        the element indices, columns concatenated; column c is
        pg_col_member(pg_col_start(c) : ...).
    
    Validation: the columns partition the group, the first column is the
    identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
    orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
    sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
    :c:data:`pg_detection_failed`: the caller of this provider must give up.
    
    These identities do NOT detect a swap of two columns that have equal
    sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
    intact): the authoritative gate for the column *labelling* is the
    projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: bool pg_detection_failed


    File : :file:`symmetry/settings.irp.f`

    .. code:: fortran

        logical	:: pg_detection_failed	


    Set when the point-group detection gave up. Monotonic: only ever set,
    never cleared, so a second attempt cannot pretend it succeeded.

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_elem_inv`
       * :c:data:`pg_elem_mul`
       * :c:data:`pg_group_name`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_irrep_filter`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_max`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`

 
.. c:var:: double pg_elem_angle


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_rot_mat	(3,3,pg_n_elem)
        integer, allocatable	:: pg_perm	(nucl_num,pg_n_elem)
        double precision, allocatable	:: pg_elem_angle	(pg_n_elem)
        double precision, allocatable	:: pg_elem_axis	(3,pg_n_elem)
        integer, allocatable	:: pg_elem_det	(pg_n_elem)
        integer, allocatable	:: pg_elem_kind	(pg_n_elem)


    Rotation matrices (input frame), induced atom permutations, angle,
    axis (plane normal for a reflection), determinant sign, and the
    integer (det,m,k) key of every group element.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_elem_inv`
       * :c:data:`pg_elem_mul`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_has_i`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_col`
       * :c:data:`pg_shell_target`

 
.. c:var:: double pg_elem_axis


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_rot_mat	(3,3,pg_n_elem)
        integer, allocatable	:: pg_perm	(nucl_num,pg_n_elem)
        double precision, allocatable	:: pg_elem_angle	(pg_n_elem)
        double precision, allocatable	:: pg_elem_axis	(3,pg_n_elem)
        integer, allocatable	:: pg_elem_det	(pg_n_elem)
        integer, allocatable	:: pg_elem_kind	(pg_n_elem)


    Rotation matrices (input frame), induced atom permutations, angle,
    axis (plane normal for a reflection), determinant sign, and the
    integer (det,m,k) key of every group element.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_elem_inv`
       * :c:data:`pg_elem_mul`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_has_i`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_col`
       * :c:data:`pg_shell_target`

 
.. c:var:: int pg_elem_det


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_rot_mat	(3,3,pg_n_elem)
        integer, allocatable	:: pg_perm	(nucl_num,pg_n_elem)
        double precision, allocatable	:: pg_elem_angle	(pg_n_elem)
        double precision, allocatable	:: pg_elem_axis	(3,pg_n_elem)
        integer, allocatable	:: pg_elem_det	(pg_n_elem)
        integer, allocatable	:: pg_elem_kind	(pg_n_elem)


    Rotation matrices (input frame), induced atom permutations, angle,
    axis (plane normal for a reflection), determinant sign, and the
    integer (det,m,k) key of every group element.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_elem_inv`
       * :c:data:`pg_elem_mul`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_has_i`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_col`
       * :c:data:`pg_shell_target`

 
.. c:var:: int pg_elem_inv


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_elem_inv	(pg_n_elem)


    Index of the inverse of every element: the inverse atom permutation with
    the same determinant. Exact, no tolerance: ``h`` is the inverse of ``g`` iff
    its permutation is the inverse one and ``det`` agrees, and the enumeration
    stores permutations, never re-matches matrices.
    
    The column matching needs it to close a column under inversion: a
    geometric class is not closed under inversion (in C3 the class {C3}
    has its inverse in the class {C3^2}), while a file column always is.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_class`
       * :c:data:`pg_n_col`

 
.. c:var:: int pg_elem_kind


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_rot_mat	(3,3,pg_n_elem)
        integer, allocatable	:: pg_perm	(nucl_num,pg_n_elem)
        double precision, allocatable	:: pg_elem_angle	(pg_n_elem)
        double precision, allocatable	:: pg_elem_axis	(3,pg_n_elem)
        integer, allocatable	:: pg_elem_det	(pg_n_elem)
        integer, allocatable	:: pg_elem_kind	(pg_n_elem)


    Rotation matrices (input frame), induced atom permutations, angle,
    axis (plane normal for a reflection), determinant sign, and the
    integer (det,m,k) key of every group element.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_elem_inv`
       * :c:data:`pg_elem_mul`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_has_i`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_col`
       * :c:data:`pg_shell_target`

 
.. c:var:: int pg_elem_mul


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer, allocatable	:: pg_elem_mul	(pg_n_elem,pg_n_elem)


    Multiplication table: ``pg_elem_mul(g,h)`` is the index of R_g R_h.
    
    Built from the atom permutations and the determinants, never from the
    matrices. A stored element is only a *representative* of a symmetry
    operation known to :c:data:`pg_tol`: a candidate direction read off the geometry
    is accepted when the atoms move by less than :c:data:`pg_tol`, which for a
    rotation about a slightly tilted axis means a matrix accurate to
    :c:data:`pg_tol` / lever arm ~ 1e-8, and a product of two such matrices differs
    from the stored representative by the same order. A matrix comparison
    with a 1e-8 threshold is then a coin flip (measured on NH3: the C3 axis
    is tilted by 4.9e-9, 18 of the 36 products miss). The (permutation,
    determinant) key is exact and is the identity criterion the enumeration
    itself uses.
    
    The table is :math:`|G|^2 integers: 64 kB at the worst case :math:`|G|` = 128. Today the
    only consumer is gate G2; it is also what a product of two irreps will
    need once the module is in the core.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`


 
.. c:var:: double pg_geom_residual


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        double precision	:: pg_geom_residual	


    max over the elements g and the atoms i of
    ``R_g (r_i - center) + center - r_{pg_perm(i,g)}``, in bohr.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`


 
.. c:var:: char pg_group_name


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        character*(16)	:: pg_group_name	


    Name of the point group from the measured invariants.
    Covers all 204 finite groups of the data set; gives up instead of
    guessing on inconsistent invariants.
    
    Verified against all 209 finite files with the closed-form construction
    of the reference oracle (refsym2.py + ident.py): the tree returns the
    file name for 204 of them. The other five (C1h, C1v, D1, D1d, D1h) are
    degenerate aliases of Cs, Cs, C2, C2h and C2v: they are the n = 1
    cases where the standard naming rules collapse two names onto one
    group, and they cannot be produced by the tree. Their files are never
    opened then; do not "fix" the tree to produce them.
    
    The two facts that make the tree short are
     * n_axes_max counts the axes carrying a rotation of order n_max, so a
       D<n> group has n_axes_max = 1 for n >= 3 (its perpendicular axes
       have order 2) and cannot be mistaken for a cubic group;
     * in the n_max = 1 branch the orientation of the single mirror is not
       usable: z is the molecular plane normal for a planar molecule but an
       inertia axis otherwise, so C1/Ci/Cs are decided by :math:`|G|` and by
       "there is exactly one reflection".

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_has_i`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_col`

 
.. c:var:: bool pg_has_i


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        logical	:: pg_has_i	
        integer	:: pg_n_sigma_h	
        integer	:: pg_n_sigma_vtot	
        integer	:: pg_n_c2_perp	


    Measured invariants used by the identification.
    
    pg_n_sigma_vtot counts *all* vertical reflections, i.e. the planes whose
    normal is not along z. It equals n_max in both D_nh and D_nd, which is
    all the tree needs (n_sigma_h is tested first, so the two never
    conflict). There is no sigma_v/sigma_d split here: the split is a
    labelling device, it is decided from the plane normal in columns.irp.f, and
    it cannot be derived reliably from the geometry alone (the prototype's
    bisector test was inverted and still "worked" because both counts were
    only ever summed).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_group_name`

 
.. c:var:: double pg_inertia_eval


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_inertia_eval	(3)
        double precision, allocatable	:: pg_inertia_evec	(3,3)


    Principal moments of inertia (ascending) and eigenvectors (columns).
    Eigenvector signs are fixed deterministically (largest component
    positive) so that repeated runs are byte-identical.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`center_of_mass`
       * :c:data:`element_name`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_detection_failed`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_is_linear`
       * :c:data:`pg_is_planar`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_max`
       * :c:data:`pg_plane_normal`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`

 
.. c:var:: double pg_inertia_evec


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_inertia_eval	(3)
        double precision, allocatable	:: pg_inertia_evec	(3,3)


    Principal moments of inertia (ascending) and eigenvectors (columns).
    Eigenvector signs are fixed deterministically (largest component
    positive) so that repeated runs are byte-identical.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`center_of_mass`
       * :c:data:`element_name`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_detection_failed`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_is_linear`
       * :c:data:`pg_is_planar`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_max`
       * :c:data:`pg_plane_normal`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`

 
.. c:var:: bool pg_irrep_allowed


    File : :file:`symmetry/settings.irp.f`

    .. code:: fortran

        logical, allocatable	:: pg_irrep_allowed	(pg_n_irrep)


    ``pg_irrep_allowed(a)`` is true when irrep ``a`` is in :c:data:`irrep_states`.
    The whole list is allowed when :c:data:`irrep_states` is empty. Names are
    compared case-insensitively; an unknown name prints the available
    names and stops: malformed user input is the one case where stopping
    is allowed.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`irrep_states`
       * :c:data:`mpi_master`
       * :c:data:`pg_n_col`
       * :c:data:`pg_n_irrep_allowed`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`

 
.. c:var:: int pg_irrep_dim


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        integer	:: pg_n_col	
        integer	:: pg_n_irrep	
        character*(32), allocatable	:: pg_irrep_name	(pg_n_elem)
        double precision, allocatable	:: pg_char	(pg_n_elem,pg_n_elem)
        character*(8), allocatable	:: pg_data_sym	(pg_n_elem)
        integer, allocatable	:: pg_col_start	(pg_n_elem)
        integer, allocatable	:: pg_col_len	(pg_n_elem)
        integer, allocatable	:: pg_col_member	(pg_n_elem)
        integer, allocatable	:: pg_irrep_dim	(pg_n_elem)
        integer, allocatable	:: pg_irrep_mult	(pg_n_elem)


    Character table of :c:data:`pg_group_name`, with every file column
    matched to the set of group elements it stands for, and validated.
    
    Main outputs, for mo_irrep.irp.f:
    
    ``pg_irrep_name(a)``
        name of irrep a, as in the file (A1, B2u, T1g, E')
    
    ``pg_irrep_dim(a)``
        chi_a(E) = nint(pg_char(a,1))
    
    ``pg_irrep_mult(a)``
        1 for a genuine real irrep, 2 when the row is the sum of a
        complex-conjugate pair (C3, C10, ...). The MO projector is
        then (d_a/mult_a) sum_col chi_a(col) U_col.
    
    ``pg_char(a,c)``
        characters. Column c is file column c: the matching does not permute
        anything, so this is the parsed matrix verbatim.
    
    ``pg_col_len(c)``
        number of ELEMENTS of column c. This replaces the prototype's
        pg_class_size everywhere.
    
    ``pg_col_member(...)``
        the element indices, columns concatenated; column c is
        pg_col_member(pg_col_start(c) : ...).
    
    Validation: the columns partition the group, the first column is the
    identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
    orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
    sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
    :c:data:`pg_detection_failed`: the caller of this provider must give up.
    
    These identities do NOT detect a swap of two columns that have equal
    sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
    intact): the authoritative gate for the column *labelling* is the
    projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: bool pg_irrep_filter


    File : :file:`symmetry/settings.irp.f`

    .. code:: fortran

        logical	:: pg_irrep_filter	


    True when :c:data:`irrep_states` is non-empty and the symmetry machinery can
    honour it. Off when :c:data:`auto_symmetry` is off (the master switch) or when
    the point-group detection failed. It deliberately does not name
    :c:func:`pg_symmetry_ready`, which lives in symmetry_mo: that would make module
    symmetry depend on symmetry_mo. The two tests are written out here.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`auto_symmetry`
       * :c:data:`mpi_master`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_n_irrep_allowed`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`

 
.. c:var:: int pg_irrep_mult


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        integer	:: pg_n_col	
        integer	:: pg_n_irrep	
        character*(32), allocatable	:: pg_irrep_name	(pg_n_elem)
        double precision, allocatable	:: pg_char	(pg_n_elem,pg_n_elem)
        character*(8), allocatable	:: pg_data_sym	(pg_n_elem)
        integer, allocatable	:: pg_col_start	(pg_n_elem)
        integer, allocatable	:: pg_col_len	(pg_n_elem)
        integer, allocatable	:: pg_col_member	(pg_n_elem)
        integer, allocatable	:: pg_irrep_dim	(pg_n_elem)
        integer, allocatable	:: pg_irrep_mult	(pg_n_elem)


    Character table of :c:data:`pg_group_name`, with every file column
    matched to the set of group elements it stands for, and validated.
    
    Main outputs, for mo_irrep.irp.f:
    
    ``pg_irrep_name(a)``
        name of irrep a, as in the file (A1, B2u, T1g, E')
    
    ``pg_irrep_dim(a)``
        chi_a(E) = nint(pg_char(a,1))
    
    ``pg_irrep_mult(a)``
        1 for a genuine real irrep, 2 when the row is the sum of a
        complex-conjugate pair (C3, C10, ...). The MO projector is
        then (d_a/mult_a) sum_col chi_a(col) U_col.
    
    ``pg_char(a,c)``
        characters. Column c is file column c: the matching does not permute
        anything, so this is the parsed matrix verbatim.
    
    ``pg_col_len(c)``
        number of ELEMENTS of column c. This replaces the prototype's
        pg_class_size everywhere.
    
    ``pg_col_member(...)``
        the element indices, columns concatenated; column c is
        pg_col_member(pg_col_start(c) : ...).
    
    Validation: the columns partition the group, the first column is the
    identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
    orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
    sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
    :c:data:`pg_detection_failed`: the caller of this provider must give up.
    
    These identities do NOT detect a swap of two columns that have equal
    sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
    intact): the authoritative gate for the column *labelling* is the
    projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: char pg_irrep_name


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        integer	:: pg_n_col	
        integer	:: pg_n_irrep	
        character*(32), allocatable	:: pg_irrep_name	(pg_n_elem)
        double precision, allocatable	:: pg_char	(pg_n_elem,pg_n_elem)
        character*(8), allocatable	:: pg_data_sym	(pg_n_elem)
        integer, allocatable	:: pg_col_start	(pg_n_elem)
        integer, allocatable	:: pg_col_len	(pg_n_elem)
        integer, allocatable	:: pg_col_member	(pg_n_elem)
        integer, allocatable	:: pg_irrep_dim	(pg_n_elem)
        integer, allocatable	:: pg_irrep_mult	(pg_n_elem)


    Character table of :c:data:`pg_group_name`, with every file column
    matched to the set of group elements it stands for, and validated.
    
    Main outputs, for mo_irrep.irp.f:
    
    ``pg_irrep_name(a)``
        name of irrep a, as in the file (A1, B2u, T1g, E')
    
    ``pg_irrep_dim(a)``
        chi_a(E) = nint(pg_char(a,1))
    
    ``pg_irrep_mult(a)``
        1 for a genuine real irrep, 2 when the row is the sum of a
        complex-conjugate pair (C3, C10, ...). The MO projector is
        then (d_a/mult_a) sum_col chi_a(col) U_col.
    
    ``pg_char(a,c)``
        characters. Column c is file column c: the matching does not permute
        anything, so this is the parsed matrix verbatim.
    
    ``pg_col_len(c)``
        number of ELEMENTS of column c. This replaces the prototype's
        pg_class_size everywhere.
    
    ``pg_col_member(...)``
        the element indices, columns concatenated; column c is
        pg_col_member(pg_col_start(c) : ...).
    
    Validation: the columns partition the group, the first column is the
    identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
    orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
    sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
    :c:data:`pg_detection_failed`: the caller of this provider must give up.
    
    These identities do NOT detect a swap of two columns that have equal
    sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
    intact): the authoritative gate for the column *labelling* is the
    projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: bool pg_is_linear


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        logical	:: pg_is_linear	


    .True. for a collinear molecule. Two tests, both needed: the moment
    about the molecular axis is the *smallest* inertia eigenvalue and it is
    exactly zero for a collinear geometry (a maxval test can never fire),
    and every nucleus must sit within half pg_tol of the line through
    pg_center along pg_inertia_evec(:,1).
    
    Half pg_tol, not one pg_tol: the group of the molecule is built from
    matrices that are exact in the canonical frame, so the only error is the
    distance d of a nucleus from the axis, and a sigma_v moves it by 2d. A
    nucleus at d = pg_tol would make every sigma_v fail pg_maps_to_self, and
    the run would end on "element does not realise its permutation" instead
    of being treated as linear.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_elem_inv`
       * :c:data:`pg_elem_mul`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_max`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`

 
.. c:var:: bool pg_is_planar


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        logical	:: pg_is_planar	


    .True. when every nucleus lies within pg_tol of a plane. The candidate
    normal is the inertia eigenvector with the largest eigenvalue (for a
    planar molecule that eigenvalue belongs to the normal).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_rot`

 
.. c:var:: int pg_n_axes_max


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer	:: pg_n_max	
        integer	:: pg_n_axes_max	


    Highest proper-rotation order, and number of distinct lines carrying an
    axis of that order.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_group_name`

 
.. c:var:: int pg_n_c2_perp


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        logical	:: pg_has_i	
        integer	:: pg_n_sigma_h	
        integer	:: pg_n_sigma_vtot	
        integer	:: pg_n_c2_perp	


    Measured invariants used by the identification.
    
    pg_n_sigma_vtot counts *all* vertical reflections, i.e. the planes whose
    normal is not along z. It equals n_max in both D_nh and D_nd, which is
    all the tree needs (n_sigma_h is tested first, so the two never
    conflict). There is no sigma_v/sigma_d split here: the split is a
    labelling device, it is decided from the plane normal in columns.irp.f, and
    it cannot be derived reliably from the geometry alone (the prototype's
    bisector test was inverted and still "worked" because both counts were
    only ever summed).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_group_name`

 
.. c:var:: int pg_n_class


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer	:: pg_n_class	
        integer, allocatable	:: pg_class	(pg_n_elem)


    Classes by conjugation. Element 1 is the identity, so class 1 is the
    identity class. Direct orbit enumeration (:math:`|G|` <= 120).
    
    The conjugation is done on the atom permutations, not on the
    matrices: the stored matrices carry the :c:data:`pg_tol`-sized tilt of the
    candidate directions, so a conjugated matrix is not within 1e-8 of any
    stored one. The permutation (with the determinant, to separate the
    identity from the reflection in the molecular plane) is the exact
    invariant used by :c:func:`pg_add_element` too. Direct orbit enumeration,
    :math:`|G|` <= 128 so :math:`|G|^2` = 16384 comparisons.
    
    Linear molecules are the exception: all their permutations are the
    identity, so the conjugation is done on the analytic matrices there
    (:c:func:`pg_find_by_mat`).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_col`

 
.. c:var:: int pg_n_col


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        integer	:: pg_n_col	
        integer	:: pg_n_irrep	
        character*(32), allocatable	:: pg_irrep_name	(pg_n_elem)
        double precision, allocatable	:: pg_char	(pg_n_elem,pg_n_elem)
        character*(8), allocatable	:: pg_data_sym	(pg_n_elem)
        integer, allocatable	:: pg_col_start	(pg_n_elem)
        integer, allocatable	:: pg_col_len	(pg_n_elem)
        integer, allocatable	:: pg_col_member	(pg_n_elem)
        integer, allocatable	:: pg_irrep_dim	(pg_n_elem)
        integer, allocatable	:: pg_irrep_mult	(pg_n_elem)


    Character table of :c:data:`pg_group_name`, with every file column
    matched to the set of group elements it stands for, and validated.
    
    Main outputs, for mo_irrep.irp.f:
    
    ``pg_irrep_name(a)``
        name of irrep a, as in the file (A1, B2u, T1g, E')
    
    ``pg_irrep_dim(a)``
        chi_a(E) = nint(pg_char(a,1))
    
    ``pg_irrep_mult(a)``
        1 for a genuine real irrep, 2 when the row is the sum of a
        complex-conjugate pair (C3, C10, ...). The MO projector is
        then (d_a/mult_a) sum_col chi_a(col) U_col.
    
    ``pg_char(a,c)``
        characters. Column c is file column c: the matching does not permute
        anything, so this is the parsed matrix verbatim.
    
    ``pg_col_len(c)``
        number of ELEMENTS of column c. This replaces the prototype's
        pg_class_size everywhere.
    
    ``pg_col_member(...)``
        the element indices, columns concatenated; column c is
        pg_col_member(pg_col_start(c) : ...).
    
    Validation: the columns partition the group, the first column is the
    identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
    orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
    sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
    :c:data:`pg_detection_failed`: the caller of this provider must give up.
    
    These identities do NOT detect a swap of two columns that have equal
    sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
    intact): the authoritative gate for the column *labelling* is the
    projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_n_elem


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer	:: pg_n_elem	


    Order of the molecular point group.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_rot`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_elem_mul`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_group_name`
       * :c:data:`pg_has_i`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_col`
       * :c:data:`pg_rot_mat`
       * :c:data:`pg_shell_target`
       * :c:data:`pg_xg`

 
.. c:var:: int pg_n_irrep


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        integer	:: pg_n_col	
        integer	:: pg_n_irrep	
        character*(32), allocatable	:: pg_irrep_name	(pg_n_elem)
        double precision, allocatable	:: pg_char	(pg_n_elem,pg_n_elem)
        character*(8), allocatable	:: pg_data_sym	(pg_n_elem)
        integer, allocatable	:: pg_col_start	(pg_n_elem)
        integer, allocatable	:: pg_col_len	(pg_n_elem)
        integer, allocatable	:: pg_col_member	(pg_n_elem)
        integer, allocatable	:: pg_irrep_dim	(pg_n_elem)
        integer, allocatable	:: pg_irrep_mult	(pg_n_elem)


    Character table of :c:data:`pg_group_name`, with every file column
    matched to the set of group elements it stands for, and validated.
    
    Main outputs, for mo_irrep.irp.f:
    
    ``pg_irrep_name(a)``
        name of irrep a, as in the file (A1, B2u, T1g, E')
    
    ``pg_irrep_dim(a)``
        chi_a(E) = nint(pg_char(a,1))
    
    ``pg_irrep_mult(a)``
        1 for a genuine real irrep, 2 when the row is the sum of a
        complex-conjugate pair (C3, C10, ...). The MO projector is
        then (d_a/mult_a) sum_col chi_a(col) U_col.
    
    ``pg_char(a,c)``
        characters. Column c is file column c: the matching does not permute
        anything, so this is the parsed matrix verbatim.
    
    ``pg_col_len(c)``
        number of ELEMENTS of column c. This replaces the prototype's
        pg_class_size everywhere.
    
    ``pg_col_member(...)``
        the element indices, columns concatenated; column c is
        pg_col_member(pg_col_start(c) : ...).
    
    Validation: the columns partition the group, the first column is the
    identity, chi_a(E) is a positive integer, mult_a is 1 or 2,
    orthogonality sum_c ``col`` chi_a chi_b = :math:`|G|` mult_a delta_ab, and
    sum_a d_a^2/mult_a = :math:`|G|`. A failed check prints its reason once and sets
    :c:data:`pg_detection_failed`: the caller of this provider must give up.
    
    These identities do NOT detect a swap of two columns that have equal
    sizes in an Abelian group (swapping sv and sd in C2v keeps every sum
    intact): the authoritative gate for the column *labelling* is the
    projector trace test G3 of mo_irrep.irp.f, which uses the actual matrices.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_elem_inv`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_elem_col`
       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: int pg_n_irrep_allowed


    File : :file:`symmetry/settings.irp.f`

    .. code:: fortran

        integer	:: pg_n_irrep_allowed	


    Number of tokens in :c:data:`irrep_states` (separators: comma, semicolon,
    blank, tab). Zero when the list is empty. No group knowledge, so this
    provider does not build the detection chain.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`irrep_states`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_irrep_allowed`
       * :c:data:`pg_irrep_filter`

 
.. c:var:: int pg_n_max


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer	:: pg_n_max	
        integer	:: pg_n_axes_max	


    Highest proper-rotation order, and number of distinct lines carrying an
    axis of that order.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_group_name`

 
.. c:var:: int pg_n_sigma_h


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        logical	:: pg_has_i	
        integer	:: pg_n_sigma_h	
        integer	:: pg_n_sigma_vtot	
        integer	:: pg_n_c2_perp	


    Measured invariants used by the identification.
    
    pg_n_sigma_vtot counts *all* vertical reflections, i.e. the planes whose
    normal is not along z. It equals n_max in both D_nh and D_nd, which is
    all the tree needs (n_sigma_h is tested first, so the two never
    conflict). There is no sigma_v/sigma_d split here: the split is a
    labelling device, it is decided from the plane normal in columns.irp.f, and
    it cannot be derived reliably from the geometry alone (the prototype's
    bisector test was inverted and still "worked" because both counts were
    only ever summed).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_group_name`

 
.. c:var:: int pg_n_sigma_vtot


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        logical	:: pg_has_i	
        integer	:: pg_n_sigma_h	
        integer	:: pg_n_sigma_vtot	
        integer	:: pg_n_c2_perp	


    Measured invariants used by the identification.
    
    pg_n_sigma_vtot counts *all* vertical reflections, i.e. the planes whose
    normal is not along z. It equals n_max in both D_nh and D_nd, which is
    all the tree needs (n_sigma_h is tested first, so the two never
    conflict). There is no sigma_v/sigma_d split here: the split is a
    labelling device, it is decided from the plane normal in columns.irp.f, and
    it cannot be derived reliably from the geometry alone (the prototype's
    bisector test was inverted and still "worked" because both counts were
    only ever summed).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_group_name`

 
.. c:var:: int pg_perm


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_rot_mat	(3,3,pg_n_elem)
        integer, allocatable	:: pg_perm	(nucl_num,pg_n_elem)
        double precision, allocatable	:: pg_elem_angle	(pg_n_elem)
        double precision, allocatable	:: pg_elem_axis	(3,pg_n_elem)
        integer, allocatable	:: pg_elem_det	(pg_n_elem)
        integer, allocatable	:: pg_elem_kind	(pg_n_elem)


    Rotation matrices (input frame), induced atom permutations, angle,
    axis (plane normal for a reflection), determinant sign, and the
    integer (det,m,k) key of every group element.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_elem_inv`
       * :c:data:`pg_elem_mul`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_has_i`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_col`
       * :c:data:`pg_shell_target`

 
.. c:var:: double pg_plane_normal


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_plane_normal	(3)


    Unit normal of the molecular plane, in the INPUT frame.
    Valid when pg_is_planar is true; defined for every molecule.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_inertia_eval`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_rot`

 
.. c:var:: double pg_rot


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_rot	(3,3)


    Rows are the canonical x, y, z unit vectors, in the INPUT frame.
    Standard coordinates of a point r are pg_rot . (r - pg_center).
    
    Orientation follows the steps below, including the tie-breaks. Cubic and
    icosahedral labels are invariant under relabelling of their equivalent
    axes, so no attempt is made to impose a "standard cube orientation".

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_is_planar`
       * :c:data:`pg_plane_normal`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_has_i`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

 
.. c:var:: double pg_rot_mat


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pg_rot_mat	(3,3,pg_n_elem)
        integer, allocatable	:: pg_perm	(nucl_num,pg_n_elem)
        double precision, allocatable	:: pg_elem_angle	(pg_n_elem)
        double precision, allocatable	:: pg_elem_axis	(3,pg_n_elem)
        integer, allocatable	:: pg_elem_det	(pg_n_elem)
        integer, allocatable	:: pg_elem_kind	(pg_n_elem)


    Rotation matrices (input frame), induced atom permutations, angle,
    axis (plane normal for a reflection), determinant sign, and the
    integer (det,m,k) key of every group element.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pg_center`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot`
       * :c:data:`pg_tol`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_elem_inv`
       * :c:data:`pg_elem_mul`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_has_i`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_col`
       * :c:data:`pg_shell_target`

 
.. c:var:: double pg_thresh


    File : :file:`symmetry/settings.irp.f`

    .. code:: fortran

        double precision	:: pg_thresh	


    Threshold on representation-matrix elements used to decide whether an MO
    belongs to a symmetry block and whether a determinant has a definite
    irrep. Default 1e-6 (a weight).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`pg_col_diag`
       * :c:data:`pg_det_irrep`
       * :c:data:`pg_mo_weight`
       * :c:data:`pg_state_irrep_weight`

 
.. c:var:: double pg_tol


    File : :file:`symmetry/settings.irp.f`

    .. code:: fortran

        double precision	:: pg_tol	


    Geometries are considered symmetric when every atom maps onto an atom of
    the same element within this distance (bohr). Default 1e-4 bohr.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_is_linear`
       * :c:data:`pg_is_planar`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_max`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void pg_add_direction(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_add_direction(v, nmax, list, n, ierr)


    Append a direction to a deduplicated (up to sign) list.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_build_candidates`

 
.. c:function:: void pg_add_element(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_add_element(r, p, nat, nmax, rot, perm, n_rot, ierr)


    Add (r,p) to the list. Identical matrices to 1e-8 are duplicates;
    beyond that, two elements are the same group element when they induce
    the same atom permutation *and* have the same determinant. If two
    orthogonal maps permute the atoms identically, their quotient fixes
    every atom, hence is the identity unless the molecule is planar, in
    which case it is the reflection in the molecular plane (the only
    orthogonal map fixing a whole plane pointwise), which differs by its
    determinant. Linear molecules, the other degenerate case, are
    rejected before enumeration.
    
    The permutation test is required, not an optimisation: a candidate
    direction read off a geometry that is symmetric only to :c:data:`pg_tol` is
    tilted by O(:c:data:`pg_tol` / lever arm), so a product of two accepted
    elements differs from the stored element by far more than 1e-8 and
    would be appended again, without bound, until the array overflows.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_enumerate`

 
.. c:function:: void pg_axis_along(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        logical function pg_axis_along(axis, k)


    True when the unit vector ``axis`` (INPUT frame) lies on the ``k``-th
    canonical axis (``k`` = 1, 2, 3 for x, y, z). The canonical direction in
    the input frame is row ``k`` of pg_rot, so this is a dot
    product, exactly as the identification invariants do it with row 3.
    
    This is a LABELLING test and it is deliberately tight (1e-4, 0.8
    degree): it must answer "is this element one of the cartesian axes",
    not "is it roughly aligned". The closest non-axis direction in the data
    set is the sigma_d normal of D32d at pi/(2*32) = 2.8 degrees from x,
    ``cos`` = 0.9988, while an element that IS an axis has ``cos`` = 1 up to the
    tilt between two candidate directions built from the same geometry: for
    a geometry symmetric to 1e-12 that tilt is 1e-9, and it stays below 1e-7
    for any :c:data:`pg_tol` up to 1e-2.
    
    The identification invariants use 1/2 instead: they must survive a
    geometry that is symmetric only within :c:data:`pg_tol`.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_rot`

 
.. c:function:: void pg_build_candidates(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_build_candidates(nat, xyz, znum, ctr, evecs, ncandmax, &
             cand, ncand, ierr)


    Build the list of candidate symmetry axes / plane normals.
    ``ierr`` = 1 when the list overflows ``ncandmax``.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_enumerate`
       * :c:func:`pg_proper_scan`
       * :c:func:`pg_select_x`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_add_direction`
       * :c:func:`pg_cross`

 
.. c:function:: void pg_build_group(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_build_group(nat, xyz, znum, ctr, evecs, frame, linear, tol, &
             nmax, rot, perm, n_rot, ierr)


    Build the group of the molecule and store it as matrices in the INPUT
    frame plus the induced atom permutations.
    
    ``frame`` is :c:data:`pg_rot` (rows = canonical axes in the input frame),
    used by the linear branch only. ``evecs`` is :c:data:`pg_inertia_evec`,
    used by the generic branch only. This dispatcher exists because IRPF90
    does not look inside a subroutine body: the two callers
    (:c:data:`pg_n_elem`, :c:data:`pg_rot_mat`) name :c:data:`pg_is_linear` and
    :c:data:`pg_rot` themselves so that the providers are ordered.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_build_linear`
       * :c:func:`pg_enumerate`

 
.. c:function:: void pg_build_linear(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_build_linear(nat, xyz, znum, ctr, frame, tol, nmax, &
             rot, perm, n_rot, ierr)


    The finite group used for a linear molecule is C31v, extended to
    C31v x {E, sigma_h} when sigma_h maps the molecule onto itself. The
    order is 62 or 124.
    
    Why not enumerate: for a collinear geometry every rotation about the
    axis and every reflection in a plane containing it induces the identity
    permutation, so the (permutation, determinant) key of
    :c:func:`pg_add_element` collapses the 31 sigma_v planes with sigma_h and
    the enumeration keeps only {E, one sigma_v}. The matrices are therefore
    built analytically from the canonical frame of :c:data:`pg_rot`: every
    angle is an exact multiple of 2 pi / 31, so the inverse, product and
    conjugation lookups by matrix (1e-8) are stable.
    
    C31v, not Cinfv, and C31v x {E,sigma_h}, not C31v x Ci: the tables
    shipped are C31v and D31h. With the inversion the improper elements
    would be i C31^k = S62^{31+2k}, k = 1..15: fifteen columns named
    S62^33 .. S62^61, a notation no other table uses, and the inversion
    itself would be S62^31, an unreduced symbol (gcd(31,62) = 31) that the
    column parser rejects. The horizontal plane keeps sigma_h and S31^k.
    The rows of the D31h table are the restrictions of the Dinfh irreps to
    this group, so Sigma_g+, Pi_u, ... keep their chemical meaning.
    
    ``ierr`` = 1 when the group does not fit in ``nmax``, 2 when a built
    element fails :c:func:`pg_maps_to_self`.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_build_group`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_linear_add`
       * :c:func:`pg_linear_to_input`
       * :c:func:`pg_maps_to_self`

 
.. c:function:: void pg_character_table_path(void)


    File : :file:`symmetry/chartab.irp.f`

    .. code:: fortran

        subroutine pg_character_table_path(group, path, ierr)


    Build the full path ``$QP_ROOT/data/Symmetry/<group>`` and check that
    it exists. ``ierr`` = 1 when QP_ROOT is not set, 2 when the file does
    not exist. ``path`` is always the path that was tried.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_col`
       * :c:func:`pg_report`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_environment_variable`

 
.. c:function:: void pg_class_any_used(void)


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        logical function pg_class_any_used(cl, used)


    True when some member of class ``cl`` has already been taken by a column.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`

 
.. c:function:: void pg_class_has_axis(void)


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        logical function pg_class_has_axis(cl, k)


    True when some member of class ``cl`` has its axis along canonical axis ``k``
    (1 = x, 2 = y, 3 = z). For a reflection the "axis" is the plane normal,
    so this is exactly "the plane contains z" when ``k`` = 2 and "the plane is
    vertical" when ``k`` = 3.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

 
.. c:function:: void pg_class_has_key(void)


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        logical function pg_class_has_key(cl, key)


    True when some member of class ``cl`` has the encoded key ``key``.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

 
.. c:function:: void pg_class_kset_subset(void)


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        logical function pg_class_kset_subset(cl, ka, kb)


    True when every key of class ``cl`` is ``ka`` or ``kb`` (the key set is a subset
    of {``ka``,``kb``}). Membership of ``ka`` itself is a separate test.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

 
.. c:function:: void pg_col_fail(void)


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        subroutine pg_col_fail(reason, s, n, syms)


    Error of the column matcher. Prints the reason, the failing symbol, the
    symbol order of the file, and the class table with the key set and the
    axis flags of every class, then sets :c:data:`pg_detection_failed` and returns:
    the caller must give up at once.
    
    Only variables of group.irp.f are read here: everything of the
    pg_n_col provider block (pg_data_sym included) must be passed in, or
    IRPF90 would insert a build call for the provider that is currently
    being built, i.e. infinite recursion.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_class_size`
       * :c:data:`pg_detection_failed`
       * :c:data:`pg_group_name`
       * :c:data:`pg_n_class`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_col`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_decode_kind`

 
.. c:function:: void pg_count_tokens(void)


    File : :file:`symmetry/chartab.irp.f`

    .. code:: fortran

        subroutine pg_count_tokens(line, ntok)


    Count the blank- or TAB-separated fields of ``line``.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_parse_table`

 
.. c:function:: void pg_cross(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_cross(a, b, c)


    Cross product.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_build_candidates`
       * :c:func:`pg_elem_props`
       * :c:data:`pg_rot`

 
.. c:function:: void pg_decode_kind(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_decode_kind(kind, det, m, k)


    Inverse of pg_encode_kind.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_col_fail`
       * :c:data:`pg_has_i`
       * :c:func:`pg_report`

 
.. c:function:: void pg_det3(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        double precision function pg_det3(r)


    Determinant of a 3x3 matrix.

 
.. c:function:: void pg_diag3(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_diag3(a, e)


    dsyev wrapper for a 3x3 symmetric matrix (destroys ``a``).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_detection_failed`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_inertia_eval`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dsyev`

 
.. c:function:: void pg_elem_props(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_elem_props(r, angle, axis, det)


    Angle and axis: for an improper element
    the angle is theta_p - pi, where theta_p is the angle of -r. The
    oriented angle is obtained from atan2 of the antisymmetric part, not
    from acos, which is why S4 gives theta_p = 3pi/2 and theta = pi/2.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_rot_mat`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_cross`
       * :c:func:`pg_fix_sign`

 
.. c:function:: void pg_encode_kind(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_encode_kind(det, m, k, kind)


    Injective integer encoding of (det,m,k): kind = 2*(64*m + k) + (1-det)/2,
    with m <= 64 and k < m. For m <= 64 the pair (m,k) is recovered as
    k = mod(kind/2, 64), m = (kind/2 - k)/64.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_col`
       * :c:data:`pg_rot_mat`

 
.. c:function:: void pg_enumerate(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_enumerate(nat, xyz, znum, ctr, evecs, tol, nmax, &
             rot, perm, n_rot, ierr)


    Build the full group from the geometry. ``ierr`` = 1 when ``nmax`` is too
    small, 2 when the candidate list overflows, 3 on a non-orthogonal
    product (should never happen).

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_build_group`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_add_element`
       * :c:func:`pg_build_candidates`
       * :c:func:`pg_maps_to_self`
       * :c:func:`pg_reflection`
       * :c:func:`pg_rotation`

 
.. c:function:: void pg_expected_order(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_expected_order(name, order, ierr)


    ``C_n`` = n, ``C_nv`` = ``C_nh`` = 2n, ``S_n`` = n, ``D_n`` = 2n,
    ``D_nh`` = ``D_nd`` = 4n, :math:`|T|` = 12, ``Th`` = ``Td`` = :math:`|O|` = 24,
    ``Oh`` = 48, :math:`|I|` = 60, ``Ih`` = 120.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_group_name`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_read_int`

 
.. c:function:: void pg_find_by_mat(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer function pg_find_by_mat(r)


    Index of the stored element whose matrix equals ``r`` to 1e-8, or 0.
    
    Only for the linear branch, where the matrices are analytic (angles are
    exact multiples of 2 pi / 31) and a transpose or a product of two stored
    matrices is within 1e-15 of a stored one. The generic branch must keep
    :c:func:`pg_find_by_perm`: there the matrices carry the tilt of a candidate
    direction, and a 1e-8 matrix comparison is a coin flip.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

 
.. c:function:: void pg_find_by_perm(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        integer function pg_find_by_perm(p, det)


    Index of the unique element that induces the atom permutation ``p`` and has
    determinant ``det``, or 0. This is the enumeration's own identity
    criterion (:c:func:`pg_add_element`): two orthogonal maps that permute the atoms
    identically and share their determinant are the same element. No matrix
    is ever compared.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_rot_mat`

 
.. c:function:: void pg_fix_sign(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_fix_sign(v)


    Deterministic orientation of a vector.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_elem_props`
       * :c:data:`pg_inertia_eval`
       * :c:func:`pg_orient_axis`
       * :c:data:`pg_plane_normal`

 
.. c:function:: void pg_linear_add(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_linear_add(nat, xyz, znum, ctr, frame, rc, tol, nmax, &
             rot, perm, n_rot, ierr)


    Transform ``rc`` to the input frame, check that it maps the molecule
    onto itself, and append it. ``ierr`` = 1 when ``nmax`` is exhausted, 2
    when the element does not realise a permutation.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_build_linear`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_linear_to_input`
       * :c:func:`pg_maps_to_self`

 
.. c:function:: void pg_linear_to_input(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_linear_to_input(frame, rc, r)


    Rows of ``frame`` are the canonical axes expressed in the input frame,
    so ``s = frame (x - center)`` are the canonical coordinates and the
    input-frame matrix of a canonical operator ``R_c`` is
    ``R = frame^T R_c frame``.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_build_linear`
       * :c:func:`pg_linear_add`

 
.. c:function:: void pg_maps_to_self(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_maps_to_self(nat, xyz, znum, ctr, rot, tol, perm, ok)


    Test rot r_i == r_{perm(i)} for every atom.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_build_linear`
       * :c:func:`pg_enumerate`
       * :c:func:`pg_linear_add`
       * :c:func:`pg_proper_scan`
       * :c:func:`pg_select_x`

 
.. c:function:: void pg_normalize_irrep_states(void)


    File : :file:`symmetry/settings.irp.f`

    .. code:: fortran

        subroutine pg_normalize_irrep_states(list)


    Empty ``list`` when it holds only the token "None" (any case).

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`irrep_states`

 
.. c:function:: void pg_orient_axis(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_orient_axis(nat, xyz, ctr, tol, axis)


    Deterministic orientation of an axis.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_proper_scan`
       * :c:data:`pg_rot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_fix_sign`

 
.. c:function:: void pg_parse_symbol(void)


    File : :file:`symmetry/columns.irp.f`

    .. code:: fortran

        subroutine pg_parse_symbol(s, det, m, k, ierr)


    Split a generic column symbol. Accepted forms: C<n>, C<n>^<k>, S<n>,
    S<n>^<k>. Returns det = +1 for C, -1 for S, the order ``m`` (2..64) and the
    power ``k`` reduced modulo ``m``. ``ierr`` = 1 on any other form, 2 when
    gcd(k,m) /= 1 (an unreduced symbol such as C6^2: the data set never has
    one, and matching it would select the wrong classes).

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_col`

 
.. c:function:: void pg_parse_table(void)


    File : :file:`symmetry/chartab.irp.f`

    .. code:: fortran

        subroutine pg_parse_table(file, group_name, n_irrep, n_class, &
             irrep_name, sym, char, ierr)


    Parse one character-table file. Pure: no EZFIO keyword, no provider, no
    module-level state, nothing printed. Returns ``ierr`` /= 0 on any structural
    or I/O error.
    
    :file:       path of the file to parse
    :group_name: group name, taken from the file name (line 1 is ignored)
    :n_irrep:    number of irreps (rows) read from line 2
    :n_class:    number of table columns (the file format calls them
                 classes; it is == ``n_irrep`` in the data set). This is NOT
                 the geometric class count: C3 has two columns and
                 three classes.
    :irrep_name: irrep names, first ``n_irrep`` entries set
    :sym:        class symbols as found in the file
    :char:       characters, ``char(a,c)`` = chi_a(C_c)
    :ierr:       0 on success, > 0 otherwise

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_col`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_count_tokens`
       * :c:func:`pg_split_tab`

 
.. c:function:: void pg_proper_scan(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_proper_scan(nat, xyz, znum, ctr, evals, evecs, tol, &
             ncandmax, n_max, n_axes_max, axis, linear, ierr)


    Enumerate proper rotations of order 2..32 among the candidate
    directions, keep the highest order and select the principal axis.
    32, not 12: C32 occurs in D32h/D32d, the largest groups of the data set.
    A linear molecule is capped at 31: its improper elements are S_n^k and
    the shipped file names the 15 improper columns S31^k, k = 1..15. A
    composite order would not keep that scheme -- some powers are not S_n^k
    at all but reduce to a lower order (S25^5 = S5), so the columns would
    mix orders and "one column per power of S_n" would no longer hold.
    31 is the largest odd order <= 32, and 4*31 = 124 still fits the
    160-element arrays; the tables shipped are C31v and D31h.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_n_max`
       * :c:data:`pg_rot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_build_candidates`
       * :c:func:`pg_maps_to_self`
       * :c:func:`pg_orient_axis`
       * :c:func:`pg_rotation`

 
.. c:function:: void pg_read_int(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_read_int(name, i, n, value)


    Read a decimal integer starting at name(i); advance i past it.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_expected_order`

 
.. c:function:: void pg_reduce_angle(void)


    File : :file:`symmetry/group.irp.f`

    .. code:: fortran

        subroutine pg_reduce_angle(angle, m, k, ierr)


    Rational reduction of an angle: ``angle - 2 pi k/m`` < 1e-6.
    With this convention a reflection is (m=1,k=0) and the inversion is
    (m=2,k=1) (S1 = sigma, S2 = i).

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_rot_mat`

 
.. c:function:: void pg_reflection(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_reflection(n, rot)


    Householder reflection R = 1 - 2 n n^T.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_enumerate`

 
.. c:function:: void pg_rotation(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_rotation(axis, theta, rot)


    Proper rotation matrix, angle in radians.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_enumerate`
       * :c:func:`pg_proper_scan`
       * :c:func:`pg_select_x`

 
.. c:function:: void pg_select_x(void)


    File : :file:`symmetry/frame.irp.f`

    .. code:: fortran

        subroutine pg_select_x(nat, xyz, znum, ctr, evecs, tol, ncandmax, &
             zaxis, planar, pnormal, xaxis, ierr)


    Choose the canonical x from the branches below.

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`pg_rot`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`pg_build_candidates`
       * :c:func:`pg_maps_to_self`
       * :c:func:`pg_rotation`

 
.. c:function:: void pg_split_tab(void)


    File : :file:`symmetry/chartab.irp.f`

    .. code:: fortran

        subroutine pg_split_tab(line, nmax, nf, fld)


    Split ``line`` on TAB and store at most ``nmax`` left-padded fields in ``fld``.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`pg_parse_table`

