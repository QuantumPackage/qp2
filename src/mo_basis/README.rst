========
mo_basis
========

Molecular orbitals are expressed as

.. math::

  \phi_k({\bf r}) = \sum_i C_{ik} \chi_k({\bf r})


where :math:`\chi_k` are *normalized* atomic basis functions.

The current set of |MOs| has a label `mo_label`.
When the orbitals are modified, the label should also be updated to keep
everything consistent.

When saving the |MOs|, the :file:`mo_basis` directory of the |EZFIO| database
is copied in the :file:`save` directory, named by the current `mo_label`. All
this is done with the script named :file:`save_current_mos.sh` in the
:file:`$QP_ROOT/scripts` directory.



