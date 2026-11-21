============
two_body_rdm
============

Two-electron density matrices (2-RDMs) in the |MO| basis.

The 2-RDMs of the :math:`\alpha\alpha`, :math:`\beta\beta` and
:math:`\alpha\beta` electron pairs are stored as arrays, with the
physicists' notation, consistent with the two-electron integrals in
the |MO| basis. The state specific matrices are
:c:data:`full_occ_2_rdm_aa_mo`, :c:data:`full_occ_2_rdm_ab_mo` and
:c:data:`full_occ_2_rdm_bb_mo`, with the spin traced variant
:c:data:`full_occ_2_rdm_spin_trace_mo`. The versions restricted to
the active orbitals are :c:data:`act_2_rdm_aa_mo`,
:c:data:`act_2_rdm_ab_mo`, :c:data:`act_2_rdm_bb_mo` and
:c:data:`act_2_rdm_spin_trace_mo`.

The `io_two_body_rdm_*` subroutines read and write the 2-RDMs in the
|EZFIO| database.

The program :ref:`test_2_rdm` is a small example of the use of these
providers.
