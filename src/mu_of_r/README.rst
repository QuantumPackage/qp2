=======
mu_of_r
=======

The on-top pair density :math:`\mu(r)` of the wave function.

For a |HF| wave function, :c:data:`mu_of_r_hf` is :math:`\mu(r)`, as
in Eq. (37) of :cite:`Giner_2018` for a |HF| reference. For a wave
function expanded in an active space, :c:data:`mu_of_r_psi_cas` is
the same quantity for a |CAS| wave function, and
:c:func:`give_mu_of_r_cas` returns :math:`\mu(r)`, :math:`f_\Psi`
and :math:`n_2^\Psi`.

The Fock space components of the two-body density matrix that enter
the expression of :math:`\mu(r)`, with the Cholesky decomposition of
the two-electron integrals, are provided by
:c:data:`f_hf_cholesky`, :c:data:`mos_times_cholesky_r1` and
:c:data:`mos_times_cholesky_r2`, and the :c:data:`f_psi_hf_ab` /
:c:data:`f_psi_cas_ab` variants.

The `read_mu_of_r` and `write_mu_of_r` routines store and restore
:math:`\mu(r)` in the |EZFIO| database, and the program
:ref:`test_proj_op` is a small example.
