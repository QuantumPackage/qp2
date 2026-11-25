===================
domain_perturbation
===================

Perturbation theory derived from the Dyson equation with domains in the space
of determinants.


The user point of view
----------------------

The :ref:`domain_pert` program uses the wave function stored in the |EZFIO|
directory as a reference, and computes the perturbative corrections to the
energy from order 2 to order :option:`domain_perturbation n_pert`, for all the
states.

* :option:`domain_perturbation n_pert` : maximum order.

* :option:`domain_perturbation n_det_domain` : maximum number of determinants in
  the domains of the external determinants. With 1, the domains are singletons
  and the second order is the Epstein-Nesbet |PT2|.

* :option:`domain_perturbation n_walkers` : average number of walkers per
  generation. With 0, all the sums are done exactly.  The cost then grows with
  the size of the space which is reached: it should be used only for small
  systems, to check the stochastic results.

* :option:`domain_perturbation n_replicas` : number of independent replicas, from
  which the statistical errors are computed.

* :option:`domain_perturbation exact_coupling` : if |true|, the coupling between
  the reference and the external determinants is stored. If |false|, it is
  sampled.

* :option:`domain_perturbation domain_seed` : seed of the random numbers. The
  random numbers do not depend on the number of threads.

* :option:`domain_perturbation domain_chebyshev` : if |true|, the series is
  evaluated at the energy :math:`E_g = E_n + E_2` and accelerated with Chebyshev
  polynomials. The plain Brillouin-Wigner series is also printed.
  If |false|, the Rayleigh-Schrodinger series is computed.

* :option:`domain_perturbation domain_cheb_auto`,
  :option:`domain_perturbation domain_cheb_min`,
  :option:`domain_perturbation domain_cheb_max` : interval of the Chebyshev
  polynomials.

Example::

  qp set determinants n_det_max 1000
  qp run fci
  qp set domain_perturbation n_pert 4
  qp set domain_perturbation n_det_domain 20
  qp run domain_pert

The total energy and its error are stored in
:option:`domain_perturbation domain_energy` and
:option:`domain_perturbation domain_energy_error`.


Theory
------

Let :math:`\mathcal{D}_0` be the set of the determinants of the wave function, with
projector :math:`P`, :math:`Q=1-P`, and :math:`\Phi_n`, :math:`E_n` the eigenstates of :math:`PHP`. The
domain of all the determinants of :math:`\mathcal{D}_0` is :math:`\mathcal{D}_0`. The domain
:math:`\mathcal{D}_j` of an external determinant :math:`j` contains :math:`j` and at most
``n_det_domain-1`` external determinants which are connected to :math:`j` by
:math:`H` and have a lower diagonal element, the lowest ones being kept.

For a state :math:`n`, with :math:`v = Q H \Phi_n` and
:math:`G_{jl}(E) = \langle j | (P_j H P_j - E)^{-1} | l \rangle`,

.. math::

   E = E_n + \sum_{p \ge 2} \Delta_p(E), \quad
   \Delta_{q+2}(E) = v^T K(E)^q s(E), \quad s = - G v,

.. math::

   K_{jk} = - \sum_{l \in \mathcal{D}_j} G_{jl} H_{lk}
   \; (j \in Q, k \notin \mathcal{D}_j), \quad
   K_{ik} = \sum_{l \in \mathcal{D}_0} R_{il} H_{lk}
   \; (i \in \mathcal{D}_0, k \in Q),

with the reduced resolvent
:math:`R(E) = \sum_{m \ne n} |\Phi_m\rangle (E-E_m)^{-1} \langle\Phi_m|`.
The order is the number of times that a path leaves a domain. The energies
are obtained by the expansion of :math:`\lambda^p \Delta_p(E(\lambda))` in powers of
:math:`\lambda` (Rayleigh-Schrodinger series).

All the quantities are Taylor series in :math:`E-E_n`, which are propagated with
the walkers. Each replica gives unbiased estimates of the Taylor
coefficients of :math:`\Delta_p`; the energies are combinations of their
averages, with a jackknife correction of the bias due to products of averages.

With Chebyshev polynomials, the sum of the terms :math:`t_q = v^T K^q s` is
replaced by :math:`\sum_q w_q t_q`. The second order is kept with a weight 1, and
the rest, :math:`v^T (1-K)^{-1} K s`, is approximated with
:math:`(1-\pi(K))/(1-K)`, where :math:`\pi` is the Chebyshev polynomial of an
interval which contains the spectrum of :math:`K(E_g)`, with
:math:`\pi(1)=1`. The energy is obtained self-consistently from
the Taylor expansion of the terms around :math:`E_g`. The interval can be
estimated with the Ritz values given by the moments :math:`t_q` (available if
:option:`domain_perturbation n_pert` is at least 5). This is valid for the
ground state with singletons, where :math:`-t_q` are the moments of a positive
measure; in the other cases it is only heuristic. For excited states the spectrum of
:math:`K` can be complex, and the Chebyshev series can then diverge while the
plain series converges. A warning is printed when the moments are
not those of a positive measure.

Limitations:

* the states are treated one by one: there is no coupling between the
  states, and the reduced resolvent :math:`R` contains the other eigenstates of
  :math:`PHP`. An external determinant which is almost degenerate with the state
  gives large corrections. The smallest :math:`|H_{jj}-E_n|` is printed.

* the external domains are not spin-adapted.


Validation
----------

The implementation was compared to the following data.

* Deterministic evaluation, reference |HF|, water molecule, 6-31G, all
  electrons, singletons, orders 2 to 10: agreement within
  6e-9 a.u. with an independent |PySCF| calculation of the Epstein-Nesbet
  series (E_2=-0.170886643, E_3=0.054297907, E_4=-0.030779358,
  E_5=0.016404078, E_6=-0.009232419, E_7=0.005316275, E_8=-0.003355395,
  E_9=0.002414737, E_10=-0.001984765). The differences probably come
  from the Bohr radius used by the two programs.

* Same molecule, frozen core, domains of at most 20 determinants,
  reference |HF|: E_2=-0.162630, E_3=0.026070, E_4=0.004267,
  which differ from a Python code by 7e-9, 4e-7 and 5e-7. The
  differences at orders 3 and 4 can come from determinants with the same
  diagonal element, for which the Python code and the module do
  not make the same choice.

* Water, STO-3G, frozen core, multideterminant reference (31 determinants),
  3 states, orders 2 to 7, domains of 1, 4 and 20 determinants: agreement
  within 1.e-8 a.u. with a dense matrix code which builds K from its
  definition. The same code checks that the Brillouin-Wigner series summed
  at the exact energy gives the exact eigenvalue.

* H_3 (doublet), cc-pVDZ, CIS reference (30 determinants), 3 states,
  orders 2 to 5, domains of 1 and 6 determinants: agreement with the same
  code to 12 digits, with the integrals of this program. These are the values of the
  automatic test.

* Second order with singletons, water, 6-31G, all electrons, CIPSI reference
  of 110 determinants: -0.040014667848, identical to the Epstein-Nesbet
  |PT2| of the :ref:`pt2` program (modified to use this formula). With 1089 and
  10833 determinants the differences are 1e-9 and 2e-10, because :ref:`pt2`
  uses only 1080 and 10678 selectors and generators.

* Stochastic evaluation, compared to the deterministic one: 150 seeds with
  20 replicas, for 3 states and orders 2 to 7, with exact_coupling=False. No
  bias is seen (largest deviation 1.9 standard deviations), and the standard
  deviation of the results is between 0.8 and 1.04 times the mean of
  the estimated error.

* Chebyshev series: agreement of the weights, of the plain series and of the
  accelerated series with the dense code (1e-8 for the first two states of
  water in STO-3G; 2e-5 or less for the third one, for which the
  Taylor expansion is limited). Water, 6-31G, frozen core, reference |HF|,
  order 12: -76.119955298 instead of -76.122465967 for the plain series
  (|FCI|: -76.119955188). All electrons, 1089 determinants, order 8:
  error 0.1 micro-hartree instead of 17. Third state of H_3 (CIS
  reference): the Chebyshev series diverges, with a warning.

The programmer's point of view
------------------------------

* ``det_table_module.f90`` : sharded hash table of determinants. The amplitudes
  are Taylor expansions.

* ``connected.irp.f`` : connected determinants, domains and local resolvents.

* ``reduced_resolvent.irp.f`` : reduced resolvent in the reference space.

* ``propagate.irp.f`` : one generation.

* ``domain_series.irp.f`` : replicas, and Rayleigh-Schrodinger series.

* ``chebyshev.irp.f`` : weights, interval and self-consistent energy.
