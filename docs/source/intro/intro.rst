========
The |qp|
========

.. image:: /_static/qp2.png
   :align: center
   :width: 200px
   :alt: Quantum Package


What it is
==========

The |qp| is an open-source **programming environment** for quantum chemistry.
It is designed from the **developer** point of view, to help design new
quantum chemistry methods,
especially for `wave function theory <https://en.wikipedia.org/wiki/Ab_initio_quantum_chemistry_methods>`_ (|WFT|).

From the **user** point of view, the |qp| proposes a self-contained route
to optimized selected configuration interaction (|sCI|) based on the
|CIPSI| algorithm, which can efficiently reach near-full configuration
interaction (|FCI|) quality for relatively large systems.
The :doc:`quick-start guide <../users_guide/quickstart>` runs a complete
|CIPSI| calculation on a small molecule.

The main goal is the development of |sCI| methods and multi-reference
perturbation theory (|MRPT|) in the determinant-driven paradigm.
The code also contains the very basics of Kohn-Sham
`density functional theory <https://en.wikipedia.org/wiki/Density_functional_theory>`_
(|KS-DFT|) and of `range-separated hybrids <https://aip.scitation.org/doi/10.1063/1.1383587>`_ (|RSH|).

The determinant-driven framework allows the programmer to include any arbitrary
set of determinants in the variational space, which gives complete freedom in
methodological development. The basic ingredients of |RSH| together with those
of the |WFT| framework available in the |qp| library allow one to develop
range-separated DFT (|RSDFT|) approaches (see for instance the plugins of the
`qp_plugins_eginer <https://gitlab.com/eginer/qp_plugins_eginer>`_ repository).

All the programs are developed with the |IRPF90| code generator, which
considerably simplifies the collaborative development, and the addition of new
features.


What it is not
==============

The |qp| is *not* a general purpose quantum chemistry program.
First of all, it is a *library* to develop new theories and algorithms in
quantum chemistry. Beside the programs of the core modules, users are expected
to write their own programs.

The |qp| has been designed for |sCI|, so none of its algorithms is meant for
SCF or DFT calculations on thousands of atoms. Currently, the systems targeted
have less than 600 molecular orbitals. This limit is due to the memory
bottleneck induced by the storage of the two-electron integrals (see the
``mo_two_e_ints`` and ``ao_two_e_ints`` modules).

The |qp| is *not* a massive production code. For conventional
methods such as Hartree-Fock, CISD or MP2, the users are recommended to use the
existing standard production codes which are designed to make these methods run
fast. Again, the role of the |qp| is to make life simple for the
developer. Once a new method is developed and tested, the developer is encouraged
to consider re-expressing it with an integral-driven formulation, and to
implement the new method in open-source production codes, such as `NWChem`_
or |GAMESS|.


A few examples of applications
==============================

Multiple methods were developed with the |qp|:

- Selected Full-CI + Epstein-Nesbet PT2 (CIPSI), through the :ref:`fci` and
  :ref:`pt2` programs
  :cite:`Caffarel_2016,Caffarel_2016.2,Loos_2018,Scemama_2018,Dash_2018`
- Hybrid stochastic/deterministic MR-PT2, through the ``mrci`` program
  :cite:`Garniron_2017.2,Loos_2018`
- Orbital optimization for open-shell systems, through the ``orb_opt`` and
  ``optimization`` programs :cite:`Giner2016Mar,Giner_2017.3`
- CIS, CISD and MP2, through the :ref:`cis`, :ref:`cisd` and ``mp2`` programs
- Selected CISD and selected CAS-SD, obtained with the :ref:`fci` program by
  restricting the excitation degree (:option:`cipsi excitation_max`) or the
  reference determinants (:option:`cipsi excitation_ref`)
- Jeziorsky-Monkhorst MR-PT2 :cite:`Giner_2017`
- Effective Hamiltonian for variational MR wave functions :cite:`Giner_2017.2`
- Selected difference-dedicated CI (DD-CI)
- Multi-Reference Coupled Cluster (MR-CCSD), through the ``lccsd`` and ``ccsd``
  programs :cite:`Giner_2016,Garniron_2017`
- Shifted-Bk with CIPSI :cite:`Garniron_2018`
- CIPSI with range-separated DFT, through the plugins of the
  `qp_plugins_eginer <https://gitlab.com/eginer/qp_plugins_eginer>`_ repository
- DFT for basis set corrections :cite:`Giner_2018`
- Trial wave functions for |QMC|, with or without pseudo-potentials, through the
  ``save_for_qmcchem``, ``save_for_qmcpack`` and ``save_for_champ`` programs of
  the QMC plugins

Most of these programs can generate ground and excited states, and spin pure
wave functions (eigenstates of |S^2|) when :option:`determinants s2_eig` is
set to |true|.
