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
It has been built from the **developper** point of view in order to help
the design of new quantum chemistry methods,
especially for `wave function theory <https://en.wikipedia.org/wiki/Ab_initio_quantum_chemistry_methods>`_ (|WFT|).

From the **user** point of view, the |qp| proposes a stand-alone path
to use optimized selected configuration interaction |sCI| based on the
|CIPSI| algorithm that can efficiently reach near-full configuration interaction
|FCI| quality for relatively large systems.
To have a simple example of how to use the |CIPSI| program, go to the `users_guide/quickstart`.


The main goal is the development of selected configuration interaction |sCI|
methods and multi-reference perturbation theory |MRPT| in the
determinant-driven paradigm. It also contains the very basics of Kohn-Sham `density functional theory <https://en.wikipedia.org/wiki/Density_functional_theory>`_ |KS-DFT| and `range-separated hybrids <https://aip.scitation.org/doi/10.1063/1.1383587>`_ |RSH|.

The determinant-driven framework allows the programmer to include any arbitrary set of
determinants in the variational space, and thus gives a complete freedom in the methodological
development. The basic ingredients of |RSH| together with those of the |WFT| framework available in the |qp| library allows one to easily develop range-separated DFT (|RSDFT|) approaches (see for instance the plugins at `<https://gitlab.com/eginer/qp_plugins_eginer>`_).

All the programs are developed with the `IRPF90`_ code generator, which considerably simplifies
the collaborative development, and the development of new features.



What it is not
==============

The |qp| is *not* a general purpose quantum chemistry program.
First of all, it is a *library* to develop new theories and algorithms in quantum chemistry.
Therefore, beside the use of the programs of the core modules, the users of the |qp| should develop their own programs.

The |qp| has been designed specifically for |sCI|, so all the
algorithms which are programmed are not adapted to run SCF or DFT calculations
on thousands of atoms. Currently, the systems targeted have less than 600
molecular orbitals. This limit is due to the memory bottleneck induced by the storring of the two-electron integrals (see ``mo_two_e_integrals`` and ``ao_two_e_integrals``).

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

Multiple programs were developed with the |qp|, such as:

- Selected Full-CI + Epstein-Nesbet PT2 (CIPSI) :cite:`Caffarel_2016,Caffarel_2016.2,Loos_2018,Scemama_2018,Dash_2018`
- Hybrid stochastic/deterministic MR-PT2 :cite:`Garniron_2017.2,Loos_2018`
- Orbital optimization for open-shell systems :cite:`Giner2016Mar,Giner_2017.3`
- CIS, CISD, MP2
- Selected CISD
- Jeziorsky-Monkhorst MR-PT2 :cite:`Giner_2017`
- Effective Hamiltonian for variational MR wave functions :cite:`Giner_2017.2`
- Selected CAS+SD
- Selected difference-dedicated CI (DD-CI)
- Multi-Reference Coupled Cluster (MR-CCSD) :cite:`Giner_2016,Garniron_2017`
- Shifted-Bk with CIPSI :cite:`Garniron_2018`
- CIPSI with range-separated DFT (plugins at `<https://gitlab.com/eginer/qp_plugins_eginer>`_)
- DFT for basis set corrections :cite:`Giner_2018`

All these programs can generate ground and excited states, and spin pure wave
functions (eigenstates of |S^2|).



