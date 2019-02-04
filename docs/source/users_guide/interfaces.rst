Interfaces
==========

.. TODO

A few interfaces to external codes are available.

\* -> |qp|
----------

`GAMESS`_ / Gaussian 
  Using the |resultsFile| Python library, the geometry and |MOs| can be read.
  This is useful to make calculations with |CASSCF| orbitals

|qp| -> \* 
----------

`Molden <http://cheminf.cmbi.ru.nl/molden>`_
  3D plots of Molecular Orbitals

FCIDUMP 
  Interface with the |FCI| - |QMC| program `NECI`_, or the semi-stochastic
  Heat-Bath |CI| program `Dice`_.

`QMCPack`_ / `CHAMP <https://www.utwente.nl/en/tnw/ccp/research/CHAMP.html>`_ /
`QMC=Chem`_
Trial wave functions can be used for |QMC|, with or without pseudo-potentials.
These interfaces are provided as `external plugins`_.


