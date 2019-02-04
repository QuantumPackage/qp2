========
ao_basis
========

This module describes the atomic orbitals basis set.

An |AO| :math:`\chi` centered on nucleus A is represented as:

.. math::

   \chi_i({\bf r}) = (x-X_A)^a (y-Y_A)^b (z-Z_A)^c \sum_k c_{ki} e^{-\gamma_{ki} |{\bf r} - {\bf R}_A|^2}


The |AO| coefficients are normalized as:

.. math::

  {\tilde c}_{ki} = \frac{c_{ki}}{ \int \left( (x-X_A)^a (y-Y_A)^b (z-Z_A)^c  e^{-\gamma_{ki} |{\bf r} - {\bf R}_A|^2} \right)^2 dr}


.. warning::

  `ao_coef` contains the |AO| coefficients given in input. These do not
  include the normalization constant of the |AO|. The `ao_coef_normalized`
  provider includes this normalization factor.


The |AOs| are also sorted by increasing exponent to accelerate the calculation of
the two electron integrals.



