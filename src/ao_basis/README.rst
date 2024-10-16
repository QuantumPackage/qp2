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



Complex Gaussian-Type Orbitals (cGTOs)
=====================================

Complex Gaussian-Type Orbitals (cGTOs) are also supported:

.. math::

   \chi_i(\mathbf{r}) = x^a y^b z^c \sum_k c_{ki} \left( e^{-\alpha_{ki} \mathbf{r}^2 - \imath \mathbf{k}_{ki} \cdot \mathbf{r} - \imath \phi_{ki}} + \text{C.C.} \right)

where:
   - :math:`\alpha \in \mathbb{C}` and :math:`\Re(\alpha) > 0` (specified by ``ao_expo`` and ``ao_expo_im_cgtos``),
   - :math:`\mathbf{k} = (k_x, k_y, k_z) \in \mathbb{R}^3` (specified by ``ao_expo_pw``),
   - :math:`\phi = \phi_x + \phi_y + \phi_z \in \mathbb{R}` (specified by ``ao_expo_phase``).
