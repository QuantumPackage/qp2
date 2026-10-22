==============
generators_cas
==============

Generators for the perturbation energy, as the determinants of a
|CAS|.

The generators are the determinants of the wave function that belong
to the |CAS|: the |MOs| of the |CAS| are the ones set as active with
the :ref:`qp_set_mo_class` command. The subroutine
:c:func:`extract_cas` replaces the total wave function by its
normalized projection on the |CAS|.

This module is intended to be included in the :file:`NEED` file to
define the generators as the |CAS| determinants, which is useful to
define post-CAS approaches (see the cassd module for instance).
