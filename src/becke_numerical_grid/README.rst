====================
becke_numerical_grid
====================

This module contains all quantities needed to build Becke's grid used in general for DFT integration. Note that it can be used for whatever integration in R^3 as long as the functions to be integrated are mostly concentrated near the atomic regions.

This grid is built as the reunion of a spherical grid around each atom. Each spherical grid contains
a certain number of radial and angular points. No pruning is done on the angular part of the grid.

The main keyword for that module is:

* :option:`becke_numerical_grid grid_type_sgn` which controls the precision of the grid according the standard **SG-n** grids. This keyword controls the two providers `n_points_integration_angular` `n_points_radial_grid`.

The main providers of that module are:

* `n_points_integration_angular` which is the number of angular integration points. WARNING: it obeys to specific rules so it cannot be any integer number. Some of the possible values are [ 50 | 74 | 170 | 194 | 266 | 302 | 590 | 1202 | 2030 | 5810 ] for instance. See :file:`angular.f` for more details.
* `n_points_radial_grid` which is the number of radial angular points. This can be any strictly positive integer. Nevertheless, a minimum of 50 is in general necessary.
* `final_grid_points` which are the (x,y,z) coordinates of the grid points.
* `final_weight_at_r_vector` which are the weights at each grid point


For a simple example of how to use the grid, see :file:`example.irp.f`.

The spherical integration uses Lebedev-Laikov grids, which was used from the code distributed through CCL (http://www.ccl.net/).
See next section for explanations and citation policies.

.. code-block:: text

       This subroutine is part of a set of subroutines that generate
       Lebedev grids [1-6] for integration on a sphere. The original
       C-code [1] was kindly provided by Dr. Dmitri N. Laikov and
       translated into fortran by Dr. Christoph van Wuellen.
       This subroutine was translated using a C to fortran77 conversion
       tool written by Dr. Christoph van Wuellen.

       Users of this code are asked to include reference [1] in their
       publications, and in the user- and programmers-manuals
       describing their codes.

       This code was distributed through CCL (http://www.ccl.net/).

       [1] V.I. Lebedev, and D.N. Laikov
           "A quadrature formula for the sphere of the 131st
            algebraic order of accuracy"
           Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.

       [2] V.I. Lebedev
           "A quadrature formula for the sphere of 59th algebraic
            order of accuracy"
           Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286.

       [3] V.I. Lebedev, and A.L. Skorokhodov
           "Quadrature formulas of orders 41, 47, and 53 for the sphere"
           Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592.

       [4] V.I. Lebedev
           "Spherical quadrature formulas exact to orders 25-29"
           Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.

       [5] V.I. Lebedev
           "Quadratures on a sphere"
           Computational Mathematics and Mathematical Physics, Vol. 16,
           1976, pp. 10-24.

       [6] V.I. Lebedev
           "Values of the nodes and weights of ninth to seventeenth
            order Gauss-Markov quadrature formulae invariant under the
            octahedron group with inversion"
           Computational Mathematics and Mathematical Physics, Vol. 15,
           1975, pp. 44-51.


