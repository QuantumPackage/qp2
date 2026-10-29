.. _module_becke_numerical_grid: 
 
.. program:: becke_numerical_grid 
 
.. default-role:: option 
 
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


 
 
 
EZFIO parameters 
---------------- 
 
.. option:: grid_type_sgn
 
    Type of grid used for the Becke's numerical grid. Can be, by increasing accuracy: [ 0 | 1 | 2 | 3 ]
 
    Default: 2
 
.. option:: n_points_final_grid
 
    Total number of grid points
 
 
.. option:: thresh_grid
 
    threshold on the weight of a given grid point
 
    Default: 1.e-20
 
.. option:: my_grid_becke
 
    if True, the number of angular and radial grid points are read from EZFIO
 
    Default: False
 
.. option:: my_n_pt_r_grid
 
    Number of radial grid points given from input
 
    Default: 300
 
.. option:: my_n_pt_a_grid
 
    Number of angular grid points given from input. Warning, this number cannot be any integer. See file list_angular_grid
 
    Default: 1202
 
.. option:: n_points_extra_final_grid
 
    Total number of extra_grid points
 
 
.. option:: final_grid_points_extra
 
    Grid points on the extra grid
 
 
.. option:: final_weight_at_r_vector_extra
 
    Weights of the points on the extra grid
 
 
.. option:: extra_grid_type_sgn
 
    Type of extra_grid used for the Becke's numerical extra_grid. Can be, by increasing accuracy: [ 0 | 1 | 2 | 3 ]
 
    Default: 0
 
.. option:: thresh_extra_grid
 
    threshold on the weight of a given extra_grid point
 
    Default: 1.e-20
 
.. option:: my_extra_grid_becke
 
    if True, the number of angular and radial extra_grid points are read from EZFIO
 
    Default: False
 
.. option:: my_n_pt_r_extra_grid
 
    Number of radial extra_grid points given from input
 
    Default: 300
 
.. option:: my_n_pt_a_extra_grid
 
    Number of angular extra_grid points given from input. Warning, this number cannot be any integer. See file list_angular_extra_grid
 
    Default: 1202
 
.. option:: rad_grid_type
 
    method used to sample the radial space. Possible choices are [KNOWLES | GILL]
 
    Default: KNOWLES
 
.. option:: extra_rad_grid_type
 
    method used to sample the radial space. Possible choices are [KNOWLES | GILL]
 
    Default: KNOWLES
 
.. option:: io_extra_grid
 
    Read/Write the extra grid from/to disk [ Write | Read | None ]
 
    Default: None
 
 
Providers 
--------- 
 
.. c:var:: double alpha_knowles


    File : :file:`becke_numerical_grid/integration_radial.irp.f`

    .. code:: fortran

        double precision, allocatable	:: alpha_knowles	(100)


    Recommended values for the alpha parameters according to the paper of Knowles (JCP, 104, 1996)
    as a function of the nuclear charge

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`grid_points_per_atom`

 
.. c:var:: double angular_quadrature_points


    File : :file:`becke_numerical_grid/angular_grid_pts.irp.f`

    .. code:: fortran

        double precision, allocatable	:: angular_quadrature_points	(n_points_integration_angular,3)
        double precision, allocatable	:: weights_angular_points	(n_points_integration_angular)


    weights and grid points for the integration on the angular variables on
    the unit sphere centered on (0,0,0)
    According to the LEBEDEV scheme

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_radial_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`

 
.. c:var:: double angular_quadrature_points_extra


    File : :file:`becke_numerical_grid/angular_extra_grid.irp.f`

    .. code:: fortran

        double precision, allocatable	:: angular_quadrature_points_extra	(n_points_extra_integration_angular,3)
        double precision, allocatable	:: weights_angular_points_extra	(n_points_extra_integration_angular)


    weights and grid points_extra for the integration on the angular variables on
    the unit sphere centered on (0,0,0)
    According to the LEBEDEV scheme

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_extra_radial_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`

 
.. c:var:: double dr_radial_extra_integral


    File : :file:`becke_numerical_grid/extra_grid.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grid_points_extra_radial	(n_points_extra_radial_grid)
        double precision	:: dr_radial_extra_integral	


    points_extra in [0,1] to map the radial integral [0,\infty]

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_extra_radial_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`

 
.. c:var:: double dr_radial_integral


    File : :file:`becke_numerical_grid/grid_becke.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grid_points_radial	(n_points_radial_grid)
        double precision	:: dr_radial_integral	


    points in [0,1] to map the radial integral [0,\infty]

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_radial_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`

 
.. c:var:: double final_grid_points


    File : :file:`becke_numerical_grid/grid_becke_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points	(3,n_points_final_grid)
        double precision, allocatable	:: final_weight_at_r_vector	(n_points_final_grid)
        integer, allocatable	:: index_final_points	(3,n_points_final_grid)
        integer, allocatable	:: index_final_points_reverse	(n_points_integration_angular,n_points_radial_grid,nucl_num)


     final_grid_points(1:3,j) = (/ x, y, z /) of the jth grid point
    
    final_weight_at_r_vector(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
    
    index_final_points(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
    
    index_final_points_reverse(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_abs_int_grid`
       * :c:data:`ao_overlap_abs_grid`
       * :c:data:`ao_prod_abs_r`
       * :c:data:`ao_prod_center`
       * :c:data:`ao_prod_dist_grid`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_lapl_in_r_array`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`final_grid_points_transp`
       * :c:data:`mo_dipole_x_becke`
       * :c:data:`mo_grad_ints`
       * :c:data:`mos_in_r_array`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

 
.. c:var:: double final_grid_points_extra


    File : :file:`becke_numerical_grid/extra_grid_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points_extra	(3,n_points_extra_final_grid)
        double precision, allocatable	:: final_weight_at_r_vector_extra	(n_points_extra_final_grid)
        integer, allocatable	:: index_final_points_extra	(3,n_points_extra_final_grid)
        integer, allocatable	:: index_final_points_extra_reverse	(n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)


     final_grid_points_extra(1:3,j) = (/ x, y, z /) of the jth grid point
    
    final_weight_at_r_vector_extra(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
    
    index_final_points_extra(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
    
    index_final_points_extra_reverse(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_extra_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`read_extra_grid`
       * :c:data:`thresh_extra_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_grad_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`mos_in_r_array_extra_omp`

 
.. c:var:: double final_grid_points_per_atom


    File : :file:`becke_numerical_grid/grid_becke_per_atom.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points_per_atom	(3,n_pts_max_per_atom,nucl_num)
        double precision, allocatable	:: final_weight_at_r_vector_per_atom	(n_pts_max_per_atom,nucl_num)
        integer, allocatable	:: index_final_points_per_atom	(3,n_pts_max_per_atom,nucl_num)
        integer, allocatable	:: index_final_points_per_atom_reverse	(n_points_integration_angular,n_points_radial_grid,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`n_points_radial_grid`
       * :c:data:`n_pts_per_atom`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`


 
.. c:var:: double final_grid_points_transp


    File : :file:`becke_numerical_grid/grid_becke_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points_transp	(n_points_final_grid,3)


    Transposed final_grid_points

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`n_points_final_grid`


 
.. c:var:: double final_weight_at_r


    File : :file:`becke_numerical_grid/grid_becke.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_weight_at_r	(n_points_integration_angular,n_points_radial_grid,nucl_num)


    Total weight on each grid point which takes into account all Lebedev, Voronoi and radial weights.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`alpha_knowles`
       * :c:data:`angular_quadrature_points`
       * :c:data:`grid_atomic_number`
       * :c:data:`grid_points_radial`
       * :c:data:`m_knowles`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`r_gill`
       * :c:data:`rad_grid_type`
       * :c:data:`weight_at_r`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_per_atom`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_pts_per_atom`

 
.. c:var:: double final_weight_at_r_extra


    File : :file:`becke_numerical_grid/extra_grid.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_weight_at_r_extra	(n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)


    Total weight on each grid point which takes into account all Lebedev, Voronoi and radial weights.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`alpha_knowles`
       * :c:data:`angular_quadrature_points_extra`
       * :c:data:`extra_rad_grid_type`
       * :c:data:`grid_atomic_number`
       * :c:data:`grid_points_extra_radial`
       * :c:data:`m_knowles`
       * :c:data:`n_points_extra_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`r_gill`
       * :c:data:`weight_at_r_extra`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points_extra`
       * :c:data:`n_points_extra_final_grid`

 
.. c:var:: double final_weight_at_r_vector


    File : :file:`becke_numerical_grid/grid_becke_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points	(3,n_points_final_grid)
        double precision, allocatable	:: final_weight_at_r_vector	(n_points_final_grid)
        integer, allocatable	:: index_final_points	(3,n_points_final_grid)
        integer, allocatable	:: index_final_points_reverse	(n_points_integration_angular,n_points_radial_grid,nucl_num)


     final_grid_points(1:3,j) = (/ x, y, z /) of the jth grid point
    
    final_weight_at_r_vector(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
    
    index_final_points(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
    
    index_final_points_reverse(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_abs_int_grid`
       * :c:data:`ao_overlap_abs_grid`
       * :c:data:`ao_prod_abs_r`
       * :c:data:`ao_prod_center`
       * :c:data:`ao_prod_dist_grid`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_lapl_in_r_array`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`final_grid_points_transp`
       * :c:data:`mo_dipole_x_becke`
       * :c:data:`mo_grad_ints`
       * :c:data:`mos_in_r_array`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

 
.. c:var:: double final_weight_at_r_vector_extra


    File : :file:`becke_numerical_grid/extra_grid_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points_extra	(3,n_points_extra_final_grid)
        double precision, allocatable	:: final_weight_at_r_vector_extra	(n_points_extra_final_grid)
        integer, allocatable	:: index_final_points_extra	(3,n_points_extra_final_grid)
        integer, allocatable	:: index_final_points_extra_reverse	(n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)


     final_grid_points_extra(1:3,j) = (/ x, y, z /) of the jth grid point
    
    final_weight_at_r_vector_extra(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
    
    index_final_points_extra(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
    
    index_final_points_extra_reverse(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_extra_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`read_extra_grid`
       * :c:data:`thresh_extra_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_grad_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`mos_in_r_array_extra_omp`

 
.. c:var:: double final_weight_at_r_vector_per_atom


    File : :file:`becke_numerical_grid/grid_becke_per_atom.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points_per_atom	(3,n_pts_max_per_atom,nucl_num)
        double precision, allocatable	:: final_weight_at_r_vector_per_atom	(n_pts_max_per_atom,nucl_num)
        integer, allocatable	:: index_final_points_per_atom	(3,n_pts_max_per_atom,nucl_num)
        integer, allocatable	:: index_final_points_per_atom_reverse	(n_points_integration_angular,n_points_radial_grid,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`n_points_radial_grid`
       * :c:data:`n_pts_per_atom`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`


 
.. c:var:: int grid_atomic_number


    File : :file:`becke_numerical_grid/atomic_number.irp.f`

    .. code:: fortran

        integer, allocatable	:: grid_atomic_number	(nucl_num)


    Atomic number used to adjust the grid

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`grid_points_per_atom`

 
.. c:var:: double grid_points_extra_per_atom


    File : :file:`becke_numerical_grid/extra_grid.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grid_points_extra_per_atom	(3,n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)


    x,y,z coordinates of grid points_extra used for integration in 3d space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`alpha_knowles`
       * :c:data:`angular_quadrature_points_extra`
       * :c:data:`extra_rad_grid_type`
       * :c:data:`grid_atomic_number`
       * :c:data:`grid_points_extra_radial`
       * :c:data:`m_knowles`
       * :c:data:`n_points_extra_radial_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`r_gill`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points_extra`
       * :c:data:`weight_at_r_extra`

 
.. c:var:: double grid_points_extra_radial


    File : :file:`becke_numerical_grid/extra_grid.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grid_points_extra_radial	(n_points_extra_radial_grid)
        double precision	:: dr_radial_extra_integral	


    points_extra in [0,1] to map the radial integral [0,\infty]

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_extra_radial_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`

 
.. c:var:: double grid_points_per_atom


    File : :file:`becke_numerical_grid/grid_becke.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grid_points_per_atom	(3,n_points_integration_angular,n_points_radial_grid,nucl_num)


    x,y,z coordinates of grid points used for integration in 3d space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`alpha_knowles`
       * :c:data:`angular_quadrature_points`
       * :c:data:`grid_atomic_number`
       * :c:data:`grid_points_radial`
       * :c:data:`m_knowles`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`r_gill`
       * :c:data:`rad_grid_type`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_per_atom`
       * :c:data:`weight_at_r`

 
.. c:var:: double grid_points_radial


    File : :file:`becke_numerical_grid/grid_becke.irp.f`

    .. code:: fortran

        double precision, allocatable	:: grid_points_radial	(n_points_radial_grid)
        double precision	:: dr_radial_integral	


    points in [0,1] to map the radial integral [0,\infty]

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_radial_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`

 
.. c:var:: int index_final_points


    File : :file:`becke_numerical_grid/grid_becke_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points	(3,n_points_final_grid)
        double precision, allocatable	:: final_weight_at_r_vector	(n_points_final_grid)
        integer, allocatable	:: index_final_points	(3,n_points_final_grid)
        integer, allocatable	:: index_final_points_reverse	(n_points_integration_angular,n_points_radial_grid,nucl_num)


     final_grid_points(1:3,j) = (/ x, y, z /) of the jth grid point
    
    final_weight_at_r_vector(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
    
    index_final_points(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
    
    index_final_points_reverse(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_abs_int_grid`
       * :c:data:`ao_overlap_abs_grid`
       * :c:data:`ao_prod_abs_r`
       * :c:data:`ao_prod_center`
       * :c:data:`ao_prod_dist_grid`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_lapl_in_r_array`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`final_grid_points_transp`
       * :c:data:`mo_dipole_x_becke`
       * :c:data:`mo_grad_ints`
       * :c:data:`mos_in_r_array`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

 
.. c:var:: int index_final_points_extra


    File : :file:`becke_numerical_grid/extra_grid_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points_extra	(3,n_points_extra_final_grid)
        double precision, allocatable	:: final_weight_at_r_vector_extra	(n_points_extra_final_grid)
        integer, allocatable	:: index_final_points_extra	(3,n_points_extra_final_grid)
        integer, allocatable	:: index_final_points_extra_reverse	(n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)


     final_grid_points_extra(1:3,j) = (/ x, y, z /) of the jth grid point
    
    final_weight_at_r_vector_extra(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
    
    index_final_points_extra(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
    
    index_final_points_extra_reverse(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_extra_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`read_extra_grid`
       * :c:data:`thresh_extra_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_grad_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`mos_in_r_array_extra_omp`

 
.. c:var:: int index_final_points_extra_reverse


    File : :file:`becke_numerical_grid/extra_grid_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points_extra	(3,n_points_extra_final_grid)
        double precision, allocatable	:: final_weight_at_r_vector_extra	(n_points_extra_final_grid)
        integer, allocatable	:: index_final_points_extra	(3,n_points_extra_final_grid)
        integer, allocatable	:: index_final_points_extra_reverse	(n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)


     final_grid_points_extra(1:3,j) = (/ x, y, z /) of the jth grid point
    
    final_weight_at_r_vector_extra(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
    
    index_final_points_extra(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
    
    index_final_points_extra_reverse(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_extra_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`read_extra_grid`
       * :c:data:`thresh_extra_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_grad_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`mos_in_r_array_extra_omp`

 
.. c:var:: int index_final_points_per_atom


    File : :file:`becke_numerical_grid/grid_becke_per_atom.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points_per_atom	(3,n_pts_max_per_atom,nucl_num)
        double precision, allocatable	:: final_weight_at_r_vector_per_atom	(n_pts_max_per_atom,nucl_num)
        integer, allocatable	:: index_final_points_per_atom	(3,n_pts_max_per_atom,nucl_num)
        integer, allocatable	:: index_final_points_per_atom_reverse	(n_points_integration_angular,n_points_radial_grid,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`n_points_radial_grid`
       * :c:data:`n_pts_per_atom`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`


 
.. c:var:: int index_final_points_per_atom_reverse


    File : :file:`becke_numerical_grid/grid_becke_per_atom.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points_per_atom	(3,n_pts_max_per_atom,nucl_num)
        double precision, allocatable	:: final_weight_at_r_vector_per_atom	(n_pts_max_per_atom,nucl_num)
        integer, allocatable	:: index_final_points_per_atom	(3,n_pts_max_per_atom,nucl_num)
        integer, allocatable	:: index_final_points_per_atom_reverse	(n_points_integration_angular,n_points_radial_grid,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`n_points_radial_grid`
       * :c:data:`n_pts_per_atom`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`


 
.. c:var:: int index_final_points_reverse


    File : :file:`becke_numerical_grid/grid_becke_vector.irp.f`

    .. code:: fortran

        double precision, allocatable	:: final_grid_points	(3,n_points_final_grid)
        double precision, allocatable	:: final_weight_at_r_vector	(n_points_final_grid)
        integer, allocatable	:: index_final_points	(3,n_points_final_grid)
        integer, allocatable	:: index_final_points_reverse	(n_points_integration_angular,n_points_radial_grid,nucl_num)


     final_grid_points(1:3,j) = (/ x, y, z /) of the jth grid point
    
    final_weight_at_r_vector(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
    
    index_final_points(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
    
    index_final_points_reverse(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_abs_int_grid`
       * :c:data:`ao_overlap_abs_grid`
       * :c:data:`ao_prod_abs_r`
       * :c:data:`ao_prod_center`
       * :c:data:`ao_prod_dist_grid`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_lapl_in_r_array`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`final_grid_points_transp`
       * :c:data:`mo_dipole_x_becke`
       * :c:data:`mo_grad_ints`
       * :c:data:`mos_in_r_array`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

 
.. c:var:: int m_knowles


    File : :file:`becke_numerical_grid/grid_becke.irp.f`

    .. code:: fortran

        integer	:: m_knowles	


    value of the "m" parameter in the equation (7) of the paper of Knowles (JCP, 104, 1996)

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`grid_points_per_atom`

 
.. c:var:: int n_points_extra_final_grid


    File : :file:`becke_numerical_grid/extra_grid_vector.irp.f`

    .. code:: fortran

        integer	:: n_points_extra_final_grid	


    Number of points_extra which are non zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r_extra`
       * :c:data:`n_points_extra_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`read_extra_grid`
       * :c:data:`thresh_extra_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`aos_grad_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`aos_in_r_array_extra_transp`
       * :c:data:`final_grid_points_extra`
       * :c:data:`mos_in_r_array_extra_omp`
       * :c:data:`mos_in_r_array_extra_transp`

 
.. c:var:: int n_points_extra_grid_per_atom


    File : :file:`becke_numerical_grid/extra_grid.irp.f`

    .. code:: fortran

        integer	:: n_points_extra_grid_per_atom	


    Number of grid points_extra per atom

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_extra_radial_grid`


 
.. c:var:: int n_points_extra_integration_angular


    File : :file:`becke_numerical_grid/extra_grid.irp.f`

    .. code:: fortran

        integer	:: n_points_extra_radial_grid	
        integer	:: n_points_extra_integration_angular	


    n_points_extra_radial_grid = number of radial grid points_extra per atom
    
    n_points_extra_integration_angular = number of angular grid points_extra per atom
    
    These numbers are automatically set by setting the grid_type_sgn parameter

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`extra_grid_type_sgn`
       * :c:data:`my_extra_grid_becke`
       * :c:data:`my_n_pt_a_extra_grid`
       * :c:data:`my_n_pt_r_extra_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`angular_quadrature_points_extra`
       * :c:data:`final_grid_points_extra`
       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`grid_points_extra_radial`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_extra_grid_per_atom`
       * :c:data:`weight_at_r_extra`

 
.. c:var:: int n_points_extra_radial_grid


    File : :file:`becke_numerical_grid/extra_grid.irp.f`

    .. code:: fortran

        integer	:: n_points_extra_radial_grid	
        integer	:: n_points_extra_integration_angular	


    n_points_extra_radial_grid = number of radial grid points_extra per atom
    
    n_points_extra_integration_angular = number of angular grid points_extra per atom
    
    These numbers are automatically set by setting the grid_type_sgn parameter

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`extra_grid_type_sgn`
       * :c:data:`my_extra_grid_becke`
       * :c:data:`my_n_pt_a_extra_grid`
       * :c:data:`my_n_pt_r_extra_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`angular_quadrature_points_extra`
       * :c:data:`final_grid_points_extra`
       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`grid_points_extra_radial`
       * :c:data:`n_points_extra_final_grid`
       * :c:data:`n_points_extra_grid_per_atom`
       * :c:data:`weight_at_r_extra`

 
.. c:var:: int n_points_final_grid


    File : :file:`becke_numerical_grid/grid_becke_vector.irp.f`

    .. code:: fortran

        integer	:: n_points_final_grid	


    Number of points which are non zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`alpha_dens_kin_in_r`
       * :c:data:`ao_abs_int_grid`
       * :c:data:`ao_overlap_abs_grid`
       * :c:data:`ao_prod_abs_r`
       * :c:data:`ao_prod_center`
       * :c:data:`ao_prod_dist_grid`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_grad_in_r_array_transp`
       * :c:data:`aos_grad_in_r_array_transp_3`
       * :c:data:`aos_grad_in_r_array_transp_bis`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_in_r_array_transp`
       * :c:data:`aos_lapl_in_r_array`
       * :c:data:`aos_lapl_in_r_array_transp`
       * :c:data:`elec_beta_num_grid_becke`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_transp`
       * :c:data:`kinetic_density_generalized`
       * :c:data:`mo_dipole_x_becke`
       * :c:data:`mo_grad_ints`
       * :c:data:`mos_grad_in_r_array`
       * :c:data:`mos_grad_in_r_array_tranp`
       * :c:data:`mos_grad_in_r_array_transp_3`
       * :c:data:`mos_grad_in_r_array_transp_bis`
       * :c:data:`mos_in_r_array`
       * :c:data:`mos_in_r_array_omp`
       * :c:data:`mos_in_r_array_transp`
       * :c:data:`mos_lapl_in_r_array`
       * :c:data:`mos_lapl_in_r_array_tranp`
       * :c:data:`one_e_dm_and_grad_alpha_in_r`

 
.. c:var:: int n_points_grid_per_atom


    File : :file:`becke_numerical_grid/grid_becke.irp.f`

    .. code:: fortran

        integer	:: n_points_grid_per_atom	


    Number of grid points per atom

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_radial_grid`


 
.. c:var:: int n_points_integration_angular


    File : :file:`becke_numerical_grid/grid_becke.irp.f`

    .. code:: fortran

        integer	:: n_points_radial_grid	
        integer	:: n_points_integration_angular	


    n_points_radial_grid = number of radial grid points per atom
    
    n_points_integration_angular = number of angular grid points per atom
    
    These numbers are automatically set by setting the grid_type_sgn parameter

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`grid_type_sgn`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`angular_quadrature_points`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_per_atom`
       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`grid_points_radial`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_points_grid_per_atom`
       * :c:data:`n_pts_per_atom`
       * :c:data:`weight_at_r`

 
.. c:var:: int n_points_radial_grid


    File : :file:`becke_numerical_grid/grid_becke.irp.f`

    .. code:: fortran

        integer	:: n_points_radial_grid	
        integer	:: n_points_integration_angular	


    n_points_radial_grid = number of radial grid points per atom
    
    n_points_integration_angular = number of angular grid points per atom
    
    These numbers are automatically set by setting the grid_type_sgn parameter

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`grid_type_sgn`
       * :c:data:`my_grid_becke`
       * :c:data:`my_n_pt_a_grid`
       * :c:data:`my_n_pt_r_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`angular_quadrature_points`
       * :c:data:`final_grid_points`
       * :c:data:`final_grid_points_per_atom`
       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`grid_points_radial`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_points_grid_per_atom`
       * :c:data:`n_pts_per_atom`
       * :c:data:`weight_at_r`

 
.. c:var:: int n_pts_max_per_atom


    File : :file:`becke_numerical_grid/grid_becke_per_atom.irp.f`

    .. code:: fortran

        integer, allocatable	:: n_pts_per_atom	(nucl_num)
        integer	:: n_pts_max_per_atom	


    Number of points which are non zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points_per_atom`

 
.. c:var:: int n_pts_per_atom


    File : :file:`becke_numerical_grid/grid_becke_per_atom.irp.f`

    .. code:: fortran

        integer, allocatable	:: n_pts_per_atom	(nucl_num)
        integer	:: n_pts_max_per_atom	


    Number of points which are non zero

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_num`
       * :c:data:`thresh_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points_per_atom`

 
.. c:var:: double r_gill


    File : :file:`becke_numerical_grid/grid_becke.irp.f`

    .. code:: fortran

        double precision	:: r_gill	



    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`grid_points_per_atom`

 
.. c:var:: double weight_at_r


    File : :file:`becke_numerical_grid/grid_becke.irp.f`

    .. code:: fortran

        double precision, allocatable	:: weight_at_r	(n_points_integration_angular,n_points_radial_grid,nucl_num)


    Weight function at grid points : w_n(r) according to the equation (22)
    of Becke original paper (JCP, 88, 1988)
    
    The "n" discrete variable represents the nucleis which in this array is
    represented by the last dimension and the points are labelled by the
    other dimensions.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`grid_points_per_atom`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_coord_transp`
       * :c:data:`nucl_dist_inv`
       * :c:data:`nucl_num`
       * :c:data:`slater_bragg_type_inter_distance_ua`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`

 
.. c:var:: double weight_at_r_extra


    File : :file:`becke_numerical_grid/extra_grid.irp.f`

    .. code:: fortran

        double precision, allocatable	:: weight_at_r_extra	(n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)


    Weight function at grid points_extra : w_n(r) according to the equation (22)
    of Becke original paper (JCP, 88, 1988)
    
    The "n" discrete variable represents the nucleis which in this array is
    represented by the last dimension and the points_extra are labelled by the
    other dimensions.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`n_points_extra_radial_grid`
       * :c:data:`nucl_coord_transp`
       * :c:data:`nucl_dist_inv`
       * :c:data:`nucl_num`
       * :c:data:`slater_bragg_type_inter_distance_ua`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r_extra`

 
.. c:var:: double weights_angular_points


    File : :file:`becke_numerical_grid/angular_grid_pts.irp.f`

    .. code:: fortran

        double precision, allocatable	:: angular_quadrature_points	(n_points_integration_angular,3)
        double precision, allocatable	:: weights_angular_points	(n_points_integration_angular)


    weights and grid points for the integration on the angular variables on
    the unit sphere centered on (0,0,0)
    According to the LEBEDEV scheme

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_radial_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`

 
.. c:var:: double weights_angular_points_extra


    File : :file:`becke_numerical_grid/angular_extra_grid.irp.f`

    .. code:: fortran

        double precision, allocatable	:: angular_quadrature_points_extra	(n_points_extra_integration_angular,3)
        double precision, allocatable	:: weights_angular_points_extra	(n_points_extra_integration_angular)


    weights and grid points_extra for the integration on the angular variables on
    the unit sphere centered on (0,0,0)
    According to the LEBEDEV scheme

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_points_extra_radial_grid`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`final_weight_at_r_extra`
       * :c:data:`grid_points_extra_per_atom`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void cell_function_becke(void)


    File : :file:`becke_numerical_grid/step_function_becke.irp.f`

    .. code:: fortran

        double precision function cell_function_becke(r, atom_number)


    atom_number :: atom on which the cell function of Becke (1988, JCP,88(4))
    r(1:3)                       :: x,y,z coordinantes of the current point

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord_transp`
       * :c:data:`nucl_dist_inv`
       * :c:data:`nucl_num`
       * :c:data:`slater_bragg_type_inter_distance_ua`

 
.. c:function:: void derivative_knowles_function(void)


    File : :file:`becke_numerical_grid/integration_radial.irp.f`

    .. code:: fortran

        double precision function derivative_knowles_function(alpha, m, x)


    Derivative of the function proposed by Knowles (JCP, 104, 1996) for distributing the radial points

 
.. c:function:: void example_becke_numerical_grid(void)


    File : :file:`becke_numerical_grid/example.irp.f`

    .. code:: fortran

        subroutine example_becke_numerical_grid


    subroutine that illustrates the main features available in becke_numerical_grid

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`final_grid_points`
       * :c:data:`final_weight_at_r`
       * :c:data:`grid_points_per_atom`
       * :c:data:`n_points_final_grid`
       * :c:data:`n_points_radial_grid`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

 
.. c:function:: void f_function_becke(void)


    File : :file:`becke_numerical_grid/step_function_becke.irp.f`

    .. code:: fortran

        double precision function f_function_becke(x)



 
.. c:function:: void knowles_function(void)


    File : :file:`becke_numerical_grid/integration_radial.irp.f`

    .. code:: fortran

        double precision function knowles_function(alpha, m, x)


    Function proposed by Knowles (JCP, 104, 1996) for distributing the radial points :
    the Log "m" function ( equation (7) in the paper )

 
.. c:function:: void step_function_becke(void)


    File : :file:`becke_numerical_grid/step_function_becke.irp.f`

    .. code:: fortran

        double precision function step_function_becke(x)


    Step function of the Becke paper (1988, JCP,88(4))

