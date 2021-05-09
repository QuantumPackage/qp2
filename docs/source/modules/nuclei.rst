.. _module_nuclei: 
 
.. program:: nuclei 
 
.. default-role:: option 
 
======
nuclei
======

This module contains data relative to the nuclei (coordinates, charge,
nuclear repulsion energy, etc).
The coordinates are expressed in atomic units.

 
 
 
EZFIO parameters 
---------------- 
 
.. option:: nucl_num
 
    Number of nuclei
 
 
.. option:: nucl_label
 
    Nuclear labels
 
 
.. option:: nucl_charge
 
    Nuclear charges
 
 
.. option:: nucl_coord
 
    Nuclear coordinates in the format (:, {x,y,z})
 
 
.. option:: io_nuclear_repulsion
 
    Read/Write Nuclear Repulsion from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: nuclear_repulsion
 
    Nuclear repulsion (Computed automaticaly or Read in the |EZFIO|)
 
 
 
Providers 
--------- 
 
.. c:var:: center_of_mass


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: center_of_mass	(3)


    Center of mass of the molecule

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`element_name`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`inertia_tensor`

 
.. c:var:: element_mass


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        character*(4), allocatable	:: element_name	(0:127)
        double precision, allocatable	:: element_mass	(0:127)


    Array of the name of element, sorted by nuclear charge (integer)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`center_of_mass`
       * :c:data:`inertia_tensor`

 
.. c:var:: element_name


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        character*(4), allocatable	:: element_name	(0:127)
        double precision, allocatable	:: element_mass	(0:127)


    Array of the name of element, sorted by nuclear charge (integer)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`center_of_mass`
       * :c:data:`inertia_tensor`

 
.. c:var:: inertia_tensor


    File : :file:`nuclei/inertia.irp.f`

    .. code:: fortran

        double precision, allocatable	:: inertia_tensor	(3,3)


    Inertia tensor

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`center_of_mass`
       * :c:data:`element_name`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`inertia_tensor_eigenvectors`

 
.. c:var:: inertia_tensor_eigenvalues


    File : :file:`nuclei/inertia.irp.f`

    .. code:: fortran

        double precision, allocatable	:: inertia_tensor_eigenvectors	(3,3)
        double precision, allocatable	:: inertia_tensor_eigenvalues	(3)


    Eigenvectors/eigenvalues of the inertia_tensor. Used to find normal orientation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`inertia_tensor`


 
.. c:var:: inertia_tensor_eigenvectors


    File : :file:`nuclei/inertia.irp.f`

    .. code:: fortran

        double precision, allocatable	:: inertia_tensor_eigenvectors	(3,3)
        double precision, allocatable	:: inertia_tensor_eigenvalues	(3)


    Eigenvectors/eigenvalues of the inertia_tensor. Used to find normal orientation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`inertia_tensor`


 
.. c:var:: nucl_coord


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nucl_coord	(nucl_num,3)


    Nuclear coordinates in the format (:, {x,y,z})

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_label`
       * :c:data:`nucl_num`
       * :c:data:`output_wall_time_0`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_deriv_1_x`
       * :c:data:`ao_dipole_x`
       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_integrals_n_e_per_atom`
       * :c:data:`ao_overlap`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_pseudo_integrals_local`
       * :c:data:`ao_pseudo_integrals_non_local`
       * :c:data:`ao_spread_x`
       * :c:data:`ao_two_e_integral_schwartz`
       * :c:data:`ao_two_e_integrals_in_map`
       * :c:data:`center_of_mass`
       * :c:data:`inertia_tensor`
       * :c:data:`nucl_coord_transp`
       * :c:data:`nucl_dist_2`
       * :c:data:`nuclear_repulsion`

 
.. c:var:: nucl_coord_transp


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nucl_coord_transp	(3,nucl_num)


    Transposed array of nucl_coord

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`


 
.. c:var:: nucl_dist


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nucl_dist_2	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_x	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_y	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_z	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist	(nucl_num,nucl_num)


    nucl_dist     : Nucleus-nucleus distances
    nucl_dist_2   : Nucleus-nucleus distances squared
    nucl_dist_vec : Nucleus-nucleus distances vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_dist_inv`

 
.. c:var:: nucl_dist_2


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nucl_dist_2	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_x	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_y	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_z	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist	(nucl_num,nucl_num)


    nucl_dist     : Nucleus-nucleus distances
    nucl_dist_2   : Nucleus-nucleus distances squared
    nucl_dist_vec : Nucleus-nucleus distances vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_dist_inv`

 
.. c:var:: nucl_dist_inv


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nucl_dist_inv	(nucl_num,nucl_num)


    Inverse of the distance between nucleus I and nucleus J

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_dist_2`
       * :c:data:`nucl_num`


 
.. c:var:: nucl_dist_vec_x


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nucl_dist_2	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_x	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_y	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_z	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist	(nucl_num,nucl_num)


    nucl_dist     : Nucleus-nucleus distances
    nucl_dist_2   : Nucleus-nucleus distances squared
    nucl_dist_vec : Nucleus-nucleus distances vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_dist_inv`

 
.. c:var:: nucl_dist_vec_y


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nucl_dist_2	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_x	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_y	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_z	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist	(nucl_num,nucl_num)


    nucl_dist     : Nucleus-nucleus distances
    nucl_dist_2   : Nucleus-nucleus distances squared
    nucl_dist_vec : Nucleus-nucleus distances vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_dist_inv`

 
.. c:var:: nucl_dist_vec_z


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nucl_dist_2	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_x	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_y	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist_vec_z	(nucl_num,nucl_num)
        double precision, allocatable	:: nucl_dist	(nucl_num,nucl_num)


    nucl_dist     : Nucleus-nucleus distances
    nucl_dist_2   : Nucleus-nucleus distances squared
    nucl_dist_vec : Nucleus-nucleus distances vectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_dist_inv`

 
.. c:var:: nuclear_repulsion


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision	:: nuclear_repulsion	


    Nuclear repulsion energy

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`io_nuclear_repulsion`
       * :c:data:`mpi_master`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`output_wall_time_0`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy`
       * :c:data:`core_energy`
       * :c:data:`psi_energy_with_nucl_rep`

 
.. c:var:: slater_bragg_radii


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_radii	(0:100)


    atomic radii in Angstrom defined in table I of JCP 41, 3199 (1964) Slater
    execpt for the Hydrogen atom where we took the value of Becke (1988, JCP)

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`slater_bragg_radii_per_atom`
       * :c:data:`slater_bragg_radii_ua`

 
.. c:var:: slater_bragg_radii_per_atom


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_radii_per_atom	(nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_num`
       * :c:data:`slater_bragg_radii`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`slater_bragg_type_inter_distance`

 
.. c:var:: slater_bragg_radii_per_atom_ua


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_radii_per_atom_ua	(nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_charge`
       * :c:data:`nucl_num`
       * :c:data:`slater_bragg_radii_ua`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`slater_bragg_type_inter_distance_ua`

 
.. c:var:: slater_bragg_radii_ua


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_radii_ua	(0:100)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`slater_bragg_radii`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`slater_bragg_radii_per_atom_ua`

 
.. c:var:: slater_bragg_type_inter_distance


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_type_inter_distance	(nucl_num,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`slater_bragg_radii_per_atom`


 
.. c:var:: slater_bragg_type_inter_distance_ua


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_type_inter_distance_ua	(nucl_num,nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`slater_bragg_radii_per_atom_ua`


