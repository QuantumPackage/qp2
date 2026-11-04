.. _module_nuclei: 
 
.. program:: nuclei 
 
.. default-role:: option 
 
======
nuclei
======

Data relative to the nuclei: coordinates, charges, and the quantities
derived from them.

The nuclear coordinates `nucl_coord`, in atomic units, and the nuclear
charges `nucl_charge` are read from the |EZFIO| database, together
with the number of nuclei :option:`nuclei nucl_num` and the element
labels `nucl_label`. From them this module derives:

* the pair distances :c:data:`nucl_dist` and :c:data:`nucl_dist_2`,
  their inverse :c:data:`nucl_dist_inv`, and the vector components
  :c:data:`nucl_dist_vec_x`, :c:data:`nucl_dist_vec_y`,
  :c:data:`nucl_dist_vec_z`
* the nuclear repulsion energy :c:data:`nuclear_repulsion`, either
  read from the |EZFIO| database or computed from the charges and the
  distances, depending on `io_nuclear_repulsion`
* the element names :c:data:`element_name` and masses
  :c:data:`element_mass`
* the center of mass :c:data:`center_of_mass` and the inertia tensor
  :c:data:`inertia_tensor`, whose diagonalization gives the rotational
  constants
* the Slater-Bragg atomic radii :c:data:`slater_bragg_radii`, and the
  per-atom and pair-distance variants
* the point charges added to the potential:
  :c:data:`pts_charge_z`, :c:data:`pts_charge_coord`, and their
  interaction energy with the electrons and with the nuclei
 
 
 
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
 
 
.. option:: is_periodic
 
    If true, the calculation uses periodic boundary conditions
 
    Default: false
 
.. option:: n_pts_charge
 
    Number of point charges to be added to the potential
 
    Default: 0
 
.. option:: pts_charge_z
 
    Charge associated to each point charge
 
 
.. option:: pts_charge_coord
 
    Coordinate of each point charge.
 
 
.. option:: point_charges
 
    If |true|, point charges (see nuclei/write_pt_charges.py) are added to the one-electron potential
 
    Default: False
 
 
Providers 
--------- 
 
.. c:var:: double center_of_mass


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
       * :c:data:`pg_center`
       * :c:data:`pg_inertia_eval`

 
.. c:var:: double element_mass


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
       * :c:data:`pg_inertia_eval`

 
.. c:var:: char element_name


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
       * :c:data:`pg_inertia_eval`

 
.. c:var:: double inertia_tensor


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

 
.. c:var:: double inertia_tensor_eigenvalues


    File : :file:`nuclei/inertia.irp.f`

    .. code:: fortran

        double precision, allocatable	:: inertia_tensor_eigenvectors	(3,3)
        double precision, allocatable	:: inertia_tensor_eigenvalues	(3)


    Eigenvectors/eigenvalues of the inertia_tensor. Used to find normal orientation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`inertia_tensor`


 
.. c:var:: double inertia_tensor_eigenvectors


    File : :file:`nuclei/inertia.irp.f`

    .. code:: fortran

        double precision, allocatable	:: inertia_tensor_eigenvectors	(3,3)
        double precision, allocatable	:: inertia_tensor_eigenvalues	(3)


    Eigenvectors/eigenvalues of the inertia_tensor. Used to find normal orientation.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`inertia_tensor`


 
.. c:var:: int n_pts_charge


    File : :file:`nuclei/point_charges.irp.f`

    .. code:: fortran

        integer	:: n_pts_charge	


    Number of point charges to be added to the potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`output_wall_time_0`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_pt_chrg`
       * :c:data:`pt_chrg_interaction`
       * :c:data:`pt_chrg_nuclei_interaction`
       * :c:data:`pts_charge_coord`
       * :c:data:`pts_charge_z`

 
.. c:var:: double nucl_coord


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

       * :c:data:`ao_2e_cgtos_schwartz`
       * :c:data:`ao_coef_norm_cgtos`
       * :c:data:`ao_deriv2_cgtos_x`
       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_deriv_1_x`
       * :c:data:`ao_dipole_x`
       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_integrals_n_e_cgtos`
       * :c:data:`ao_integrals_n_e_per_atom`
       * :c:data:`ao_integrals_pt_chrg`
       * :c:data:`ao_overlap`
       * :c:data:`ao_overlap_abs`
       * :c:data:`ao_overlap_cgtos`
       * :c:data:`ao_pseudo_integrals_local`
       * :c:data:`ao_pseudo_integrals_non_local`
       * :c:data:`ao_spread_x`
       * :c:data:`ao_two_e_integral_erf_schwartz`
       * :c:data:`ao_two_e_integrals_erf_in_map`
       * :c:data:`aos_grad_in_r_array`
       * :c:data:`aos_grad_in_r_array_extra`
       * :c:data:`aos_in_r_array`
       * :c:data:`aos_in_r_array_extra`
       * :c:data:`aos_lapl_in_r_array`
       * :c:data:`center_of_mass`
       * :c:data:`cholesky_ao_num`
       * :c:data:`grid_points_extra_per_atom`
       * :c:data:`grid_points_per_atom`
       * :c:data:`inertia_tensor`
       * :c:data:`multi_s_deriv_1`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`nucl_coord_transp`
       * :c:data:`nucl_dist_2`
       * :c:data:`nuclear_repulsion`
       * :c:data:`pg_geom_residual`
       * :c:data:`pg_inertia_eval`
       * :c:data:`pg_is_linear`
       * :c:data:`pg_is_planar`
       * :c:data:`pg_n_elem`
       * :c:data:`pg_n_max`
       * :c:data:`pg_rot`
       * :c:data:`pg_rot_mat`
       * :c:data:`pt_chrg_nuclei_interaction`
       * :c:data:`z_dipole_moment`

 
.. c:var:: double nucl_coord_transp


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nucl_coord_transp	(3,nucl_num)


    Transposed array of nucl_coord

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`weight_at_r`
       * :c:data:`weight_at_r_extra`

 
.. c:var:: double nucl_dist


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

 
.. c:var:: double nucl_dist_2


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

 
.. c:var:: double nucl_dist_inv


    File : :file:`nuclei/nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: nucl_dist_inv	(nucl_num,nucl_num)


    Inverse of the distance between nucleus I and nucleus J

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_dist_2`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`weight_at_r`
       * :c:data:`weight_at_r_extra`

 
.. c:var:: double nucl_dist_vec_x


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

 
.. c:var:: double nucl_dist_vec_y


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

 
.. c:var:: double nucl_dist_vec_z


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

 
.. c:var:: double nuclear_repulsion


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
       * :c:data:`point_charges`
       * :c:data:`pt_chrg_nuclei_interaction`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy`
       * :c:data:`ci_energy_no_diag`
       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`psi_energy_with_nucl_rep`

 
.. c:var:: double pt_chrg_interaction


    File : :file:`nuclei/point_charges.irp.f`

    .. code:: fortran

        double precision	:: pt_chrg_interaction	


    Interaction between the point charges

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_pts_charge`
       * :c:data:`pts_charge_coord`
       * :c:data:`pts_charge_z`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`pt_chrg_nuclei_interaction`

 
.. c:var:: double pt_chrg_nuclei_interaction


    File : :file:`nuclei/point_charges.irp.f`

    .. code:: fortran

        double precision	:: pt_chrg_nuclei_interaction	


    repulsion between the point charges and the nuclei

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_pts_charge`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`point_charges`
       * :c:data:`pt_chrg_interaction`
       * :c:data:`pts_charge_coord`
       * :c:data:`pts_charge_z`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`nuclear_repulsion`

 
.. c:var:: double pts_charge_coord


    File : :file:`nuclei/point_charges.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pts_charge_coord	(n_pts_charge,3)


    Coordinates of each point charge.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`n_pts_charge`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_pt_chrg`
       * :c:data:`pt_chrg_interaction`
       * :c:data:`pt_chrg_nuclei_interaction`

 
.. c:var:: double pts_charge_z


    File : :file:`nuclei/point_charges.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pts_charge_z	(n_pts_charge)


    Charge associated to each point charge.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`n_pts_charge`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_pt_chrg`
       * :c:data:`pt_chrg_interaction`
       * :c:data:`pt_chrg_nuclei_interaction`

 
.. c:var:: double slater_bragg_radii


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

 
.. c:var:: double slater_bragg_radii_per_atom


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_radii_per_atom	(nucl_num)


    Slater-Bragg atomic radius in angstroms of each nucleus: for
    i = 1 to nucl_num, the table value slater_bragg_radii(Z_i),
    Z_i being int(nucl_charge(i)) floored at 1.

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

 
.. c:var:: double slater_bragg_radii_per_atom_ua


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_radii_per_atom_ua	(nucl_num)


    As slater_bragg_radii_per_atom, but with the radii in atomic
    units (a0) taken from slater_bragg_radii_ua.

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

 
.. c:var:: double slater_bragg_radii_ua


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_radii_ua	(0:100)


    The values of slater_bragg_radii converted from angstroms to
    atomic units (a0), by a factor of 1.889725989.
    As slater_bragg_radii, indexed by atomic number (1 to 100).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`slater_bragg_radii`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`slater_bragg_radii_per_atom_ua`

 
.. c:var:: double slater_bragg_type_inter_distance


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_type_inter_distance	(nucl_num,nucl_num)


    For i < j, the entry is u / (u*u - 1), where u = (x - 1) /
    (x + 1), x = slater_bragg_radii_per_atom(i) /
    slater_bragg_radii_per_atom(j); other entries are zero.
    Used in cell_function_becke of becke_numerical_grid (Becke
    1988, JCP 88(4)).
    This is the parameter a_ij of Eq. (A5) of A. D. Becke,
    J. Chem. Phys. 88, 2547 (1988): in the coordinate change
    nu_ij = mu_ij + a_ij (1 - mu_ij**2) (Eq. (A2)), it shifts
    the fuzzy cell boundary off the perpendicular bisector, so
    it crosses the internuclear axis where the distance ratio
    to the two nuclei equals the Bragg-Slater radii ratio X
    (Eqs. (A4)-(A6)).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`slater_bragg_radii_per_atom`


 
.. c:var:: double slater_bragg_type_inter_distance_ua


    File : :file:`nuclei/atomic_radii.irp.f`

    .. code:: fortran

        double precision, allocatable	:: slater_bragg_type_inter_distance_ua	(nucl_num,nucl_num)


    As slater_bragg_type_inter_distance, but the ratio x is built
    from the radii in atomic units (a0),
    slater_bragg_radii_per_atom_ua, and each value is clipped to
    the interval [-0.5, 0.5].
    The limit 0.5 is the monotonicity bound |a_ij| < 1/2 of
    Eq. (A3) of A. D. Becke, J. Chem. Phys. 88, 2547 (1988):
    when the computed a_ij falls outside that range, the paper
    prescribes assigning the endpoint value +/- 1/2 (this
    happens when the Bragg-Slater radii differ by a factor
    above about 2.4).

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`slater_bragg_radii_per_atom_ua`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`weight_at_r`
       * :c:data:`weight_at_r_extra`

