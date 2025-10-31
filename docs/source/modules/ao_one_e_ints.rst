.. _module_ao_one_e_ints: 
 
.. program:: ao_one_e_ints 
 
.. default-role:: option 
 
==================
ao_one_e_integrals
==================

All the one-electron integrals in the |AO| basis are here.

The most important providers for usual quantum-chemistry calculation are:

* `ao_kinetic_integrals` which are the kinetic operator integrals on the |AO| basis 
* `ao_integrals_n_e` which are the nuclear-elctron operator integrals on the |AO| basis
* `ao_one_e_integrals` which are the the h_core operator integrals on the |AO| basis


 
 
 
EZFIO parameters 
---------------- 
 
.. option:: ao_integrals_n_e
 
    Nucleus-electron integrals in |AO| basis set
 
 
.. option:: ao_integrals_n_e_imag
 
    Imaginary part of the nucleus-electron integrals in |AO| basis set
 
 
.. option:: io_ao_integrals_n_e
 
    Read/Write |AO| nucleus-electron attraction integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: ao_integrals_kinetic
 
    Kinetic energy integrals in |AO| basis set
 
 
.. option:: ao_integrals_kinetic_imag
 
    Imaginary part of the kinetic energy integrals in |AO| basis set
 
 
.. option:: io_ao_integrals_kinetic
 
    Read/Write |AO| kinetic integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: ao_integrals_pseudo
 
    Pseudopotential integrals in |AO| basis set
 
 
.. option:: ao_integrals_pseudo_imag
 
    Imaginary part of the pseudopotential integrals in |AO| basis set
 
 
.. option:: io_ao_integrals_pseudo
 
    Read/Write |AO| pseudopotential integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: ao_integrals_overlap
 
    Overlap integrals in |AO| basis set
 
 
.. option:: ao_integrals_overlap_imag
 
    Imaginary part of the overlap integrals in |AO| basis set
 
 
.. option:: io_ao_integrals_overlap
 
    Read/Write |AO| overlap integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: ao_one_e_integrals
 
    Combined integrals in |AO| basis set
 
 
.. option:: ao_one_e_integrals_imag
 
    Imaginary part of the combined integrals in |AO| basis set
 
 
.. option:: io_ao_one_e_integrals
 
    Read/Write |AO| one-electron integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: lin_dep_cutoff
 
    Remove linear dependencies when the eigenvalues of the overlap matrix are below this value
 
    Default: 1.e-6
 
.. option:: ao_one_e_integrals_threshold
 
    If | (p|q) | < `ao_one_e_integrals_threshold` then (p|q) is zero
 
    Default: 1.e-15
 
 
Providers 
--------- 
 
.. c:var:: ao_cart_to_sphe_coef


    File : :file:`ao_one_e_ints/ao_ortho_canonical.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_cart_to_sphe_coef	(ao_num,ao_num)
        integer	:: ao_cart_to_sphe_num	


    Coefficients to go from cartesian to spherical coordinates in the current
    basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`cart_to_sphe_1`
       * :c:data:`cart_to_sphe_2`
       * :c:data:`cart_to_sphe_3`
       * :c:data:`cart_to_sphe_4`
       * :c:data:`cart_to_sphe_5`
       * :c:data:`cart_to_sphe_6`
       * :c:data:`cart_to_sphe_7`
       * :c:data:`cart_to_sphe_8`
       * :c:data:`cart_to_sphe_9`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_inv`
       * :c:data:`ao_cart_to_sphe_overlap`
       * :c:data:`ao_ortho_canonical_coef`

 
.. c:var:: ao_cart_to_sphe_inv


    File : :file:`ao_one_e_ints/ao_ortho_canonical.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_cart_to_sphe_inv	(ao_cart_to_sphe_num,ao_num)


    Inverse of :c:data:`ao_cart_to_sphe_coef`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`
       * :c:data:`ao_num`
       * :c:data:`lin_dep_cutoff`


 
.. c:var:: ao_cart_to_sphe_num


    File : :file:`ao_one_e_ints/ao_ortho_canonical.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_cart_to_sphe_coef	(ao_num,ao_num)
        integer	:: ao_cart_to_sphe_num	


    Coefficients to go from cartesian to spherical coordinates in the current
    basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_l`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`cart_to_sphe_1`
       * :c:data:`cart_to_sphe_2`
       * :c:data:`cart_to_sphe_3`
       * :c:data:`cart_to_sphe_4`
       * :c:data:`cart_to_sphe_5`
       * :c:data:`cart_to_sphe_6`
       * :c:data:`cart_to_sphe_7`
       * :c:data:`cart_to_sphe_8`
       * :c:data:`cart_to_sphe_9`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_inv`
       * :c:data:`ao_cart_to_sphe_overlap`
       * :c:data:`ao_ortho_canonical_coef`

 
.. c:var:: ao_cart_to_sphe_overlap


    File : :file:`ao_one_e_ints/ao_ortho_canonical.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_cart_to_sphe_overlap	(ao_cart_to_sphe_num,ao_cart_to_sphe_num)


    |AO| overlap matrix in the spherical basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_canonical_coef`

 
.. c:var:: ao_coef_cgtos_norm_ord_transp


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_cgtos_norm_ord_transp	(ao_prim_num_max,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos_ord`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_2e_cgtos_schwartz`
       * :c:data:`ao_deriv2_cgtos_x`
       * :c:data:`ao_integrals_n_e_cgtos`
       * :c:data:`ao_overlap_cgtos`

 
.. c:var:: ao_coef_norm_cgtos


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_norm_cgtos	(ao_num,ao_prim_num_max)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_expo_im`
       * :c:data:`ao_expo_phase`
       * :c:data:`ao_expo_pw`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`nucl_coord`
       * :c:data:`primitives_normalized`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos_ord`

 
.. c:var:: ao_coef_norm_cgtos_ord


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_norm_cgtos_ord	(ao_num,ao_prim_num_max)
        complex*16, allocatable	:: ao_expo_cgtos_ord	(ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_expo_pw_ord	(4,ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_expo_phase_ord	(4,ao_num,ao_prim_num_max)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos`
       * :c:data:`ao_expo`
       * :c:data:`ao_expo_im`
       * :c:data:`ao_expo_phase`
       * :c:data:`ao_expo_pw`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`

 
.. c:var:: ao_coef_norm_cosgtos


    File : :file:`ao_one_e_ints/aos_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_norm_cosgtos	(ao_num,ao_prim_num_max)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef`
       * :c:data:`ao_expo`
       * :c:data:`ao_expoim_cosgtos`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`
       * :c:data:`primitives_normalized`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_cosgtos`

 
.. c:var:: ao_coef_norm_ord_cosgtos


    File : :file:`ao_one_e_ints/aos_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_norm_ord_cosgtos	(ao_num,ao_prim_num_max)
        complex*16, allocatable	:: ao_expo_ord_cosgtos	(ao_num,ao_prim_num_max)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cosgtos`
       * :c:data:`ao_expo`
       * :c:data:`ao_expoim_cosgtos`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`

 
.. c:var:: ao_coef_norm_ord_transp_cosgtos


    File : :file:`ao_one_e_ints/aos_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_norm_ord_transp_cosgtos	(ao_prim_num_max,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_cosgtos`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_2e_cosgtos_schwartz`
       * :c:data:`ao_deriv2_cosgtos_x`
       * :c:data:`ao_integrals_n_e_cosgtos`
       * :c:data:`ao_overlap_cosgtos`

 
.. c:var:: ao_deriv2_cgtos_x


    File : :file:`ao_one_e_ints/one_e_kin_integrals_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv2_cgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_kinetic_integrals_cgtos`

 
.. c:var:: ao_deriv2_cgtos_y


    File : :file:`ao_one_e_ints/one_e_kin_integrals_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv2_cgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_kinetic_integrals_cgtos`

 
.. c:var:: ao_deriv2_cgtos_z


    File : :file:`ao_one_e_ints/one_e_kin_integrals_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv2_cgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_kinetic_integrals_cgtos`

 
.. c:var:: ao_deriv2_cosgtos_x


    File : :file:`ao_one_e_ints/one_e_kin_integrals_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv2_cosgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cosgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cosgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_kinetic_integrals_cosgtos`

 
.. c:var:: ao_deriv2_cosgtos_y


    File : :file:`ao_one_e_ints/one_e_kin_integrals_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv2_cosgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cosgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cosgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_kinetic_integrals_cosgtos`

 
.. c:var:: ao_deriv2_cosgtos_z


    File : :file:`ao_one_e_ints/one_e_kin_integrals_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv2_cosgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cosgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_cosgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_kinetic_integrals_cosgtos`

 
.. c:var:: ao_deriv2_x


    File : :file:`ao_one_e_ints/kin_ao_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv2_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_z	(ao_num,ao_num)


    Second derivative matrix elements in the |AO| basis.
    
    .. math::
    
      {\tt ao\_deriv2\_x} =
      \langle \chi_i(x,y,z) | \frac{\partial^2}{\partial x^2} |\chi_j (x,y,z) \rangle
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_deriv2_cgtos_x`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`
       * :c:data:`use_cgtos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_kinetic_integrals`

 
.. c:var:: ao_deriv2_y


    File : :file:`ao_one_e_ints/kin_ao_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv2_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_z	(ao_num,ao_num)


    Second derivative matrix elements in the |AO| basis.
    
    .. math::
    
      {\tt ao\_deriv2\_x} =
      \langle \chi_i(x,y,z) | \frac{\partial^2}{\partial x^2} |\chi_j (x,y,z) \rangle
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_deriv2_cgtos_x`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`
       * :c:data:`use_cgtos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_kinetic_integrals`

 
.. c:var:: ao_deriv2_z


    File : :file:`ao_one_e_ints/kin_ao_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv2_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv2_z	(ao_num,ao_num)


    Second derivative matrix elements in the |AO| basis.
    
    .. math::
    
      {\tt ao\_deriv2\_x} =
      \langle \chi_i(x,y,z) | \frac{\partial^2}{\partial x^2} |\chi_j (x,y,z) \rangle
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_deriv2_cgtos_x`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`
       * :c:data:`use_cgtos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_kinetic_integrals`

 
.. c:var:: ao_deriv_1_x


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv_1_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv_1_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv_1_z	(ao_num,ao_num)


    * array of the integrals of AO_i * d/dx  AO_j
    
    * array of the integrals of AO_i * d/dy  AO_j
    
    * array of the integrals of AO_i * d/dz  AO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_deriv_1_x`

 
.. c:var:: ao_deriv_1_y


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv_1_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv_1_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv_1_z	(ao_num,ao_num)


    * array of the integrals of AO_i * d/dx  AO_j
    
    * array of the integrals of AO_i * d/dy  AO_j
    
    * array of the integrals of AO_i * d/dz  AO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_deriv_1_x`

 
.. c:var:: ao_deriv_1_z


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_deriv_1_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv_1_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_deriv_1_z	(ao_num,ao_num)


    * array of the integrals of AO_i * d/dx  AO_j
    
    * array of the integrals of AO_i * d/dy  AO_j
    
    * array of the integrals of AO_i * d/dz  AO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_deriv_1_x`

 
.. c:var:: ao_dipole_x


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_dipole_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_dipole_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_dipole_z	(ao_num,ao_num)


    * array of the integrals of AO_i * x AO_j
    
    * array of the integrals of AO_i * y AO_j
    
    * array of the integrals of AO_i * z AO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_dipole_x`

 
.. c:var:: ao_dipole_y


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_dipole_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_dipole_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_dipole_z	(ao_num,ao_num)


    * array of the integrals of AO_i * x AO_j
    
    * array of the integrals of AO_i * y AO_j
    
    * array of the integrals of AO_i * z AO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_dipole_x`

 
.. c:var:: ao_dipole_z


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_dipole_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_dipole_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_dipole_z	(ao_num,ao_num)


    * array of the integrals of AO_i * x AO_j
    
    * array of the integrals of AO_i * y AO_j
    
    * array of the integrals of AO_i * z AO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_dipole_x`

 
.. c:var:: ao_expo_cgtos_ord


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_norm_cgtos_ord	(ao_num,ao_prim_num_max)
        complex*16, allocatable	:: ao_expo_cgtos_ord	(ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_expo_pw_ord	(4,ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_expo_phase_ord	(4,ao_num,ao_prim_num_max)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos`
       * :c:data:`ao_expo`
       * :c:data:`ao_expo_im`
       * :c:data:`ao_expo_phase`
       * :c:data:`ao_expo_pw`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`

 
.. c:var:: ao_expo_cgtos_ord_transp


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        complex*16, allocatable	:: ao_expo_cgtos_ord_transp	(ao_prim_num_max,ao_num)
        double precision, allocatable	:: ao_expo_pw_ord_transp	(4,ao_prim_num_max,ao_num)
        double precision, allocatable	:: ao_expo_phase_ord_transp	(4,ao_prim_num_max,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos_ord`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_2e_cgtos_schwartz`
       * :c:data:`ao_deriv2_cgtos_x`
       * :c:data:`ao_integrals_n_e_cgtos`
       * :c:data:`ao_overlap_cgtos`
       * :c:data:`use_pw`

 
.. c:var:: ao_expo_ord_cosgtos


    File : :file:`ao_one_e_ints/aos_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_norm_ord_cosgtos	(ao_num,ao_prim_num_max)
        complex*16, allocatable	:: ao_expo_ord_cosgtos	(ao_num,ao_prim_num_max)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cosgtos`
       * :c:data:`ao_expo`
       * :c:data:`ao_expoim_cosgtos`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`

 
.. c:var:: ao_expo_ord_transp_cosgtos


    File : :file:`ao_one_e_ints/aos_cosgtos.irp.f`

    .. code:: fortran

        complex*16, allocatable	:: ao_expo_ord_transp_cosgtos	(ao_prim_num_max,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_cosgtos`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_2e_cosgtos_schwartz`
       * :c:data:`ao_deriv2_cosgtos_x`
       * :c:data:`ao_integrals_n_e_cosgtos`
       * :c:data:`ao_overlap_cosgtos`

 
.. c:var:: ao_expo_phase_ord


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_norm_cgtos_ord	(ao_num,ao_prim_num_max)
        complex*16, allocatable	:: ao_expo_cgtos_ord	(ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_expo_pw_ord	(4,ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_expo_phase_ord	(4,ao_num,ao_prim_num_max)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos`
       * :c:data:`ao_expo`
       * :c:data:`ao_expo_im`
       * :c:data:`ao_expo_phase`
       * :c:data:`ao_expo_pw`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`

 
.. c:var:: ao_expo_phase_ord_transp


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        complex*16, allocatable	:: ao_expo_cgtos_ord_transp	(ao_prim_num_max,ao_num)
        double precision, allocatable	:: ao_expo_pw_ord_transp	(4,ao_prim_num_max,ao_num)
        double precision, allocatable	:: ao_expo_phase_ord_transp	(4,ao_prim_num_max,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos_ord`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_2e_cgtos_schwartz`
       * :c:data:`ao_deriv2_cgtos_x`
       * :c:data:`ao_integrals_n_e_cgtos`
       * :c:data:`ao_overlap_cgtos`
       * :c:data:`use_pw`

 
.. c:var:: ao_expo_pw_ord


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_coef_norm_cgtos_ord	(ao_num,ao_prim_num_max)
        complex*16, allocatable	:: ao_expo_cgtos_ord	(ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_expo_pw_ord	(4,ao_num,ao_prim_num_max)
        double precision, allocatable	:: ao_expo_phase_ord	(4,ao_num,ao_prim_num_max)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos`
       * :c:data:`ao_expo`
       * :c:data:`ao_expo_im`
       * :c:data:`ao_expo_phase`
       * :c:data:`ao_expo_pw`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`

 
.. c:var:: ao_expo_pw_ord_transp


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        complex*16, allocatable	:: ao_expo_cgtos_ord_transp	(ao_prim_num_max,ao_num)
        double precision, allocatable	:: ao_expo_pw_ord_transp	(4,ao_prim_num_max,ao_num)
        double precision, allocatable	:: ao_expo_phase_ord_transp	(4,ao_prim_num_max,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_cgtos_ord`
       * :c:data:`ao_num`
       * :c:data:`ao_prim_num_max`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_2e_cgtos_schwartz`
       * :c:data:`ao_deriv2_cgtos_x`
       * :c:data:`ao_integrals_n_e_cgtos`
       * :c:data:`ao_overlap_cgtos`
       * :c:data:`use_pw`

 
.. c:var:: ao_integrals_n_e


    File : :file:`ao_one_e_ints/pot_ao_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_integrals_n_e	(ao_num,ao_num)


    Nucleus-electron interaction, in the |AO| basis set.
    
    :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
    
    These integrals also contain the pseudopotential integrals.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_integrals_n_e_cgtos`
       * :c:data:`ao_integrals_pt_chrg`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`ao_pseudo_integrals`
       * :c:data:`do_pseudo`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`point_charges`
       * :c:data:`read_ao_integrals_n_e`
       * :c:data:`use_cgtos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_ortho_canonical_nucl_elec_integrals`
       * :c:data:`ao_ortho_lowdin_nucl_elec_integrals`
       * :c:data:`hf_kinetic_energy`
       * :c:data:`mo_integrals_n_e`

 
.. c:var:: ao_integrals_n_e_cgtos


    File : :file:`ao_one_e_ints/one_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_integrals_n_e_cgtos	(ao_num,ao_num)


    
    Nucleus-electron interaction, in the cgtos |AO| basis set.
    
    :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`use_pw`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`

 
.. c:var:: ao_integrals_n_e_cosgtos


    File : :file:`ao_one_e_ints/one_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_integrals_n_e_cosgtos	(ao_num,ao_num)


    
    Nucleus-electron interaction, in the cosgtos |AO| basis set.
    
    :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`

 
.. c:var:: ao_integrals_n_e_imag


    File : :file:`ao_one_e_ints/pot_ao_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_integrals_n_e_imag	(ao_num,ao_num)


    Nucleus-electron interaction, in the |AO| basis set.
    
    :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`read_ao_integrals_n_e`


 
.. c:var:: ao_integrals_n_e_per_atom


    File : :file:`ao_one_e_ints/pot_ao_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_integrals_n_e_per_atom	(ao_num,ao_num,nucl_num)


    Nucleus-electron interaction in the |AO| basis set, per atom A.
    
    :math:`\langle \chi_i | -\frac{1}{|r-R_A|} | \chi_j \rangle`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_n_e_per_atom`

 
.. c:var:: ao_integrals_pt_chrg


    File : :file:`ao_one_e_ints/pot_pt_charges.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_integrals_pt_chrg	(ao_num,ao_num)


     Point charge-electron interaction, in the |AO| basis set.
    
     :math:`\langle \chi_i | -\sum_charge charge * \frac{1}{|r-R_charge|} | \chi_j \rangle`
    
    Notice the minus sign convention as it is supposed to be for electrons.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`n_pts_charge`
       * :c:data:`nucl_coord`
       * :c:data:`pts_charge_coord`
       * :c:data:`pts_charge_z`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`

 
.. c:var:: ao_kinetic_integrals


    File : :file:`ao_one_e_ints/kin_ao_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_kinetic_integrals	(ao_num,ao_num)


    Kinetic energy integrals in the |AO| basis.
    
    :math:`\langle \chi_i |\hat{T}| \chi_j \rangle` 
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_x`
       * :c:data:`ao_num`
       * :c:data:`read_ao_integrals_kinetic`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_one_e_integrals`
       * :c:data:`hf_kinetic_energy`
       * :c:data:`mo_kinetic_integrals`

 
.. c:var:: ao_kinetic_integrals_cgtos


    File : :file:`ao_one_e_ints/one_e_kin_integrals_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_kinetic_integrals_cgtos	(ao_num,ao_num)


    
    Kinetic energy integrals in the cgtos |AO| basis.
    
    :math:`\langle \chi_i |\hat{T}| \chi_j \rangle` 
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_cgtos_x`
       * :c:data:`ao_num`


 
.. c:var:: ao_kinetic_integrals_cosgtos


    File : :file:`ao_one_e_ints/one_e_kin_integrals_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_kinetic_integrals_cosgtos	(ao_num,ao_num)


    
    Kinetic energy integrals in the cosgtos |AO| basis.
    
    :math:`\langle \chi_i |\hat{T}| \chi_j \rangle` 
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv2_cosgtos_x`
       * :c:data:`ao_num`


 
.. c:var:: ao_kinetic_integrals_imag


    File : :file:`ao_one_e_ints/kin_ao_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_kinetic_integrals_imag	(ao_num,ao_num)


    Kinetic energy integrals in the |AO| basis.
    
    :math:`\langle \chi_i |\hat{T}| \chi_j \rangle` 
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`read_ao_integrals_kinetic`


 
.. c:var:: ao_one_e_integrals


    File : :file:`ao_one_e_ints/ao_one_e_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_one_e_integrals	(ao_num,ao_num)
        double precision, allocatable	:: ao_one_e_integrals_diag	(ao_num)


    One-electron Hamiltonian in the |AO| basis.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_kinetic_integrals`
       * :c:data:`ao_num`
       * :c:data:`read_ao_one_e_integrals`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`hf_energy`
       * :c:data:`scf_energy`

 
.. c:var:: ao_one_e_integrals_diag


    File : :file:`ao_one_e_ints/ao_one_e_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_one_e_integrals	(ao_num,ao_num)
        double precision, allocatable	:: ao_one_e_integrals_diag	(ao_num)


    One-electron Hamiltonian in the |AO| basis.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_kinetic_integrals`
       * :c:data:`ao_num`
       * :c:data:`read_ao_one_e_integrals`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`hf_energy`
       * :c:data:`scf_energy`

 
.. c:var:: ao_one_e_integrals_imag


    File : :file:`ao_one_e_ints/ao_one_e_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_one_e_integrals_imag	(ao_num,ao_num)


    One-electron Hamiltonian in the |AO| basis.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`read_ao_one_e_integrals`


 
.. c:var:: ao_ortho_canonical_coef


    File : :file:`ao_one_e_ints/ao_ortho_canonical.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_ortho_canonical_coef	(ao_num,ao_num)
        integer	:: ao_ortho_canonical_num	


    matrix of the coefficients of the mos generated by the
    orthonormalization by the S^{-1/2} canonical transformation of the aos
    ao_ortho_canonical_coef(i,j) = coefficient of the ith ao on the jth ao_ortho_canonical orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`
       * :c:data:`ao_cart_to_sphe_overlap`
       * :c:data:`ao_cartesian`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`lin_dep_cutoff`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_canonical_coef_inv`
       * :c:data:`ao_ortho_canonical_nucl_elec_integrals`
       * :c:data:`ao_ortho_canonical_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_aux`
       * :c:data:`mo_num`

 
.. c:var:: ao_ortho_canonical_coef_inv


    File : :file:`ao_one_e_ints/ao_ortho_canonical.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_ortho_canonical_coef_inv	(ao_num,ao_num)


    ao_ortho_canonical_coef^(-1)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_ortho_canonical_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_coef_in_ao_ortho_basis`

 
.. c:var:: ao_ortho_canonical_num


    File : :file:`ao_one_e_ints/ao_ortho_canonical.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_ortho_canonical_coef	(ao_num,ao_num)
        integer	:: ao_ortho_canonical_num	


    matrix of the coefficients of the mos generated by the
    orthonormalization by the S^{-1/2} canonical transformation of the aos
    ao_ortho_canonical_coef(i,j) = coefficient of the ith ao on the jth ao_ortho_canonical orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_coef`
       * :c:data:`ao_cart_to_sphe_overlap`
       * :c:data:`ao_cartesian`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`lin_dep_cutoff`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_ortho_canonical_coef_inv`
       * :c:data:`ao_ortho_canonical_nucl_elec_integrals`
       * :c:data:`ao_ortho_canonical_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_aux`
       * :c:data:`mo_num`

 
.. c:var:: ao_ortho_canonical_overlap


    File : :file:`ao_one_e_ints/ao_ortho_canonical.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_ortho_canonical_overlap	(ao_ortho_canonical_num,ao_ortho_canonical_num)


    overlap matrix of the ao_ortho_canonical.
    Expected to be the Identity

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ao_overlap`


 
.. c:var:: ao_overlap


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_z	(ao_num,ao_num)


    Overlap between atomic basis functions:
    
    :math:`\int \chi_i(r) \chi_j(r) dr`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_cgtos`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`
       * :c:data:`read_ao_integrals_overlap`
       * :c:data:`use_cgtos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_overlap`
       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ao_ortho_canonical_overlap`
       * :c:data:`ao_ortho_lowdin_coef`
       * :c:data:`ao_ortho_lowdin_overlap`
       * :c:data:`ao_overlap_complex`
       * :c:data:`fps_spf_matrix_ao`
       * :c:data:`mo_overlap`
       * :c:data:`s_half`
       * :c:data:`s_half_inv`
       * :c:data:`s_inv`
       * :c:data:`s_mo_coef`

 
.. c:var:: ao_overlap_abs


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_abs	(ao_num,ao_num)


    Overlap between absolute values of atomic basis functions:
    
    :math:`\int |\chi_i(r)| |\chi_j(r)| dr`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_complex`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`is_periodic`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`cholesky_ao_num`

 
.. c:var:: ao_overlap_cgtos


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_cgtos	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap`

 
.. c:var:: ao_overlap_cgtos_x


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_cgtos	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap`

 
.. c:var:: ao_overlap_cgtos_y


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_cgtos	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap`

 
.. c:var:: ao_overlap_cgtos_z


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_cgtos	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_cgtos_norm_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap`

 
.. c:var:: ao_overlap_complex


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        complex*16, allocatable	:: ao_overlap_complex	(ao_num,ao_num)


    Overlap for complex AOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`ao_overlap_imag`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap_abs`
       * :c:data:`s_inv_complex`

 
.. c:var:: ao_overlap_cosgtos


    File : :file:`ao_one_e_ints/aos_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_cosgtos	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap`

 
.. c:var:: ao_overlap_cosgtos_x


    File : :file:`ao_one_e_ints/aos_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_cosgtos	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap`

 
.. c:var:: ao_overlap_cosgtos_y


    File : :file:`ao_one_e_ints/aos_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_cosgtos	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap`

 
.. c:var:: ao_overlap_cosgtos_z


    File : :file:`ao_one_e_ints/aos_cosgtos.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_cosgtos	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_cosgtos_z	(ao_num,ao_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_norm_ord_transp_cosgtos`
       * :c:data:`ao_expo_ord_transp_cosgtos`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap`

 
.. c:var:: ao_overlap_imag


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap_imag	(ao_num,ao_num)


    Imaginary part of the overlap

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap_complex`

 
.. c:var:: ao_overlap_x


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_z	(ao_num,ao_num)


    Overlap between atomic basis functions:
    
    :math:`\int \chi_i(r) \chi_j(r) dr`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_cgtos`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`
       * :c:data:`read_ao_integrals_overlap`
       * :c:data:`use_cgtos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_overlap`
       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ao_ortho_canonical_overlap`
       * :c:data:`ao_ortho_lowdin_coef`
       * :c:data:`ao_ortho_lowdin_overlap`
       * :c:data:`ao_overlap_complex`
       * :c:data:`fps_spf_matrix_ao`
       * :c:data:`mo_overlap`
       * :c:data:`s_half`
       * :c:data:`s_half_inv`
       * :c:data:`s_inv`
       * :c:data:`s_mo_coef`

 
.. c:var:: ao_overlap_y


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_z	(ao_num,ao_num)


    Overlap between atomic basis functions:
    
    :math:`\int \chi_i(r) \chi_j(r) dr`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_cgtos`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`
       * :c:data:`read_ao_integrals_overlap`
       * :c:data:`use_cgtos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_overlap`
       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ao_ortho_canonical_overlap`
       * :c:data:`ao_ortho_lowdin_coef`
       * :c:data:`ao_ortho_lowdin_overlap`
       * :c:data:`ao_overlap_complex`
       * :c:data:`fps_spf_matrix_ao`
       * :c:data:`mo_overlap`
       * :c:data:`s_half`
       * :c:data:`s_half_inv`
       * :c:data:`s_inv`
       * :c:data:`s_mo_coef`

 
.. c:var:: ao_overlap_z


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_overlap	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_overlap_z	(ao_num,ao_num)


    Overlap between atomic basis functions:
    
    :math:`\int \chi_i(r) \chi_j(r) dr`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_overlap_cgtos`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`
       * :c:data:`read_ao_integrals_overlap`
       * :c:data:`use_cgtos`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_cart_to_sphe_overlap`
       * :c:data:`ao_ortho_canonical_coef`
       * :c:data:`ao_ortho_canonical_overlap`
       * :c:data:`ao_ortho_lowdin_coef`
       * :c:data:`ao_ortho_lowdin_overlap`
       * :c:data:`ao_overlap_complex`
       * :c:data:`fps_spf_matrix_ao`
       * :c:data:`mo_overlap`
       * :c:data:`s_half`
       * :c:data:`s_half_inv`
       * :c:data:`s_inv`
       * :c:data:`s_mo_coef`

 
.. c:var:: ao_pseudo_integrals


    File : :file:`ao_one_e_ints/pot_ao_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_pseudo_integrals	(ao_num,ao_num)


    Pseudo-potential integrals in the |AO| basis set.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_pseudo_integrals_local`
       * :c:data:`ao_pseudo_integrals_non_local`
       * :c:data:`do_pseudo`
       * :c:data:`pseudo_klocmax`
       * :c:data:`pseudo_kmax`
       * :c:data:`read_ao_integrals_pseudo`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`
       * :c:data:`mo_pseudo_integrals`

 
.. c:var:: ao_pseudo_integrals_local


    File : :file:`ao_one_e_ints/pot_ao_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_pseudo_integrals_local	(ao_num,ao_num)


    Local pseudo-potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pseudo_v_k_transp`
       * :c:data:`pseudo_klocmax`
       * :c:data:`pseudo_v_k_transp`
       * :c:data:`pseudo_v_k_transp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_pseudo_integrals`
       * :c:data:`mo_pseudo_integrals_local`

 
.. c:var:: ao_pseudo_integrals_non_local


    File : :file:`ao_one_e_ints/pot_ao_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_pseudo_integrals_non_local	(ao_num,ao_num)


    Non-local pseudo-potential

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`pseudo_v_kl_transp`
       * :c:data:`pseudo_kmax`
       * :c:data:`pseudo_lmax`
       * :c:data:`pseudo_v_kl_transp`
       * :c:data:`pseudo_v_kl_transp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_pseudo_integrals`
       * :c:data:`mo_pseudo_integrals_non_local`

 
.. c:var:: ao_spread_x


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_spread_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_spread_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_spread_z	(ao_num,ao_num)


    * array of the integrals of AO_i * x^2 AO_j
    
    * array of the integrals of AO_i * y^2 AO_j
    
    * array of the integrals of AO_i * z^2 AO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_spread_x`

 
.. c:var:: ao_spread_y


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_spread_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_spread_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_spread_z	(ao_num,ao_num)


    * array of the integrals of AO_i * x^2 AO_j
    
    * array of the integrals of AO_i * y^2 AO_j
    
    * array of the integrals of AO_i * z^2 AO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_spread_x`

 
.. c:var:: ao_spread_z


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_spread_x	(ao_num,ao_num)
        double precision, allocatable	:: ao_spread_y	(ao_num,ao_num)
        double precision, allocatable	:: ao_spread_z	(ao_num,ao_num)


    * array of the integrals of AO_i * x^2 AO_j
    
    * array of the integrals of AO_i * y^2 AO_j
    
    * array of the integrals of AO_i * z^2 AO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_num`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`nucl_coord`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_spread_x`

 
.. c:function:: give_cpolynomial_mult_center_one_e:


    File : :file:`ao_one_e_ints/one_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        subroutine give_cpolynomial_mult_center_one_e(A_center, B_center, alpha, beta, &
                                              power_A, power_B, C_center, n_pt_in, d, n_pt_out)


    Returns the explicit polynomial in terms of the "t" variable of the following
    
    $I_{x1}(a_x, d_x,p,q) \times I_{x1}(a_y, d_y,p,q) \times I_{x1}(a_z, d_z,p,q)$.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`nai_pol_mult_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e_cgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x1_pol_mult_one_e:


    File : :file:`ao_one_e_ints/pot_ao_ints.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_one_e(a,c,R1x,R1xp,R2x,d,nd,n_pt_in)


    Recursive routine involved in the electron-nucleus potential

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_polynomial_mult_center_one_e`
       * :c:func:`give_polynomial_mult_center_one_e_erf`
       * :c:func:`give_polynomial_mult_center_one_e_erf_opt`
       * :c:func:`i_x1_pol_mult_one_e`
       * :c:func:`i_x2_pol_mult_one_e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e`
       * :c:func:`i_x2_pol_mult_one_e`
       * :c:func:`multiply_poly_c2`

 
.. c:function:: i_x1_pol_mult_one_e_cgtos:


    File : :file:`ao_one_e_ints/one_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_one_e_cgtos(a, c, R1x, R1xp, R2x, d, nd, n_pt_in)


    Recursive routine involved in the electron-nucleus potential

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_cpolynomial_mult_center_one_e`
       * :c:func:`i_x1_pol_mult_one_e_cgtos`
       * :c:func:`i_x2_pol_mult_one_e_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e_cgtos`
       * :c:func:`i_x2_pol_mult_one_e_cgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x1_pol_mult_one_e_cosgtos:


    File : :file:`ao_one_e_ints/one_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x1_pol_mult_one_e_cosgtos(a, c, R1x, R1xp, R2x, d, nd, n_pt_in)


    Recursive routine involved in the electron-nucleus potential

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`give_cpolynomial_mult_center_one_e`
       * :c:func:`i_x1_pol_mult_one_e_cosgtos`
       * :c:func:`i_x2_pol_mult_one_e_cosgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e_cosgtos`
       * :c:func:`i_x2_pol_mult_one_e_cosgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x2_pol_mult_one_e:


    File : :file:`ao_one_e_ints/pot_ao_ints.irp.f`

    .. code:: fortran

        recursive subroutine I_x2_pol_mult_one_e(c,R1x,R1xp,R2x,d,nd,dim)


    Recursive routine involved in the electron-nucleus potential

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e`
       * :c:func:`multiply_poly_c2`

 
.. c:function:: i_x2_pol_mult_one_e_cgtos:


    File : :file:`ao_one_e_ints/one_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x2_pol_mult_one_e_cgtos(c, R1x, R1xp, R2x, d, nd, dim)


    Recursive routine involved in the electron-nucleus potential

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e_cgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e_cgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: i_x2_pol_mult_one_e_cosgtos:


    File : :file:`ao_one_e_ints/one_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        recursive subroutine I_x2_pol_mult_one_e_cosgtos(c, R1x, R1xp, R2x, d, nd, dim)


    Recursive routine involved in the electron-nucleus potential

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e_cosgtos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e_cosgtos`
       * :c:func:`multiply_cpoly`

 
.. c:function:: nai_pol_mult_cgtos:


    File : :file:`ao_one_e_ints/one_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        complex*16 function NAI_pol_mult_cgtos(Ae_center, Be_center, power_A, power_B, alpha, beta, &
                                       Ap_center, Bp_center, C_center, n_pt_in)


    
    Computes the electron-nucleus attraction with two primitves cgtos.
    
    :math:`\langle g_i | \frac{1}{|r-R_c|} | g_j \rangle`
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`use_pw`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_cpolynomial_mult_center_one_e`

 
.. c:function:: nai_pol_mult_erf_with1s:


    File : :file:`ao_one_e_ints/pot_ao_erf_ints.irp.f`

    .. code:: fortran

        double precision function NAI_pol_mult_erf_with1s( A1_center, A2_center, power_A1, power_A2, alpha1, alpha2 &
                                                 , beta, B_center, C_center, n_pt_in, mu_in )


    
    Computes the following integral :
    
    .. math::
    
      \int dx (x - A1_x)^a_1 (x - B1_x)^a_2 \exp(-\alpha_1 (x - A1_x)^2 - \alpha_2 (x - A2_x)^2)
      \int dy (y - A1_y)^b_1 (y - B1_y)^b_2 \exp(-\alpha_1 (y - A1_y)^2 - \alpha_2 (y - A2_y)^2)
      \int dz (x - A1_z)^c_1 (z - B1_z)^c_2 \exp(-\alpha_1 (z - A1_z)^2 - \alpha_2 (z - A2_z)^2)
      \exp(-\beta (r - B)^2)
      \frac{\erf(\mu |r - R_C|)}{|r - R_C|}$.
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_polynomial_mult_center_one_e_erf_opt`

 
.. c:var:: pseudo_dz_k_transp


    File : :file:`ao_one_e_ints/pot_ao_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pseudo_v_k_transp	(pseudo_klocmax,nucl_num)
        integer, allocatable	:: pseudo_n_k_transp	(pseudo_klocmax,nucl_num)
        double precision, allocatable	:: pseudo_dz_k_transp	(pseudo_klocmax,nucl_num)


    Transposed arrays for pseudopotentials

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pseudo_dz_k`
       * :c:data:`pseudo_klocmax`
       * :c:data:`pseudo_n_k`
       * :c:data:`pseudo_v_k`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_pseudo_integrals_local`

 
.. c:var:: pseudo_dz_kl_transp


    File : :file:`ao_one_e_ints/pot_ao_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pseudo_v_kl_transp	(pseudo_kmax,0:pseudo_lmax,nucl_num)
        integer, allocatable	:: pseudo_n_kl_transp	(pseudo_kmax,0:pseudo_lmax,nucl_num)
        double precision, allocatable	:: pseudo_dz_kl_transp	(pseudo_kmax,0:pseudo_lmax,nucl_num)


    Transposed arrays for pseudopotentials

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pseudo_dz_kl`
       * :c:data:`pseudo_kmax`
       * :c:data:`pseudo_lmax`
       * :c:data:`pseudo_n_kl`
       * :c:data:`pseudo_v_kl`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_pseudo_integrals_non_local`

 
.. c:var:: pseudo_n_k_transp


    File : :file:`ao_one_e_ints/pot_ao_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pseudo_v_k_transp	(pseudo_klocmax,nucl_num)
        integer, allocatable	:: pseudo_n_k_transp	(pseudo_klocmax,nucl_num)
        double precision, allocatable	:: pseudo_dz_k_transp	(pseudo_klocmax,nucl_num)


    Transposed arrays for pseudopotentials

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pseudo_dz_k`
       * :c:data:`pseudo_klocmax`
       * :c:data:`pseudo_n_k`
       * :c:data:`pseudo_v_k`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_pseudo_integrals_local`

 
.. c:var:: pseudo_n_kl_transp


    File : :file:`ao_one_e_ints/pot_ao_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pseudo_v_kl_transp	(pseudo_kmax,0:pseudo_lmax,nucl_num)
        integer, allocatable	:: pseudo_n_kl_transp	(pseudo_kmax,0:pseudo_lmax,nucl_num)
        double precision, allocatable	:: pseudo_dz_kl_transp	(pseudo_kmax,0:pseudo_lmax,nucl_num)


    Transposed arrays for pseudopotentials

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pseudo_dz_kl`
       * :c:data:`pseudo_kmax`
       * :c:data:`pseudo_lmax`
       * :c:data:`pseudo_n_kl`
       * :c:data:`pseudo_v_kl`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_pseudo_integrals_non_local`

 
.. c:var:: pseudo_v_k_transp


    File : :file:`ao_one_e_ints/pot_ao_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pseudo_v_k_transp	(pseudo_klocmax,nucl_num)
        integer, allocatable	:: pseudo_n_k_transp	(pseudo_klocmax,nucl_num)
        double precision, allocatable	:: pseudo_dz_k_transp	(pseudo_klocmax,nucl_num)


    Transposed arrays for pseudopotentials

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pseudo_dz_k`
       * :c:data:`pseudo_klocmax`
       * :c:data:`pseudo_n_k`
       * :c:data:`pseudo_v_k`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_pseudo_integrals_local`

 
.. c:var:: pseudo_v_kl_transp


    File : :file:`ao_one_e_ints/pot_ao_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: pseudo_v_kl_transp	(pseudo_kmax,0:pseudo_lmax,nucl_num)
        integer, allocatable	:: pseudo_n_kl_transp	(pseudo_kmax,0:pseudo_lmax,nucl_num)
        double precision, allocatable	:: pseudo_dz_kl_transp	(pseudo_kmax,0:pseudo_lmax,nucl_num)


    Transposed arrays for pseudopotentials

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`nucl_num`
       * :c:data:`pseudo_dz_kl`
       * :c:data:`pseudo_kmax`
       * :c:data:`pseudo_lmax`
       * :c:data:`pseudo_n_kl`
       * :c:data:`pseudo_v_kl`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_pseudo_integrals_non_local`

 
.. c:var:: s_half


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: s_half	(ao_num,ao_num)


    :math:`S^{1/2}`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`


 
.. c:var:: s_half_inv


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: s_half_inv	(AO_num,AO_num)


    :math:`X = S^{-1/2}` obtained by SVD

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigenvalues_fock_matrix_ao`

 
.. c:var:: s_inv


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: s_inv	(ao_num,ao_num)


    Inverse of the overlap matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`lin_dep_cutoff`


 
.. c:var:: s_inv_complex


    File : :file:`ao_one_e_ints/ao_overlap.irp.f`

    .. code:: fortran

        complex*16, allocatable	:: s_inv_complex	(ao_num,ao_num)


    Inverse of the overlap matrix

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap_complex`
       * :c:data:`lin_dep_cutoff`


 
.. c:var:: use_pw


    File : :file:`ao_one_e_ints/aos_cgtos.irp.f`

    .. code:: fortran

        logical	:: use_pw	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_expo_cgtos_ord_transp`
       * :c:data:`ao_expo_cgtos_ord_transp`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_2e_cgtos_schwartz`
       * :c:data:`ao_integrals_n_e_cgtos`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: ao_one_e_integral_zero:


    File : :file:`ao_one_e_ints/screening.irp.f`

    .. code:: fortran

        logical function ao_one_e_integral_zero(i,k)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_one_e_integrals_threshold`
       * :c:data:`ao_overlap_abs`
       * :c:data:`io_ao_integrals_overlap`
       * :c:data:`is_periodic`
       * :c:data:`use_cgtos`

 
.. c:function:: give_all_erf_kl_ao:


    File : :file:`ao_one_e_ints/pot_ao_erf_ints.irp.f`

    .. code:: fortran

        subroutine give_all_erf_kl_ao(integrals_ao,mu_in,C_center)


    Subroutine that returns all integrals over $r$ of type
    $\frac{ \erf(\mu * | r - R_C | ) }{ | r - R_C | }$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`

 
.. c:function:: give_polynomial_mult_center_one_e:


    File : :file:`ao_one_e_ints/pot_ao_ints.irp.f`

    .. code:: fortran

        subroutine give_polynomial_mult_center_one_e(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out)


    Returns the explicit polynomial in terms of the "t" variable of the following
    
    $I_{x1}(a_x, d_x,p,q) \times I_{x1}(a_y, d_y,p,q) \times I_{x1}(a_z, d_z,p,q)$.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`nai_pol_mult`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e`
       * :c:func:`multiply_poly`

 
.. c:function:: give_polynomial_mult_center_one_e_erf:


    File : :file:`ao_one_e_ints/pot_ao_erf_ints.irp.f`

    .. code:: fortran

        subroutine give_polynomial_mult_center_one_e_erf(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out,mu_in)


    Returns the explicit polynomial in terms of the $t$ variable of the
    following polynomial:
    
    $I_{x1}(a_x, d_x,p,q) \times I_{x1}(a_y, d_y,p,q) \times I_{x1}(a_z, d_z,p,q)$.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e`
       * :c:func:`multiply_poly`

 
.. c:function:: give_polynomial_mult_center_one_e_erf_opt:


    File : :file:`ao_one_e_ints/pot_ao_erf_ints.irp.f`

    .. code:: fortran

        subroutine give_polynomial_mult_center_one_e_erf_opt(A_center, B_center, power_A, power_B, C_center, n_pt_in, d, n_pt_out, p_inv_2, p_new, P_center)


    Returns the explicit polynomial in terms of the $t$ variable of the
    following polynomial:
    
    $I_{x1}(a_x, d_x,p,q) \times I_{x1}(a_y, d_y,p,q) \times I_{x1}(a_z, d_z,p,q)$.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`nai_pol_mult_erf`
       * :c:func:`nai_pol_mult_erf_v`
       * :c:func:`nai_pol_mult_erf_with1s`
       * :c:func:`nai_pol_mult_erf_with1s_v`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_x1_pol_mult_one_e`
       * :c:func:`multiply_poly`

 
.. c:function:: int_gaus_pol:


    File : :file:`ao_one_e_ints/pot_ao_ints.irp.f`

    .. code:: fortran

        double precision function int_gaus_pol(alpha,n)


    Computes the integral:
    
    $\int_{-\infty}^{\infty} x^n \exp(-\alpha x^2) dx$.

 
.. c:function:: nai_pol_mult:


    File : :file:`ao_one_e_ints/pot_ao_ints.irp.f`

    .. code:: fortran

        double precision function NAI_pol_mult(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in)


    Computes the electron-nucleus attraction with two primitves.
    
    :math:`\langle g_i | \frac{1}{|r-R_c|} | g_j \rangle`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_polynomial_mult_center_one_e`

 
.. c:function:: nai_pol_mult_cosgtos:


    File : :file:`ao_one_e_ints/one_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        complex*16 function NAI_pol_mult_cosgtos(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in)


    
    Computes the electron-nucleus attraction with two primitves cosgtos.
    
    :math:`\langle g_i | \frac{1}{|r-R_c|} | g_j \rangle`
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_cpolynomial_mult_center_one_e`

 
.. c:function:: nai_pol_mult_erf:


    File : :file:`ao_one_e_ints/pot_ao_erf_ints.irp.f`

    .. code:: fortran

        double precision function NAI_pol_mult_erf(A_center, B_center, power_A, power_B, alpha, beta, C_center, n_pt_in, mu_in)


    
    Computes the following integral :
    
    .. math::
    
      \int dr (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
      \frac{\erf(\mu |r - R_C |)}{| r - R_C |}$.
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_polynomial_mult_center_one_e_erf_opt`

 
.. c:function:: nai_pol_mult_erf_ao:


    File : :file:`ao_one_e_ints/pot_ao_erf_ints.irp.f`

    .. code:: fortran

        double precision function NAI_pol_mult_erf_ao(i_ao, j_ao, mu_in, C_center)


    
    Computes the following integral :
    $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) \frac{\erf(\mu |r - R_C|)}{|r - R_C|}$.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`

 
.. c:function:: nai_pol_mult_erf_ao_with1s:


    File : :file:`ao_one_e_ints/pot_ao_erf_ints.irp.f`

    .. code:: fortran

        double precision function NAI_pol_mult_erf_ao_with1s(i_ao, j_ao, beta, B_center, mu_in, C_center)


    
    Computes the following integral :
    $\int_{-\infty}^{infty} dr \chi_i(r) \chi_j(r) e^{-\beta (r - B_center)^2} \frac{\erf(\mu |r - R_C|)}{|r - R_C|}$.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_coef_normalized_ordered_transp`
       * :c:data:`ao_expo_ordered_transp`
       * :c:data:`ao_nucl`
       * :c:data:`ao_power`
       * :c:data:`ao_prim_num`
       * :c:data:`n_pt_max_integrals`
       * :c:data:`nucl_coord`

 
.. c:function:: nai_pol_mult_erf_v:


    File : :file:`ao_one_e_ints/pot_ao_erf_ints.irp.f`

    .. code:: fortran

        subroutine NAI_pol_mult_erf_v(A_center, B_center, power_A, power_B, alpha, beta, C_center, LD_C, n_pt_in, mu_in, res_v, LD_resv, n_points)


    
    Computes the following integral :
    
    .. math::
    
      \int dr (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
      \frac{\erf(\mu |r - R_C |)}{| r - R_C |}$.
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_polynomial_mult_center_one_e_erf_opt`

 
.. c:function:: nai_pol_mult_erf_with1s_v:


    File : :file:`ao_one_e_ints/pot_ao_erf_ints.irp.f`

    .. code:: fortran

        subroutine NAI_pol_mult_erf_with1s_v(A1_center, A2_center, power_A1, power_A2, alpha1, alpha2, beta, B_center, LD_B, C_center, LD_C, n_pt_in, mu_in, res_v, LD_resv, n_points)


    
    Computes the following integral :
    
    .. math                      ::
    
      \int dx (x - A1_x)^a_1 (x - B1_x)^a_2 \exp(-\alpha_1 (x - A1_x)^2 - \alpha_2 (x - A2_x)^2)
      \int dy (y - A1_y)^b_1 (y - B1_y)^b_2 \exp(-\alpha_1 (y - A1_y)^2 - \alpha_2 (y - A2_y)^2)
      \int dz (x - A1_z)^c_1 (z - B1_z)^c_2 \exp(-\alpha_1 (z - A1_z)^2 - \alpha_2 (z - A2_z)^2)
      \exp(-\beta (r - B)^2)
      \frac{\erf(\mu |r - R_C|)}{|r - R_C|}$.
    

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`give_polynomial_mult_center_one_e_erf_opt`

 
.. c:function:: overlap_bourrin_deriv_x:


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        subroutine overlap_bourrin_deriv_x(i_component,A_center,B_center,alpha,beta,power_A,power_B,dx,lower_exp_val,overlap_x,nx)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv_1_x`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_bourrin_x`

 
.. c:function:: overlap_bourrin_dipole:


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        subroutine overlap_bourrin_dipole(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,lower_exp_val,dx,nx)



    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_dipole_x`

 
.. c:function:: overlap_bourrin_spread:


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        subroutine overlap_bourrin_spread(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,lower_exp_val,dx,nx)


    Computes the following integral :
     int [-infty ; +infty] of [(x-A_center)^(power_A) * (x-B_center)^power_B * exp(-alpha(x-A_center)^2) * exp(-beta(x-B_center)^2) * x^2 ]
     needed for the dipole and those things

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_spread_x`

 
.. c:function:: overlap_bourrin_x:


    File : :file:`ao_one_e_ints/spread_dipole_ao.irp.f`

    .. code:: fortran

        subroutine overlap_bourrin_x(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,lower_exp_val,dx,nx)



    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`overlap_bourrin_deriv_x`

 
.. c:function:: v_n_e:


    File : :file:`ao_one_e_ints/pot_ao_ints.irp.f`

    .. code:: fortran

        double precision function V_n_e(a_x,a_y,a_z,b_x,b_y,b_z,alpha,beta)


    Primitve nuclear attraction between the two primitves centered on the same atom.
    
    $p_1 = x^{a_x} y^{a_y} z^{a_z} \exp(-\alpha r^2)$
    
    $p_2 = x^{b_x} y^{b_y} z^{b_z} \exp(-\beta  r^2)$

 
.. c:function:: v_n_e_cgtos:


    File : :file:`ao_one_e_ints/one_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        complex*16 function V_n_e_cgtos(a_x, a_y, a_z, b_x, b_y, b_z, alpha, beta)


    Primitve nuclear attraction between the two primitves centered on the same atom.
    
    $p_1 = x^{a_x} y^{a_y} z^{a_z} \exp(-\alpha r^2)$
    
    $p_2 = x^{b_x} y^{b_y} z^{b_z} \exp(-\beta  r^2)$

 
.. c:function:: v_n_e_cosgtos:


    File : :file:`ao_one_e_ints/one_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        complex*16 function V_n_e_cosgtos(a_x, a_y, a_z, b_x, b_y, b_z, alpha, beta)


    Primitve nuclear attraction between the two primitves centered on the same atom.
    
    $p_1 = x^{a_x} y^{a_y} z^{a_z} \exp(-\alpha r^2)$
    
    $p_2 = x^{b_x} y^{b_y} z^{b_z} \exp(-\beta  r^2)$

 
.. c:function:: v_r:


    File : :file:`ao_one_e_ints/pot_ao_ints.irp.f`

    .. code:: fortran

        double precision function V_r(n,alpha)


    Computes the radial part of the nuclear attraction integral:
    
    $\int_{0}^{\infty} r^n  \exp(-\alpha  r^2)  dr$
    

 
.. c:function:: v_r_cgtos:


    File : :file:`ao_one_e_ints/one_e_coul_integrals_cgtos.irp.f`

    .. code:: fortran

        complex*16 function V_r_cgtos(n, alpha)


    Computes the radial part of the nuclear attraction integral:
    
    $\int_{0}^{\infty} r^n  \exp(-\alpha  r^2)  dr$
    

 
.. c:function:: v_r_cosgtos:


    File : :file:`ao_one_e_ints/one_e_Coul_integrals_cosgtos.irp.f`

    .. code:: fortran

        complex*16 function V_r_cosgtos(n, alpha)


    Computes the radial part of the nuclear attraction integral:
    
    $\int_{0}^{\infty} r^n  \exp(-\alpha  r^2)  dr$
    

