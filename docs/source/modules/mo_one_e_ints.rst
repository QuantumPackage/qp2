.. _module_mo_one_e_ints: 
 
.. program:: mo_one_e_ints 
 
.. default-role:: option 
 
==================
mo_one_e_integrals
==================

All the one-electron integrals in |MO| basis are defined here.

The most important providers for usual quantum-chemistry calculation are:

* `mo_kinetic_integrals` which are the kinetic operator integrals on the |AO| basis (see :file:`kin_mo_ints.irp.f`)
* `mo_integrals_n_e` which are the nuclear-elctron operator integrals on the |AO| basis (see :file:`pot_mo_ints.irp.f`)
* `mo_one_e_integrals` which are the the h_core operator integrals on the |AO| basis (see :file:`mo_mono_ints.irp.f`)

Note that you can find other interesting integrals related to the position operator in :file:`spread_dipole_mo.irp.f`.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: mo_integrals_e_n
 
    Nucleus-electron integrals in |MO| basis set
 
 
.. option:: io_mo_integrals_e_n
 
    Read/Write |MO| electron-nucleus attraction integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: mo_integrals_kinetic
 
    Kinetic energy integrals in |MO| basis set
 
 
.. option:: io_mo_integrals_kinetic
 
    Read/Write |MO| one-electron kinetic integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: mo_integrals_pseudo
 
    Pseudopotential integrals in |MO| basis set
 
 
.. option:: io_mo_integrals_pseudo
 
    Read/Write |MO| pseudopotential integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
.. option:: mo_one_e_integrals
 
    One-electron integrals in |MO| basis set
 
 
.. option:: io_mo_one_e_integrals
 
    Read/Write |MO| one-electron integrals from/to disk [ Write | Read | None ]
 
    Default: None
 
 
Providers 
--------- 
 
.. c:var:: mo_dipole_x


    File : :file:`mo_one_e_ints/spread_dipole_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_dipole_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_dipole_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_dipole_z	(mo_num,mo_num)


    array of the integrals of MO_i * x MO_j
    array of the integrals of MO_i * y MO_j
    array of the integrals of MO_i * z MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_dipole_x`
       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: mo_dipole_y


    File : :file:`mo_one_e_ints/spread_dipole_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_dipole_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_dipole_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_dipole_z	(mo_num,mo_num)


    array of the integrals of MO_i * x MO_j
    array of the integrals of MO_i * y MO_j
    array of the integrals of MO_i * z MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_dipole_x`
       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: mo_dipole_z


    File : :file:`mo_one_e_ints/spread_dipole_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_dipole_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_dipole_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_dipole_z	(mo_num,mo_num)


    array of the integrals of MO_i * x MO_j
    array of the integrals of MO_i * y MO_j
    array of the integrals of MO_i * z MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_dipole_x`
       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: mo_integrals_n_e


    File : :file:`mo_one_e_ints/pot_mo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_integrals_n_e	(mo_num,mo_num)


    Nucleus-electron interaction on the |MO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e`
       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`read_mo_integrals_n_e`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_one_e_integrals`
       * :c:data:`ref_bitmask_energy`

 
.. c:var:: mo_integrals_n_e_per_atom


    File : :file:`mo_one_e_ints/pot_mo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_integrals_n_e_per_atom	(mo_num,mo_num,nucl_num)


    mo_integrals_n_e_per_atom(i,j,k) =
    :math:`\langle \phi_i| -\frac{1}{|r-R_k|} | \phi_j \rangle` .
    where R_k is the coordinate of the k-th nucleus.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_integrals_n_e_per_atom`
       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`nucl_num`


 
.. c:var:: mo_kinetic_integrals


    File : :file:`mo_one_e_ints/kin_mo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_kinetic_integrals	(mo_num,mo_num)


    Kinetic energy integrals in the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_kinetic_integrals`
       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`read_mo_integrals_kinetic`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_one_e_integrals`
       * :c:data:`ref_bitmask_energy`

 
.. c:var:: mo_one_e_integrals


    File : :file:`mo_one_e_ints/mo_one_e_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_one_e_integrals	(mo_num,mo_num)


    array of the one-electron Hamiltonian on the |MO| basis :
    sum of the kinetic and nuclear electronic potentials (and pseudo potential if needed)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`read_mo_one_e_integrals`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`psi_energy_h_core`
       * :c:data:`ref_bitmask_energy`

 
.. c:var:: mo_overlap


    File : :file:`mo_one_e_ints/mo_overlap.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_overlap	(mo_num,mo_num)


    Provider to check that the MOs are indeed orthonormal.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: mo_pseudo_integrals


    File : :file:`mo_one_e_ints/pot_mo_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_pseudo_integrals	(mo_num,mo_num)


    Pseudopotential integrals in |MO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_pseudo_integrals`
       * :c:data:`do_pseudo`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`read_mo_integrals_pseudo`


 
.. c:var:: mo_spread_x


    File : :file:`mo_one_e_ints/spread_dipole_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_spread_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_z	(mo_num,mo_num)


    array of the integrals of MO_i * x^2 MO_j
    array of the integrals of MO_i * y^2 MO_j
    array of the integrals of MO_i * z^2 MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_spread_x`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: mo_spread_y


    File : :file:`mo_one_e_ints/spread_dipole_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_spread_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_z	(mo_num,mo_num)


    array of the integrals of MO_i * x^2 MO_j
    array of the integrals of MO_i * y^2 MO_j
    array of the integrals of MO_i * z^2 MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_spread_x`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: mo_spread_z


    File : :file:`mo_one_e_ints/spread_dipole_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_spread_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_z	(mo_num,mo_num)


    array of the integrals of MO_i * x^2 MO_j
    array of the integrals of MO_i * y^2 MO_j
    array of the integrals of MO_i * z^2 MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_spread_x`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: s_mo_coef


    File : :file:`mo_one_e_ints/ao_to_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: s_mo_coef	(ao_num,mo_num)


    Product S.C where S is the overlap matrix in the AO basis and C the mo_coef matrix.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: mo_to_ao:


    File : :file:`mo_one_e_ints/ao_to_mo.irp.f`

    .. code:: fortran

        subroutine mo_to_ao(A_mo,LDA_mo,A_ao,LDA_ao)


    Transform A from the MO basis to the AO basis
    
    $(S.C).A_{mo}.(S.C)^\dagger$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_num`
       * :c:data:`s_mo_coef`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: mo_to_ao_no_overlap:


    File : :file:`mo_one_e_ints/ao_to_mo.irp.f`

    .. code:: fortran

        subroutine mo_to_ao_no_overlap(A_mo,LDA_mo,A_ao,LDA_ao)


    $C.A_{mo}.C^\dagger$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`one_e_dm_alpha_ao_for_dft_no_core`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: orthonormalize_mos:


    File : :file:`mo_one_e_ints/orthonormalize.irp.f`

    .. code:: fortran

        subroutine orthonormalize_mos



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`lin_dep_cutoff`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_overlap`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`save_natural_mos`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`nullify_small_elements`
       * :c:func:`ortho_lowdin`
       * :c:func:`restore_symmetry`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_coef`

