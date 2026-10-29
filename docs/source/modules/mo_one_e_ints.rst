.. _module_mo_one_e_ints: 
 
.. program:: mo_one_e_ints 
 
.. default-role:: option 
 
==================
mo_one_e_ints
==================

All the one-electron integrals in the |MO| basis.

They are obtained from the |AO| integrals of the
:ref:`module_ao_one_e_ints` module with the |MO| transformation
(:c:func:`mo_to_ao` and its variants, and
:c:data:`ao_one_e_integrals_from_mo`):

* :c:data:`mo_kinetic_integrals`: the kinetic operator integrals
  (`kin_mo_ints.irp.f`)
* :c:data:`mo_integrals_n_e`: the nuclear-electron integrals
  (`pot_mo_ints.irp.f`)
* :c:data:`mo_one_e_integrals`: the core Hamiltonian integrals
  (`mo_one_e_ints.irp.f`)
* :c:data:`mo_pseudo_integrals`, with the local and non local parts
  (`pot_mo_pseudo_ints.irp.f`)

Other integrals related to the position operator, dipole moments and
spreads, are in `spread_dipole_mo.irp.f`, the orthonormalization of
the |MOs| in `orthonormalize.irp.f`, and the reading and writing of
the integrals in the |EZFIO| database are provided by the
`read_mo_integrals_*` and `write_mo_integrals_*` subroutines.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: mo_integrals_n_e
 
    Nucleus-electron integrals in |MO| basis set
 
 
.. option:: io_mo_integrals_n_e
 
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
 
.. c:var:: double ao_one_e_integrals_from_mo


    File : :file:`mo_one_e_ints/ao_to_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_one_e_integrals_from_mo	(ao_num,ao_num)


    Integrals of the one e hamiltonian obtained from the integrals on the MO basis
    
    WARNING : this is equal to ao_one_e_integrals only if the AO and MO basis have the same number of functions

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`s_mo_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_sphe_one_e_integrals_from_mo`

 
.. c:var:: double ao_sphe_one_e_integrals_from_mo


    File : :file:`mo_one_e_ints/ao_to_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ao_sphe_one_e_integrals_from_mo	(ao_sphe_num,ao_sphe_num)


    Integrals of the one e hamiltonian obtained from the integrals on the MO basis
    
    WARNING : this is equal to ao_one_e_integrals only if the AO and MO basis have the same number of functions

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_one_e_integrals_from_mo`
       * :c:data:`ao_sphe_num`
       * :c:data:`mo_num`
       * :c:data:`mo_one_e_integrals`
       * :c:data:`s_mo_sphe_coef`


 
.. c:var:: double mo_dipole_x


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_prop_dipole_x`
       * :c:data:`mo_spread_centered_x`
       * :c:data:`z_dipole_moment`

 
.. c:var:: double mo_dipole_y


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_prop_dipole_x`
       * :c:data:`mo_spread_centered_x`
       * :c:data:`z_dipole_moment`

 
.. c:var:: double mo_dipole_z


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_prop_dipole_x`
       * :c:data:`mo_spread_centered_x`
       * :c:data:`z_dipole_moment`

 
.. c:var:: double mo_integrals_n_e


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
       * :c:data:`ref_bitmask_n_e_energy`
       * :c:data:`v_ne_psi_energy`

 
.. c:var:: double mo_integrals_n_e_per_atom


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


 
.. c:var:: double mo_kinetic_integrals


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
       * :c:data:`ref_bitmask_kinetic_energy`

 
.. c:var:: double mo_one_e_integrals


    File : :file:`mo_one_e_ints/mo_one_e_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_one_e_integrals	(mo_num,mo_num)


    array of the one-electron Hamiltonian on the |MO| basis :
    sum of the kinetic and nuclear electronic potentials (and pseudo potential if needed)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`mo_coef`
       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`read_ao_one_e_integrals`
       * :c:data:`read_mo_one_e_integrals`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`ao_one_e_integrals_from_mo`
       * :c:data:`ao_sphe_one_e_integrals_from_mo`
       * :c:data:`core_energy`
       * :c:data:`etwo`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`h_core_ri`
       * :c:data:`one_ints_no`
       * :c:data:`psi_energy_h_core`
       * :c:data:`ref_bitmask_one_e_energy`

 
.. c:var:: double mo_overlap


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


 
.. c:var:: double mo_pseudo_integrals


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


 
.. c:var:: double mo_pseudo_integrals_local


    File : :file:`mo_one_e_ints/pot_mo_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_pseudo_integrals_local	(mo_num,mo_num)


    Pseudopotential integrals in |MO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_pseudo_integrals_local`
       * :c:data:`do_pseudo`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: double mo_pseudo_integrals_non_local


    File : :file:`mo_one_e_ints/pot_mo_pseudo_ints.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_pseudo_integrals_non_local	(mo_num,mo_num)


    Pseudopotential integrals in |MO| basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_pseudo_integrals_non_local`
       * :c:data:`do_pseudo`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: double mo_spread_centered_x


    File : :file:`mo_one_e_ints/spread_dipole_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_spread_centered_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_centered_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_centered_z	(mo_num,mo_num)


    array of the integrals of MO_i * (x^2 - <MO_i|x|MO_j>^2) MO_j = MO_i x^2 MO_j - (MO_i x MO_j)^2
    array of the integrals of MO_i * (y^2 - <MO_i|y|MO_j>^2) MO_j = MO_i y^2 MO_j - (MO_i y MO_j)^2
    array of the integrals of MO_i * (z^2 - <MO_i|z|MO_j>^2) MO_j = MO_i z^2 MO_j - (MO_i z MO_j)^2

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_dipole_x`
       * :c:data:`mo_num`
       * :c:data:`mo_spread_x`


 
.. c:var:: double mo_spread_centered_y


    File : :file:`mo_one_e_ints/spread_dipole_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_spread_centered_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_centered_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_centered_z	(mo_num,mo_num)


    array of the integrals of MO_i * (x^2 - <MO_i|x|MO_j>^2) MO_j = MO_i x^2 MO_j - (MO_i x MO_j)^2
    array of the integrals of MO_i * (y^2 - <MO_i|y|MO_j>^2) MO_j = MO_i y^2 MO_j - (MO_i y MO_j)^2
    array of the integrals of MO_i * (z^2 - <MO_i|z|MO_j>^2) MO_j = MO_i z^2 MO_j - (MO_i z MO_j)^2

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_dipole_x`
       * :c:data:`mo_num`
       * :c:data:`mo_spread_x`


 
.. c:var:: double mo_spread_centered_z


    File : :file:`mo_one_e_ints/spread_dipole_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_spread_centered_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_centered_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_spread_centered_z	(mo_num,mo_num)


    array of the integrals of MO_i * (x^2 - <MO_i|x|MO_j>^2) MO_j = MO_i x^2 MO_j - (MO_i x MO_j)^2
    array of the integrals of MO_i * (y^2 - <MO_i|y|MO_j>^2) MO_j = MO_i y^2 MO_j - (MO_i y MO_j)^2
    array of the integrals of MO_i * (z^2 - <MO_i|z|MO_j>^2) MO_j = MO_i z^2 MO_j - (MO_i z MO_j)^2

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_dipole_x`
       * :c:data:`mo_num`
       * :c:data:`mo_spread_x`


 
.. c:var:: double mo_spread_x


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_spread_centered_x`

 
.. c:var:: double mo_spread_y


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_spread_centered_x`

 
.. c:var:: double mo_spread_z


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_spread_centered_x`

 
.. c:var:: double s_mo_coef


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

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_one_e_integrals_from_mo`
       * :c:data:`fock_matrix_ao`
       * :c:data:`s_mo_sphe_coef`

 
.. c:var:: double s_mo_sphe_coef


    File : :file:`mo_one_e_ints/ao_to_mo.irp.f`

    .. code:: fortran

        double precision, allocatable	:: s_mo_sphe_coef	(ao_sphe_num,mo_num)


    Product S.C where S is the overlap matrix in the AO basis and C the mo_sphe_coef matrix.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_overlap`
       * :c:data:`ao_sphe_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`mo_sphe_coef`
       * :c:data:`s_mo_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_sphe_one_e_integrals_from_mo`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void mo_to_ao(void)


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

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_one_e_integrals_from_mo`
       * :c:data:`fock_matrix_ao`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void mo_to_ao_no_overlap(void)


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

       * :c:data:`d0tu_alpha_ao`
       * :c:data:`one_e_dm_alpha_ao_for_dft`
       * :c:data:`one_e_dm_alpha_ao_for_dft_no_core`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void mo_to_ao_sphe(void)


    File : :file:`mo_one_e_ints/ao_to_mo.irp.f`

    .. code:: fortran

        subroutine mo_to_ao_sphe(A_mo,LDA_mo,A_ao,LDA_ao)


    Transform A from the MO basis to the AO basis
    
    $(S.C).A_{mo}.(S.C)^\dagger$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_sphe_num`
       * :c:data:`mo_num`
       * :c:data:`s_mo_sphe_coef`

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_sphe_one_e_integrals_from_mo`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void mo_to_ao_sphe_no_overlap(void)


    File : :file:`mo_one_e_ints/ao_to_mo.irp.f`

    .. code:: fortran

        subroutine mo_to_ao_sphe_no_overlap(A_mo,LDA_mo,A_ao,LDA_ao)


    $C.A_{mo}.C^\dagger$

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_sphe_num`
       * :c:data:`mo_num`
       * :c:data:`mo_sphe_coef`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`

 
.. c:function:: void orthonormalize_mos(void)


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

       * :c:func:`huckel_guess`
       * :c:func:`roothaan_hall_scf`
       * :c:func:`save_natural_mos`
       * :c:func:`save_natural_mos_canon_label`
       * :c:func:`save_natural_mos_no_ov_rot`
       * :c:func:`scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ortho_lowdin`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`mo_coef`

