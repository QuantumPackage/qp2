.. _module_bitmask: 
 
.. program:: bitmask 
 
.. default-role:: option 
 
==============
bitmask module
==============

The central part of this module is the :file:`bitmasks_module.f90` file. It contains
the constants that will be used to define on which kind of integer the bitmasks
will be defined.

In the program, to represent a determinant as a pair of bitstrings,
the determinant should be defined as

.. code-block:: fortran

  use bitmasks
  integer(bit_kind)  :: determinant(N_int,2)


:file:`bitmasks_routines.irp.f` contains helper routines to manipulate bitmask, like
transforming a bit string to a list of integers for example.


`bit_kind_shift`, `bit_kind_size` and `bit_kind` are supposed to be consistent::

   2**bit_kind_shift = bit_kind_size
   bit_kind = bit_kind_size / 8


For an example of how to use the bitmaks, see the file :file:`example.irp.f`.
 
 
 
Providers 
--------- 
 
.. c:var:: act_bitmask


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: cas_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: cas_bitmask	(N_int,2,N_cas_bitmask)


    Bitmasks for CAS reference determinants. (N_int, alpha/beta, CAS reference)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`full_ijkl_bitmask`
       * :c:data:`generators_bitmask_restart`
       * :c:data:`hf_bitmask`
       * :c:data:`mpi_master`
       * :c:data:`n_cas_bitmask`
       * :c:data:`n_generators_bitmask`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`closed_shell_ref_bitmask`
       * :c:data:`psi_cas`
       * :c:data:`reunion_of_bitmask`

 
.. c:var:: closed_shell_ref_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: closed_shell_ref_bitmask	(N_int,2)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`n_int`
       * :c:data:`ref_bitmask`


 
.. c:var:: core_bitmask


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: core_inact_act_bitmask_4


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: core_inact_act_bitmask_4	(N_int,4)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_act_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_in_map`

 
.. c:var:: core_inact_virt_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: inact_virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: core_inact_virt_bitmask	(N_int,2)


    Reunion of the inactive and virtual bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_int`


 
.. c:var:: del_bitmask


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: dim_list_act_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_core_orb	
        integer	:: dim_list_inact_orb	
        integer	:: dim_list_virt_orb	
        integer	:: dim_list_act_orb	
        integer	:: dim_list_del_orb	


    dimensions for the allocation of list_inact, list_virt, list_core and list_act
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`

 
.. c:var:: dim_list_core_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_core_orb	
        integer	:: dim_list_inact_orb	
        integer	:: dim_list_virt_orb	
        integer	:: dim_list_act_orb	
        integer	:: dim_list_del_orb	


    dimensions for the allocation of list_inact, list_virt, list_core and list_act
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`

 
.. c:var:: dim_list_del_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_core_orb	
        integer	:: dim_list_inact_orb	
        integer	:: dim_list_virt_orb	
        integer	:: dim_list_act_orb	
        integer	:: dim_list_del_orb	


    dimensions for the allocation of list_inact, list_virt, list_core and list_act
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`

 
.. c:var:: dim_list_inact_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_core_orb	
        integer	:: dim_list_inact_orb	
        integer	:: dim_list_virt_orb	
        integer	:: dim_list_act_orb	
        integer	:: dim_list_del_orb	


    dimensions for the allocation of list_inact, list_virt, list_core and list_act
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`

 
.. c:var:: dim_list_virt_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: dim_list_core_orb	
        integer	:: dim_list_inact_orb	
        integer	:: dim_list_virt_orb	
        integer	:: dim_list_act_orb	
        integer	:: dim_list_del_orb	


    dimensions for the allocation of list_inact, list_virt, list_core and list_act
    it is at least 1

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`

 
.. c:var:: full_ijkl_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: full_ijkl_bitmask	(N_int)


    Bitmask to include all possible MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`full_ijkl_bitmask_4`
       * :c:data:`generators_bitmask`
       * :c:data:`generators_bitmask_restart`

 
.. c:var:: full_ijkl_bitmask_4


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: full_ijkl_bitmask_4	(N_int,4)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`full_ijkl_bitmask`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`

 
.. c:var:: generators_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: generators_bitmask	(N_int,2,6,N_generators_bitmask)


    Bitmasks for generator determinants.
    (N_int, alpha/beta, hole/particle, generator).
    
    3rd index is :
    
    * 1 : hole     for single exc
    
    * 2 : particle for single exc
    
    * 3 : hole     for 1st exc of double
    
    * 4 : particle for 1st exc of double
    
    * 5 : hole     for 2nd exc of double
    
    * 6 : particle for 2nd exc of double
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`full_ijkl_bitmask`
       * :c:data:`mpi_master`
       * :c:data:`n_generators_bitmask`
       * :c:data:`n_int`


 
.. c:var:: generators_bitmask_restart


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: generators_bitmask_restart	(N_int,2,6,N_generators_bitmask_restart)


    Bitmasks for generator determinants.
    (N_int, alpha/beta, hole/particle, generator).
    
    3rd index is :
    
    * 1 : hole     for single exc
    
    * 2 : particle for single exc
    
    * 3 : hole     for 1st exc of double
    
    * 4 : particle for 1st exc of double
    
    * 5 : hole     for 2nd exc of double
    
    * 6 : particle for 2nd exc of double
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`full_ijkl_bitmask`
       * :c:data:`mpi_master`
       * :c:data:`n_generators_bitmask`
       * :c:data:`n_generators_bitmask_restart`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`

 
.. c:var:: hf_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: hf_bitmask	(N_int,2)


    Hartree Fock bit mask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`degree_max_generators`
       * :c:data:`double_exc_bitmask`
       * :c:data:`max_degree_exc`
       * :c:data:`psi_cas`
       * :c:data:`psi_det`
       * :c:data:`ref_bitmask`
       * :c:data:`single_exc_bitmask`
       * :c:data:`unpaired_alpha_electrons`

 
.. c:var:: i_bitmask_gen


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer	:: i_bitmask_gen	


    Current bitmask for the generators


 
.. c:var:: inact_bitmask


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: inact_virt_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: inact_virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: core_inact_virt_bitmask	(N_int,2)


    Reunion of the inactive and virtual bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_int`


 
.. c:var:: index_holes_bitmask


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        integer, allocatable	:: index_holes_bitmask	(3)


    Index of the holes in the generators_bitmasks


 
.. c:var:: index_particl_bitmask


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        integer, allocatable	:: index_particl_bitmask	(3)


    Index of the holes in the generators_bitmasks


 
.. c:var:: list_act


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: list_act_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: list_core


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: list_core_inact_act


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_core_inact_act	(n_core_inact_act_orb)
        integer, allocatable	:: list_core_inact_act_reverse	(n_core_inact_act_orb)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_act_bitmask`


 
.. c:var:: list_core_inact_act_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_core_inact_act	(n_core_inact_act_orb)
        integer, allocatable	:: list_core_inact_act_reverse	(n_core_inact_act_orb)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_act_bitmask`


 
.. c:var:: list_core_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: list_del


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: list_del_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: list_inact


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: list_inact_act


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact_act	(n_inact_act_orb)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_act_orb`


 
.. c:var:: list_inact_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: list_virt


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: list_virt_reverse


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: mpi_bit_kind


    File : :file:`bitmask/mpi.irp.f`

    .. code:: fortran

        integer	:: mpi_bit_kind	


    MPI bit kind type


 
.. c:var:: n_act_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_core_orb	
        integer	:: n_inact_orb	
        integer	:: n_act_orb	
        integer	:: n_virt_orb	
        integer	:: n_del_orb	


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`dim_list_core_orb`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb_allocate`
       * :c:data:`n_inact_act_orb`
       * :c:data:`n_inact_orb_allocate`
       * :c:data:`n_virt_orb_allocate`
       * :c:data:`pt2_f`

 
.. c:var:: n_cas_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_cas_bitmask	


    Number of bitmasks for CAS

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`psi_cas`

 
.. c:var:: n_core_inact_act_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_core_inact_act_orb	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_core_inact_act`

 
.. c:var:: n_core_inact_orb


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_core_inact_orb	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_bitmask`


 
.. c:var:: n_core_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_core_orb	
        integer	:: n_inact_orb	
        integer	:: n_act_orb	
        integer	:: n_virt_orb	
        integer	:: n_del_orb	


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`dim_list_core_orb`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb_allocate`
       * :c:data:`n_inact_act_orb`
       * :c:data:`n_inact_orb_allocate`
       * :c:data:`n_virt_orb_allocate`
       * :c:data:`pt2_f`

 
.. c:var:: n_core_orb_allocate


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_core_orb_allocate	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`


 
.. c:var:: n_del_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_core_orb	
        integer	:: n_inact_orb	
        integer	:: n_act_orb	
        integer	:: n_virt_orb	
        integer	:: n_del_orb	


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`dim_list_core_orb`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb_allocate`
       * :c:data:`n_inact_act_orb`
       * :c:data:`n_inact_orb_allocate`
       * :c:data:`n_virt_orb_allocate`
       * :c:data:`pt2_f`

 
.. c:var:: n_generators_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_generators_bitmask	


    Number of bitmasks for generators

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`generators_bitmask`
       * :c:data:`generators_bitmask_restart`

 
.. c:var:: n_generators_bitmask_restart


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_generators_bitmask_restart	


    Number of bitmasks for generators

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask_restart`

 
.. c:var:: n_inact_act_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_inact_act_orb	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact_act`

 
.. c:var:: n_inact_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_core_orb	
        integer	:: n_inact_orb	
        integer	:: n_act_orb	
        integer	:: n_virt_orb	
        integer	:: n_del_orb	


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`dim_list_core_orb`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb_allocate`
       * :c:data:`n_inact_act_orb`
       * :c:data:`n_inact_orb_allocate`
       * :c:data:`n_virt_orb_allocate`
       * :c:data:`pt2_f`

 
.. c:var:: n_inact_orb_allocate


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_inact_orb_allocate	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`


 
.. c:var:: n_int


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_int	


    Number of 64-bit integers needed to represent determinants as binary strings

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`ci_electronic_energy`
       * :c:data:`closed_shell_ref_bitmask`
       * :c:data:`coef_hf_selector`
       * :c:data:`core_inact_act_bitmask_4`
       * :c:data:`degree_max_generators`
       * :c:data:`det_to_occ_pattern`
       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`double_exc_bitmask`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`fock_operator_closed_shell_ref_bitmask`
       * :c:data:`fock_wee_closed_shell`
       * :c:data:`full_ijkl_bitmask`
       * :c:data:`full_ijkl_bitmask_4`
       * :c:data:`generators_bitmask`
       * :c:data:`generators_bitmask_restart`
       * :c:data:`global_selection_buffer`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`h_matrix_all_dets`
       * :c:data:`h_matrix_cas`
       * :c:data:`hf_bitmask`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact`
       * :c:data:`max_degree_exc`
       * :c:data:`mo_two_e_integrals_erf_in_map`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`n_cas_bitmask`
       * :c:data:`n_core_inact_orb`
       * :c:data:`n_generators_bitmask`
       * :c:data:`n_generators_bitmask_restart`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_cas`
       * :c:data:`psi_cas_sorted_bit`
       * :c:data:`psi_det`
       * :c:data:`psi_det_alpha`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_det_hii`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det_sorted_bit`
       * :c:data:`psi_det_sorted_gen`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy_two_e`
       * :c:data:`psi_non_cas`
       * :c:data:`psi_non_cas_sorted_bit`
       * :c:data:`psi_occ_pattern`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_diag_h_mat`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`
       * :c:data:`ref_closed_shell_bitmask`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`s2_matrix_all_dets`
       * :c:data:`s2_values`
       * :c:data:`single_exc_bitmask`
       * :c:data:`singles_alpha_csc`
       * :c:data:`singles_alpha_csc_idx`
       * :c:data:`singles_beta_csc`
       * :c:data:`singles_beta_csc_idx`
       * :c:data:`unpaired_alpha_electrons`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: n_virt_orb


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer	:: n_core_orb	
        integer	:: n_inact_orb	
        integer	:: n_act_orb	
        integer	:: n_virt_orb	
        integer	:: n_del_orb	


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`mpi_master`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`dim_list_core_orb`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`n_core_inact_act_orb`
       * :c:data:`n_core_orb_allocate`
       * :c:data:`n_inact_act_orb`
       * :c:data:`n_inact_orb_allocate`
       * :c:data:`n_virt_orb_allocate`
       * :c:data:`pt2_f`

 
.. c:var:: n_virt_orb_allocate


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer	:: n_virt_orb_allocate	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`


 
.. c:var:: ref_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: ref_bitmask	(N_int,2)


    Reference bit mask, used in Slater rules, chosen as Hartree-Fock bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`closed_shell_ref_bitmask`
       * :c:data:`coef_hf_selector`
       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`psi_det_hii`
       * :c:data:`psi_selectors_diag_h_mat`
       * :c:data:`ref_bitmask_energy`
       * :c:data:`ref_closed_shell_bitmask`

 
.. c:var:: reunion_of_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: reunion_of_bitmask	(N_int,2)


    Reunion of the inactive, active and virtual bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`list_inact`
       * :c:data:`n_int`


 
.. c:var:: reunion_of_cas_inact_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: reunion_of_cas_inact_bitmask	(N_int,2)


    Reunion of the inactive, active and virtual bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_int`


 
.. c:var:: reunion_of_core_inact_act_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: reunion_of_core_inact_act_bitmask	(N_int,2)


    Reunion of the core, inactive and active bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_int`
       * :c:data:`reunion_of_core_inact_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_inact_act_bitmask_4`
       * :c:data:`list_core_inact_act`

 
.. c:var:: reunion_of_core_inact_bitmask


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: reunion_of_core_inact_bitmask	(N_int,2)


    Reunion of the core and inactive and virtual bitmasks

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_inact_orb`
       * :c:data:`reunion_of_core_inact_act_bitmask`

 
.. c:var:: unpaired_alpha_electrons


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: unpaired_alpha_electrons	(N_int)


    Bitmask reprenting the unpaired alpha electrons in the HF_bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_int`


 
.. c:var:: virt_bitmask


    File : :file:`bitmask/core_inact_act_virt.irp.f`

    .. code:: fortran

        integer, allocatable	:: list_inact	(dim_list_inact_orb)
        integer, allocatable	:: list_virt	(dim_list_virt_orb)
        integer, allocatable	:: list_inact_reverse	(mo_num)
        integer, allocatable	:: list_virt_reverse	(mo_num)
        integer, allocatable	:: list_del_reverse	(mo_num)
        integer, allocatable	:: list_del	(mo_num)
        integer, allocatable	:: list_core	(dim_list_core_orb)
        integer, allocatable	:: list_core_reverse	(mo_num)
        integer, allocatable	:: list_act	(dim_list_act_orb)
        integer, allocatable	:: list_act_reverse	(mo_num)
        integer(bit_kind), allocatable	:: core_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: inact_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: act_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: virt_bitmask	(N_int,2)
        integer(bit_kind), allocatable	:: del_bitmask	(N_int,2)


    inact_bitmask : Bitmask of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    n_inact_orb   : Number of inactive orbitals
    virt_bitmask  : Bitmaks of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    n_virt_orb    : Number of virtual orbitals
    list_inact : List of the inactive orbitals which are supposed to be doubly excited
    in post CAS methods
    list_virt  : List of vritual orbitals which are supposed to be recieve electrons
    in post CAS methods
    list_inact_reverse : reverse list of inactive orbitals
    list_inact_reverse(i) = 0 ::> not an inactive
    list_inact_reverse(i) = k ::> IS the kth inactive
    list_virt_reverse : reverse list of virtual orbitals
    list_virt_reverse(i) = 0 ::> not an virtual
    list_virt_reverse(i) = k ::> IS the kth virtual
    list_act(i) = index of the ith active orbital
    
    list_act_reverse : reverse list of active orbitals
    list_act_reverse(i) = 0 ::> not an active
    list_act_reverse(i) = k ::> IS the kth active orbital

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dim_list_core_orb`
       * :c:data:`mo_class`
       * :c:data:`mo_num`
       * :c:data:`n_core_orb`
       * :c:data:`n_int`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`core_energy`
       * :c:data:`core_energy_erf`
       * :c:data:`core_fock_operator`
       * :c:data:`core_fock_operator_erf`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_mo`
       * :c:data:`inact_virt_bitmask`
       * :c:data:`list_core_inact_act`
       * :c:data:`list_inact_act`
       * :c:data:`mo_two_e_integrals_in_map`
       * :c:data:`mo_two_e_integrals_vv_from_ao`
       * :c:data:`reunion_of_bitmask`
       * :c:data:`reunion_of_cas_inact_bitmask`
       * :c:data:`reunion_of_core_inact_act_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`virt_bitmask_4`

 
.. c:var:: virt_bitmask_4


    File : :file:`bitmask/bitmasks.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: virt_bitmask_4	(N_int,4)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_int`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: bitstring_to_hexa:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine bitstring_to_hexa( output, string, Nint )


    Transform a bit string to a string in hexadecimal format for printing

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`
       * :c:func:`debug_spindet`

 
.. c:function:: bitstring_to_list:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine bitstring_to_list( string, list, n_elements, Nint)


    Gives the inidices(+1) of the bits set to 1 in the bit string

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map`
       * :c:func:`add_integrals_to_map_erf`
       * :c:func:`add_integrals_to_map_no_exit_34`
       * :c:func:`add_integrals_to_map_three_indices`
       * :c:func:`create_microlist`
       * :c:func:`example_bitmask`
       * :c:func:`getmobiles`
       * :c:data:`list_core_inact_act`
       * :c:data:`ref_bitmask_energy`

 
.. c:function:: bitstring_to_str:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine bitstring_to_str( output, string, Nint )


    Transform a bit string to a string for printing

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`add_integrals_to_map`
       * :c:func:`add_integrals_to_map_erf`
       * :c:func:`add_integrals_to_map_three_indices`
       * :c:func:`example_bitmask`
       * :c:func:`print_det`
       * :c:func:`print_spindet`

 
.. c:function:: broadcast_chunks_bit_kind:


    File : :file:`bitmask/mpi.irp.f`

    .. code:: fortran

        subroutine broadcast_chunks_bit_kind(A, LDA)


    Broadcast with chunks of ~2GB

 
.. c:function:: clear_bit_to_integer:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine clear_bit_to_integer(i_physical,key,Nint)


    set to 0 the bit number i_physical in the bitstring key

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`example_bitmask`
       * :c:data:`ref_closed_shell_bitmask`

 
.. c:function:: debug_det:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine debug_det(string,Nint)


    Subroutine to print the content of a determinant in '+-' notation and
    hexadecimal representation.

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`build_fock_tmp`
       * :c:func:`example_determinants`
       * :c:func:`get_excitation_degree_vector_single_or_exchange_verbose`
       * :c:func:`number_of_holes_verbose`
       * :c:func:`number_of_particles_verbose`
       * :c:func:`routine_example_psi_det`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_hexa`
       * :c:func:`print_det`

 
.. c:function:: debug_spindet:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine debug_spindet(string,Nint)


    Subroutine to print the content of a determinant in '+-' notation and
    hexadecimal representation.

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_hexa`
       * :c:func:`print_spindet`

 
.. c:function:: example_bitmask:


    File : :file:`bitmask/example.irp.f`

    .. code:: fortran

        subroutine example_bitmask


    subroutine that illustrates the main features available in bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_core_orb`
       * :c:data:`list_inact`
       * :c:data:`n_int`
       * :c:data:`mo_num`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_list`
       * :c:func:`bitstring_to_str`
       * :c:func:`clear_bit_to_integer`
       * :c:func:`set_bit_to_integer`

 
.. c:function:: initialize_bitmask_to_restart_ones:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine initialize_bitmask_to_restart_ones


    Initialization of the generators_bitmask to the restart bitmask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask_restart`
       * :c:data:`n_generators_bitmask`
       * :c:data:`generators_bitmask`
       * :c:data:`n_int`

 
.. c:function:: is_a_1h:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_1h(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_1h1p:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_1h1p(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_1h2p:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_1h2p(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_1p:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_1p(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_2h:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_2h(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_2h1p:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_2h1p(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_2p:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_2p(key_in)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_a_two_holes_two_particles:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_a_two_holes_two_particles(key_in)


    logical function that returns True if the determinant 'key_in'
    belongs to the 2h-2p excitation class of the DDCI space
    this is calculated using the CAS_bitmask that defines the active
    orbital space, the inact_bitmasl that defines the inactive oribital space
    and the virt_bitmask that defines the virtual orbital space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`list_inact`
       * :c:data:`n_int`

 
.. c:function:: is_i_in_virtual:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        logical function is_i_in_virtual(i)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`list_inact`
       * :c:data:`n_int`

 
.. c:function:: is_the_hole_in_det:


    File : :file:`bitmask/find_hole.irp.f`

    .. code:: fortran

        logical function is_the_hole_in_det(key_in,ispin,i_hole)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: is_the_particl_in_det:


    File : :file:`bitmask/find_hole.irp.f`

    .. code:: fortran

        logical function is_the_particl_in_det(key_in,ispin,i_particl)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_int`

 
.. c:function:: list_to_bitstring:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine list_to_bitstring( string, list, n_elements, Nint)


    Returns the physical string "string(N_int,2)" from the array of
    occupations "list(N_int*bit_kind_size,2)

    Called by:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`list_inact`

 
.. c:function:: modify_bitmasks_for_hole:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine modify_bitmasks_for_hole(i_hole)


    modify the generators_bitmask in order that one can only excite
    the electrons occupying i_hole

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_generators_bitmask`
       * :c:data:`generators_bitmask`
       * :c:data:`n_int`
       * :c:data:`index_holes_bitmask`

 
.. c:function:: modify_bitmasks_for_hole_in_out:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine modify_bitmasks_for_hole_in_out(i_hole)


    modify the generators_bitmask in order that one can only excite
    the electrons occupying i_hole

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_generators_bitmask`
       * :c:data:`generators_bitmask`
       * :c:data:`index_holes_bitmask`

 
.. c:function:: modify_bitmasks_for_particl:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine modify_bitmasks_for_particl(i_part)


    modify the generators_bitmask in order that one can only excite
    the electrons to the orbital i_part

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`index_particl_bitmask`
       * :c:data:`n_generators_bitmask`
       * :c:data:`generators_bitmask`
       * :c:data:`n_int`

 
.. c:function:: number_of_holes:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        integer function number_of_holes(key_in)


    Function that returns the number of holes in the inact space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`n_int`

 
.. c:function:: number_of_holes_verbose:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        integer function number_of_holes_verbose(key_in)


    function that returns the number of holes in the inact space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`reunion_of_core_inact_bitmask`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`

 
.. c:function:: number_of_particles:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        integer function number_of_particles(key_in)


    function that returns the number of particles in the virtual space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`list_inact`
       * :c:data:`n_int`

 
.. c:function:: number_of_particles_verbose:


    File : :file:`bitmask/bitmask_cas_routines.irp.f`

    .. code:: fortran

        integer function number_of_particles_verbose(key_in)


    function that returns the number of particles in the inact space

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`cas_bitmask`
       * :c:data:`list_inact`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`

 
.. c:function:: print_det:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine print_det(string,Nint)


    Subroutine to print the content of a determinant using the '+-' notation

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_det`
       * :c:func:`example_determinants`
       * :c:func:`print_generators_bitmasks_holes`
       * :c:func:`print_generators_bitmasks_holes_for_one_generator`
       * :c:func:`print_generators_bitmasks_particles`
       * :c:func:`print_generators_bitmasks_particles_for_one_generator`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_str`

 
.. c:function:: print_generators_bitmasks_holes:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine print_generators_bitmasks_holes



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`n_int`
       * :c:data:`index_holes_bitmask`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_det`

 
.. c:function:: print_generators_bitmasks_holes_for_one_generator:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine print_generators_bitmasks_holes_for_one_generator(i_gen)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`
       * :c:data:`n_int`
       * :c:data:`index_holes_bitmask`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_det`

 
.. c:function:: print_generators_bitmasks_particles:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine print_generators_bitmasks_particles



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`index_particl_bitmask`
       * :c:data:`generators_bitmask`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_det`

 
.. c:function:: print_generators_bitmasks_particles_for_one_generator:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine print_generators_bitmasks_particles_for_one_generator(i_gen)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`index_particl_bitmask`
       * :c:data:`generators_bitmask`
       * :c:data:`n_int`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_det`

 
.. c:function:: print_spindet:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine print_spindet(string,Nint)


    Subroutine to print the content of a determinant using the '+-' notation

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`debug_spindet`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`bitstring_to_str`

 
.. c:function:: set_bit_to_integer:


    File : :file:`bitmask/bitmasks_routines.irp.f`

    .. code:: fortran

        subroutine set_bit_to_integer(i_physical,key,Nint)


    set to 1 the bit number i_physical in the bitstring key

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`example_bitmask`

 
.. c:function:: set_bitmask_hole_as_input:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine set_bitmask_hole_as_input(input_bimask)


    set the generators_bitmask for the holes
    as the input_bimask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_generators_bitmask`
       * :c:data:`generators_bitmask`
       * :c:data:`n_int`
       * :c:data:`index_holes_bitmask`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`

 
.. c:function:: set_bitmask_particl_as_input:


    File : :file:`bitmask/modify_bitmasks.irp.f`

    .. code:: fortran

        subroutine set_bitmask_particl_as_input(input_bimask)


    set the generators_bitmask for the particles
    as the input_bimask

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`index_particl_bitmask`
       * :c:data:`n_generators_bitmask`
       * :c:data:`generators_bitmask`
       * :c:data:`n_int`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`generators_bitmask`

