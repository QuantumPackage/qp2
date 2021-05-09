.. _module_electrons: 
 
.. program:: electrons 
 
.. default-role:: option 
 
=========
electrons
=========

Describes the electrons. For the moment, only the number of alpha
and beta electrons are provided by this module.


Assumptions
===========

* `elec_num` >= 0
* `elec_alpha_num` >= 0
* `elec_beta_num` >= 0
* `elec_alpha_num` >= `elec_beta_num`


 
 
 
EZFIO parameters 
---------------- 
 
.. option:: elec_alpha_num
 
    Numbers of electrons alpha ("up")
 
 
.. option:: elec_beta_num
 
    Numbers of electrons beta ("down")
 
 
.. option:: elec_num
 
    Numbers total of electrons (alpha + beta)
 
    Default: = electrons.elec_alpha_num + electrons.elec_beta_num
 
 
Providers 
--------- 
 
.. c:var:: elec_num


    File : :file:`electrons/electrons.irp.f`

    .. code:: fortran

        integer	:: elec_num	
        integer, allocatable	:: elec_num_tab	(2)


    Numbers of alpha ("up") , beta ("down") and total electrons

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`ezfio_filename`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`psi_det_hii`

 
.. c:var:: elec_num_tab


    File : :file:`electrons/electrons.irp.f`

    .. code:: fortran

        integer	:: elec_num	
        integer, allocatable	:: elec_num_tab	(2)


    Numbers of alpha ("up") , beta ("down") and total electrons

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`ezfio_filename`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`diagonal_h_matrix_on_psi_det`
       * :c:data:`psi_det_hii`

