.. _module_extra_nuclei: 
 
.. program:: extra_nuclei 
 
.. default-role:: option 
 
============
extra_nuclei
============

Extra nuclei added to the system, some of them fictitious.

The number, coordinates, charges and labels of the extra nuclei are
read from the |EZFIO| database, and the module distinguishes the real
ones, :c:data:`extra_nucl_real_list`, from the fictitious ones,
:c:data:`extra_nucl_fictious_list`. The derived quantities are the
pair distances :c:data:`extra_nucl_dist` and the center of mass
:c:data:`extra_center_of_mass`.

The program :ref:`extra_nuclei` is a small example of the use of
these providers.
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: extra_nucl_num
 
    Number of nuclei
 
 
.. option:: extra_nucl_pouet
 
    Number of nuclei
 
    Default: 1
 
.. option:: extra_nucl_label
 
    Nuclear labels
 
 
.. option:: extra_nucl_charge
 
    Nuclear charges
 
 
.. option:: extra_nucl_coord
 
    Nuclear coordinates in the format (:, {x,y,z})
 
 
.. option:: extra_nucl_real_num
 
    Number of real nuclei
 
 
.. option:: extra_nucl_fictious_num
 
    Number of fictious nuclei
 
 
.. option:: extra_nucl_real_fictious_list
 
    List of real nuclei to which fictious nuclei are attached to
 
 
.. option:: extra_nucl_fictious_list
 
    List of fictious nuclei
 
 
.. option:: extra_nucl_real_list
 
    List of real nuclei
 
 
 
Programs 
-------- 
 
 * :ref:`extra_nuclei` 
 
Providers 
--------- 
 
.. c:var:: double extra_center_of_mass


    File : :file:`nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: extra_center_of_mass	(3)


    Center of mass of the molecule

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`element_name`
       * :c:data:`extra_nucl_charge`
       * :c:data:`extra_nucl_coord`
       * :c:data:`extra_nucl_num`


 
.. c:var:: double extra_nucl_coord


    File : :file:`nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: extra_nucl_coord	(extra_nucl_num,3)


    Nuclear coordinates in the format (:, {x,y,z})

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`extra_nucl_charge`
       * :c:data:`extra_nucl_label`
       * :c:data:`extra_nucl_num`
       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`output_wall_time_0`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`extra_center_of_mass`
       * :c:data:`extra_nucl_coord_transp`
       * :c:data:`extra_nucl_dist`

 
.. c:var:: double extra_nucl_coord_transp


    File : :file:`nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: extra_nucl_coord_transp	(3,extra_nucl_num)


    Transposed array of extra_nucl_coord

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`extra_nucl_coord`
       * :c:data:`extra_nucl_num`


 
.. c:var:: double extra_nucl_dist


    File : :file:`nuclei.irp.f`

    .. code:: fortran

        double precision, allocatable	:: extra_nucl_dist	(extra_nucl_num,extra_nucl_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`extra_nucl_coord`
       * :c:data:`extra_nucl_num`


 
 
Subroutines / functions 
----------------------- 
