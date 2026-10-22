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
