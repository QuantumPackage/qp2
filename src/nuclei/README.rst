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
