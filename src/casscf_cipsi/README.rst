======
casscf
======

|CASSCF| program with the CIPSI algorithm.

Example of inputs
-----------------

a) Small active space : standard CASSCF 
---------------------------------------
Let's do O2 (triplet) in aug-cc-pvdz with the following geometry (xyz format, Bohr units)
3

 O           0.0000000000        0.0000000000       -1.1408000000
 O           0.0000000000        0.0000000000        1.1408000000

# Create the ezfio folder 
qp create_ezfio -b aug-cc-pvdz O2.xyz -m 3 -a -o O2_avdz

# Start with an ROHF guess 
qp run scf | tee ${EZFIO_FILE}.rohf.out

# Get the ROHF energy for check 
qp get hartree_fock energy # should be -149.4684509

# Define the full valence active space: the two 1s are doubly occupied, the other 8 valence orbitals are active 
# CASSCF(12e,10orb) 
qp set_mo_class -c "[1-2]" -a "[3-10]" -v "[11-46]"

# Specify that you want an near exact CASSCF, i.e. the CIPSI selection will stop at pt2_max = 10^-10
qp set casscf_cipsi small_active_space True 
# RUN THE CASSCF 
qp run casscf | tee ${EZFIO_FILE}.casscf.out
# you should find around -149.7243542


b) Large active space : Exploit the selected CI in the active space 
-------------------------------------------------------------------
#Let us start from the small active space calculation orbitals and add another 10 virtuals: CASSCF(12e,20orb)
qp set_mo_class -c "[1-2]" -a "[3-20]" -v "[21-46]"
# As this active space is larger, you unset the small_active_space feature 
qp set casscf_cipsi small_active_space False
# As it is a large active space, the energy convergence thereshold is set to be 0.0001
qp run casscf | tee ${EZFIO_FILE}.casscf_large.out
# you should find around -149.9046

