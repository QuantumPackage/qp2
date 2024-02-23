# This is an example for MULTI STATE CALCULATION STATE AVERAGE CASSCF 
# We will compute 3 states on the O2 molecule 
# The Ground state and 2 degenerate excited states 
# Please follow carefully the tuto :) 

##### PREPARING THE EZFIO 
# Set the path to your QP2 directory 
QP_ROOT=my_fancy_path 
source ${QP_ROOT}/quantum_package.rc 
# Create the EZFIO folder
qp create_ezfio -b aug-cc-pvdz O2.xyz -m 3 -a -o O2_avdz_multi_state
# Start with ROHF orbitals 
qp run scf # ROHF energy : -149.619992871398
# Freeze the 1s orbitals of the two oxygen
qp set_frozen_core 

##### PREPARING THE ORBITALS WITH NATURAL ORBITALS OF A CIS 
# Tell that you want 3 states in your WF
qp set determinants n_states 3 
# Run a CIS wave function to start your calculation 
qp run cis | tee ${EZFIO_FILE}.cis_3_states.out # -149.6652601409258 -149.4714726176746 -149.4686165431939
# Save the STATE AVERAGE natural orbitals for having a balanced description 
# This will also order the orbitals according to their occupation number 
# Which makes the active space selection easyer !
qp run save_natorb | tee ${EZFIO_FILE}.natorb_3states.out 

##### PREPARING A CIS GUESS WITHIN THE ACTIVE SPACE 
# Set an active space which has the most of important excitations 
# and that maintains symmetry : the ACTIVE ORBITALS are from """6 to 13"""

# YOU FIRST FREEZE THE VIRTUALS THAT ARE NOT IN THE ACTIVE SPACE 
# !!!!! WE SET TO "-D" for DELETED !!!!
qp set_mo_class -c "[1-5]" -a "[6-13]" -d "[14-46]" 
# You create a guess of CIS type WITHIN THE ACTIVE SPACE 
qp run cis | tee ${EZFIO_FILE}.cis_3_states_active_space.out # -149.6515472533511 -149.4622878024821 -149.4622878024817
# You tell to read the WFT stored (i.e. the guess we just created) 
qp set determinants read_wf True

##### DOING THE CASSCF 
### SETTING PROPERLY THE ACTIVE SPACE FOR CASSCF 
# You set the active space WITH THE VIRTUAL ORBITALS !!!
# !!!!! NOW WE SET TO "-v" for VIRTUALS  !!!!!
qp set_mo_class -c "[1-5]" -a "[6-13]" -v "[14-46]" 

# You tell that it is a small actice space so the CIPSI can take all Slater determinants 
qp set casscf_cipsi small_active_space True 
# You specify the output file 
output=${EZFIO_FILE}.casscf_3states.out 
# You run the CASSCF calculation 
qp run casscf  | tee ${output} # -149.7175867510 -149.5059010227 -149.5059010226

# Some grep in order to get some numbers useful to check convergence 
# State average energy 
grep "State-average CAS-SCF energy =" $output | cut -d "=" -f 2 > data_e_average
# Delta E anticipated for State-average energy, only usefull to check convergence
grep "Predicted energy improvement =" $output | cut -d "=" -f 2 > data_improve
# Ground state energy
grep "state  1 E + PT2 energy"        $output | cut -d "=" -f 2 > data_1
# First excited state energy
grep "state  2 E + PT2 energy"        $output | cut -d "=" -f 2 > data_2
# First excitation energy 
grep "state  2 Delta E+PT2"           $output | cut -d "=" -f 2 > data_delta_E2
# Second excited state energy
grep "state  3 E + PT2 energy"        $output | cut -d "=" -f 2  > data_3
# Second excitation energy 
grep "state  3 Delta E+PT2"           $output | cut -d "=" -f 2 > data_delta_E3
