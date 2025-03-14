source ~/qp2/quantum_package.rc
## Example of how to generate an additional h2o molecule, stored as a extra basis/nuclei etc .. to an He

sys_B=h2o.xyz
basis_B=sto-3g
output_B=${sys_B%.xyz}_${basis_B}

sys_A=He_A.xyz
basis_A=cc-pvtz
output_A=${sys_A%.xyz}_${basis_A}_extra_${output_B}

# we create the system "B" that will be attached as an "extra system" to the syste "A"
qp create_ezfio -b $basis_B $sys_B -o ${output_B}
# we perform an HF calculation to obtain the AO density matrix 
qp run scf 
# we save the density matrix in the EZFIO
qp run save_one_e_dm 
# we create the system "A" 
qp create_ezfio -b $basis_A $sys_A -o ${output_A}
# We perform an SCF calculation
qp run scf 
# we copy the system "B" information as extra nuclei/basis etc in the EZFIO of system "A"
qp_copy_extra_basis ${output_B} ${output_A}

# we execute an example of progra that prints a lot of useful integrals/information on the A-B interaction
qp run test_extra_basis | tee ${output_A}.test_extra_basis
