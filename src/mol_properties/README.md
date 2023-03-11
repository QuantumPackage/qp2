# Molecular properties

Available quantities:  
- Electric dipole moment  
- Electric transition dipole moment  
- Oscillator strength  

They are not computed by default. To compute them:
```
qp set mol_properties calc_dipole_moment true  
qp set mol_properties calc_tr_dipole_moment true  
qp set mol_properties calc_osc_str true  
```
If you are interested in transitions between two excited states:  
```
qp set mol_properties print_all_transitions true
```
They can be obtained by running
```
qp run properties
```
or at each step of a cipsi calculation with
```
qp run fci
```
