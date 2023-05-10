# Orbital optimization

## Methods  
Different methods are available:  
- full hessian  
``` 
qp set orbital_optimization optimization_method full  
```  
- diagonal hessian  
``` 
qp set orbital_optimization optimization_method diag  
``` 
- identity matrix  
``` 
qp set orbital_optimization optimization_method none  
``` 

After the optimization the ezfio contains the optimized orbitals
 
## For a fixed number of determinants
To optimize the MOs for the actual determinants:  
``` 
qp run orb_opt
``` 
 
## For a complete optimization, i.e, with a larger and larger wave function
To optimize the MOs with a larger and larger wave function:  
``` 
qp run optimization  
``` 

The results are stored in the EZFIO in "mo_optimization/result_opt",
with the following format:  
(1) (2) (3) (4)  
1: Number of determinants in the wf,  
2: Cispi energy before the optimization,   
3: Cipsi energy after the optimization,  
4: Energy difference between (2) and (3).  
 
The optimization process if the following: 
- we do a first cipsi step to obtain a small number of determinants in the wf 
- we run an orbital optimization for this wf 
- we do a new cipsi step to double the number of determinants in the wf 
- we run an orbital optimization for this wf 
- ... 
- we do that until the energy difference between (2) and (3) is  
  smaller than the targeted accuracy for the cispi (targeted_accuracy_cipsi in qp edit) 
  or the wf is larger than a given size (n_det_max_opt in qp_edit) 
- after that you can reset your determinants (qp reset -d) and run a clean Cispi calculation  
  
### End of the optimization
You can choos the number of determinants after what the 
optimization will stop:
```
qp set orbital_optimization n_det_max_opt 1e5 # or any number
```
## Weight of the states
You can change the weights of the differents states directly in qp edit.  
It will affect ths weights used in the orbital optimization.

# Tests
To run the tests:  
``` 
qp test
``` 

# Org files
The org files are stored in the directory org in order to avoid overwriting on user changes.
The org files can be modified, to export the change to the source code, run
```
./TANGLE_org_mode.sh
mv *.irp.f ../.
```

