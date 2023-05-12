# Orbital localisation
To localize the MOs:  
```
qp run localization  
```
By default, the different otbital classes are automatically set by splitting  
the orbitales in the following classes:  
- Core -> Core  
- Active, doubly occupied -> Inactive  
- Active, singly occupied -> Active  
- Active, empty -> Virtual  
- Deleted -> Deleted  
The orbitals will be localized among each class, excpect the deleted ones.
If you want to choose another splitting, you can set
```
qp set mo_localization auto_mo_class false
```
and define the classes with
```
qp set_mo_class -c [] -a [] -v [] -i [] -d [] 
```
for more information
```
qp set_mo_class -q
```
We don't care about the name of the   
mo classes. The algorithm just localizes all the MOs of  
a given class between them, for all the classes, except the deleted MOs.  
If you are using the last option don't forget to reset the initial mo classes  
after the localization.

Before the localization, a kick is done for each mo class  
(except the deleted ones) to break the MOs. This is done by   
doing a given rotation between the MOs.
This feature can be removed by setting:
```
qp set localization kick_in_mos false
```
and the default angle for the rotation can be changed with:
```
qp set localization angle_pre_rot 1e-3 # or something else
```

After the localization, the MOs of each class (except the deleted ones)  
can be sorted between them using the diagonal elements of  
the fock matrix with:
```
qp set localization sort_mos_by_e true
```

You can check the Hartree-Fock energy before/during/after the localization  
by putting (only for debugging):
```
qp set localization debug_hf true 
```

## Foster-Boys & Pipek-Mezey
Foster-Boys:  
``` 
qp set localization localization_method boys 
``` 
 
Pipek-Mezey:  
``` 
qp set localization localization_method pipek 
``` 

# Break the spatial symmetry of the MOs
This program work exactly as the localization.  
To break the spatial symmetry of the MOs:   
```
qp run break_spatial_sym
```
The default angle for the rotations is too big for this kind of
application, a value between 1e-3 and 1e-6 should break the spatial
symmetry with just a small change in the energy:
```
qp set localization angle_pre_rot 1e-3
``` 

# With or without hessian + trust region
With hessian +  trust region
```
qp set localization localisation_use_hessian true
```
It uses the trust region algorithm with the diagonal of the hessian of the 
localization criterion with respect to the MO rotations.   

Without the hessian and the trust region
```
qp set localization localisation_use_hessian false
```
By doing so it does not require to store the hessian but the
convergence is not easy, in particular for virtual MOs.
It seems that it not possible to converge with Pipek-Mezey
localization with this approach.

# Parameters
Some other parameters are available for the localization (qp edit for more details).

# Tests
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
 
