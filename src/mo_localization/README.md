# mo_localization

Some parameters can be changed with qp edit in the mo_localization section
(cf below). Similarly for the trust region parameters in the 
utils_trust_region section. The localization without the trust region 
is not available for the moment.  

The irf.f files can be generated from the org ones using emacs. 
If you modify the .org files, don't forget to do (you need emacs):  
``` 
./TANGLE_org_mode.sh  
ninja  
```  

# Orbital localisation
To localize the MOs:  
```
qp run localization  
```
After that the ezfio directory contains the localized MOs  
 
But to do so the mo_class must be defined before, run 
```
qp set_mo_class -q
```
for more information or  
```
qp set_mo_class -c [] -a [] -v [] -i [] -d [] 
```
to set the mo classes. We don't care about the name of the   
mo classes. The algorithm just localizes all the MOs of  
a given class between them, for all the classes, except the deleted MOs.  

If you just on kind of mo class to localize all the MOs between them  
you have to put:
```
qp set mo_localization security_mo_class false
```

Before the localization, a kick is done for each mo class  
(except the deleted ones) to break the MOs. This is done by   
doing a rotation between the MOs.
This feature can be removed by setting:
```
qp set mo_localization kick_in_mos false
```
and the default angle for the rotation can be changed with:
```
qp set mo_localization angle_pre_rot 1e-3 # or something else
```

After the localization, the MOs of each class (except the deleted ones)  
can be sorted between them using the diagonal elements of  
the fock matrix with:
```
qp set mo_localization sort_mos_by_e true
```

You can check the Hartree-Fock energy before/during/after the localization  
by putting (only for debugging):
```
qp set mo_localization debug_hf true 
```

## Foster-Boys & Pipek-Mezey
Foster-Boys:  
``` 
qp set mo_localization localization_method boys 
``` 
 
Pipek-Mezey:  
``` 
qp set mo_localization localization_method pipek 
``` 

# Break the spatial symmetry of the MOs
To break the spatial symmetry of the MOs:   
```
qp run break_spatial_sym
```
The default angle for the rotations is too big for this kind of
application, a value between 1e-3 and 1e-6 should break the spatial
symmetry with just a small change in the energy:
```
qp set mo_localization angle_pre_rot 1e-3
``` 

# With or without hessian + trust region
With hessian +  trust region
```
qp set mo_localization localisation_use_hessian true
```
It uses the trust region algorithm with the diagonal of the hessian of the
localization criterion with respect to the MO rotations.

Without the hessian and the trust region
```
qp set mo_localization localisation_use_hessian false
```
By doing so it does not require to store the hessian but the
convergence is not easy, in particular for virtual MOs.
It seems that it not possible to converge with Pipek-Mezey
localization with this approach.

# Further improvements: 
- Cleaner repo 
- Correction of the errors in the documentations 
- option with/without trust region 
