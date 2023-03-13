# CCSD in spin orbitals and spatial orbitals

CCSD and CCSD(T) in spin orbitals for open and closed shell systems.
CCSD and CCSD(T) in spatial orbitals for closed shell systems.

## Calculations
The program will automatically choose the version in spin or spatial orbitals
To run the general program:  
```
qp run ccsd
```
Nevertheless, you can enforce the run in spin orbitals with
```
qp run ccsd_spin_orb
```

## Settings
The settings can be changed with:
```
qp set utils_cc cc_#param #val
```
For more informations on the settings, look at the module utils_cc and its documentation. 

## Org files
The org files are stored in the directory org in order to avoid overwriting on user changes.
The org files can be modified, to export the change to the source code, run
```
./TANGLE_org_mode.sh and
mv *.irp.f ../.
```

