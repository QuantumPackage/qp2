---
title: "QP terminal"
date: Wed Jan 23 22:35:00 CET 2019
draft: false
---


---------------------

You can try *Quantum Package* in the terminal below.
To configure the terminal for your favorite text editor,
set the ``EDITOR`` environment variable:

```
export EDITOR=vim
```


Here is an example of a few commands you can run to
get the Full-CI energy of the HCN molecule.

First create a file named `hcn.xyz` containing the *xyz* coordinates.

``` 
$ cat << EOF > hcn.xyz
3
HCN molecule
C    0.0    0.0    0.0
H    0.0    0.0    1.064
N    0.0    0.0    -1.156
EOF
```

Create the EZFIO database as follows:

```
qp create_ezfio -b 6-31g hcn.xyz -o hcn
```

Run a Hartree-Fock calculation:

```
qp run SCF | tee scf.out
```

The MOs are saved in the EZFIO database. Now freeze the core electrons:

```
qp set_frozen_core
```

And run the CIPSI calculation in the valence full CI space:

```
qp run fci | tee fci.out
```

That's it!

<iframe id="shellframe" src="http://localhost:8080" width="100%" height="500" frameBorder="0" scrolling="no">Browser not compatible.</iframe>
</body>




