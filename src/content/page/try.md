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

**You may need to tell your web browser you accept to load insecure scripts to
see the terminal**.

<iframe id="shellframe" src="https://lcpq.ups-tlse.fr/qpterm/" width="800" height="600" frameBorder="1" scrolling="no">Browser not compatible.</iframe>
</body>


Here is an example of a few commands you can run to
get the Full-CI energy of the HCN molecule.

First create a file named `be.zmt` containing the z-matrix of a Beryllium atom.

``` 
echo be > be.zmt
```

Create the EZFIO database as follows:

```
qp create_ezfio -b cc-pvtz be.zmt -o be
```

Run a Hartree-Fock calculation:

```
qp run scf | tee scf.out
```

The MOs are saved in the EZFIO database. 

Now run the CIPSI calculation in the full CI space:

```
qp run fci | tee fci.out
```

That's it!



