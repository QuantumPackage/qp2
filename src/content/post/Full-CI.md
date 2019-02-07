---
title: "Full-CI calculation"
date: Tue Feb  5 16:39:50 CET 2019
draft: false
---


---------------------



{{< youtube 4nmdCAPkZlc >}}


1. Enter in qpsh mode for auto completion

```bash
~/your_path_to_qp/bin/qpsh
```

2. Create an EZFIO DATA BASE from .xyz

```bash
qp create_ezfio -b 6-31g hcn.xyz
# output:: hcn.ezfio DATA BASE
```

3. Run a HF calculation on the EZFIO

```bash
qp run scf | tee output_file.scf.out
# output:: create MO basis in the EZFIO
```

4. Freeze the *1s* orbitals in the EZFIO

```bash
qp set_frozen_core
```

5. Run the CIPSI algorithm on the EZFIO

```bash
qp run fci | tee output_file.fci.out
```

6. Check the CONVERGENCE of the ENERGY

```bash
qp_e_conv_fci hcn.ezfio
# output:: hcn.ezfio.1.conv and pdf
```




