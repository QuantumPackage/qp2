---
title: "Video Tutorials"
date: Tue Feb  4 16:39:50 CET 2019
draft: false
---


---------------------


# Hartree-Fock calculation

{{< youtube KIl_xq-NRLo >}}

1. Enter in qpsh mode for auto completion

```bash
~/your_path_to_qp/bin/qpsh
```

2. Create an EZFIO DATA BASE from .xyz

```bash
qp create_ezfio -b 6-31g hcn.xyz
```

3. Run a HF calculation on the EZFIO  and put the output in a file

```bash
qp run scf | tee hcn.ezfio.scf.out
```




