---
title: "Installation"
date: Tue Feb  3 16:39:50 CET 2019
draft: false
---


---------------------


{{< youtube Un69atP2-30 >}}

1. Download the source files

```bash
git clone https://github.com/QuantumPackage/qp2.git
```

2. Install the dependencies

```bash
./configure   # tells you what to install. See in INSTALL.rst
```

3. *USE AT YOUR OWN RISK, NO SUPPORT WILL BE PROVIDED*

```bash
./configure --install something
```

4. The following libraries are needed for the ocaml package

  * zlib1g-dev
  * libncurses5-dev
  * pkg-config 
  * libgmp3-dev 
  * m4

5. Once the smiling cow appears, load the environment variables

```bash
source quantum_package.rc
```

6. Define your compilation configuration file

```bash
./configure -c ./config/a_file.cfg
```

7. Compile with Ninja

```bash
ninja
```



