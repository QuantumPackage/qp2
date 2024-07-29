=========
gpu_intel
=========

Intel implementation of GPU routines. Uses MKL and SYCL.
```bash
icpx -fsycl gpu.cxx -c -qmkl=sequential
```
