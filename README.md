# CIS NTO Module


## Purpose of Module

The CIS NTO module contains programs for calculating the natural transition orbitals for CIS type wave functions and using them for calculating the overlaps or Dyson orbitals between two sets of CIS wave functions.

## Installation and technical information

* Language
  * Fortran 2008
* License
  * MIT license (MIT)
* Prerequisites
  * cmake
  * make
  * BLAS/LAPACK
    * OpenBLAS
    * MKL (with the [MKL Fortran 95 Interface](https://software.intel.com/en-us/mkl-linux-developer-guide-fortran-95-interfaces-to-lapack-and-blas))

Commands to downolad and install the module:

```
git clone --recursive https://github.com/marin-sapunar/cis_nto.git
cd cis_nto
mkdir build
cd build
cmake ..
make
```

after running  make, the executables will be in the bin directory. 
For instructions/options run the executables with the --help command line argument.


## Testing


## Source Code

The source code is available at: [GitHub][Git]


[Git]: https://github.com/msapunar/cis_nto

