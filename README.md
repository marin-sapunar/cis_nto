# CIS NTO Module


## Purpose of Module

The main purpose of the module is calculating overlap matrices between sets of wave functions calculated at different geometries, using different atomic basis sets or different electronic structure methods. Wave functions are read from single point energy calculations performed by a quantum chemistry code. Geometries, basis sets, molecular orbital coefficients and wave function coefficients need to be read to calculate the overlaps. Only wave functions (approximately) expanded in CIS form are presently supported and overlaps between them are very efficiently calculated either by expanding the wave functions in terms of excitations between natural transition orbitals (NTOs) or by expanding the overlap determinants into level 2 minors (L2M).

If using the code, please cite the following reference:
 * M. Sapunar, T. Piteša, D. Davidović and N. Došlić, [J. Chem. Theory Comput. (2019.)][NTOpaper]
 
which describes the main methods implemented in the code.

For wave function overlap calculations, the code is currently interfaced only to Turbomole, but additional interfaces can easily be added.


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
    * MKL (with the [MKL Fortran 95 Interface][F95])

Commands to download and compile the code:

```
git clone https://github.com/marin-sapunar/cis_nto.git
cd cis_nto
mkdir build
cd build
cmake ..
make
```

after running  make, the executables will be in the bin directory.
For instructions/options run the executables with the --help command line argument.

## Testing

Some tests are available by running test.py in the 'test/test_cases' directory. 

## Source Code

The source code is available at: [GitHub][Git]


[Git]: https://github.com/msapunar/cis_nto
[F95]: https://software.intel.com/en-us/mkl-linux-developer-guide-fortran-95-interfaces-to-lapack-and-blas
[NTOpaper]: https://pubs.acs.org/doi/10.1021/acs.jctc.9b00235


