# CIS Overlap


## Purpose of Module

The main purpose of the module is calculating overlap matrices between sets of wave functions calculated at different geometries, using different atomic basis sets or different electronic structure methods. Wave functions are read from single point energy calculations performed by a quantum chemistry code. Geometries, basis sets, molecular orbital coefficients and wave function coefficients need to be read to calculate the overlaps. Only wave functions (approximately) expanded in CIS form are presently supported and overlaps between them are very efficiently calculated either by expanding the wave functions in terms of excitations between natural transition orbitals (NTOs) or by expanding the overlap determinants into level 2 minors (L2M).

The algorithms implemented in the code are described in the following references:

 * &nbsp;M. Sapunar, T. Piteša, D. Davidović and N. Došlić, Highly Efficient Algorithms for CIS Type Excited State Wave Function Overlaps, [*J. Chem. Theory Comput.* (2019)][NTOpaper]
 * &nbsp;P. Alonso-Jordá, D. Davidović, M. Sapunar, J. R.Herrero and E. S. Quintana-Ortí., Efficient update of determinants for many-electron wave function overlaps, [*Comput. Phys. Commun.* (2021)][L2Mpaper]
 
For wave function overlap calculations, the code is currently interfaced only to Turbomole, but additional interfaces can easily be added.


## Installation and technical information

* Language
  * Fortran 2008, C
* License
  * MIT license (MIT)
* Prerequisites
  * cmake
  * make
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

after running  make, the executables will be in the bin directory. The main executable is cis_overlap.exe.
For instructions/options run the executables with the \--help command line argument.

## Testing

Some tests are available by running test.py in the 'test/test_cases' directory. 

## Contributors

The code has been developed by (listed alphabetically):

* Pedro Alonso-Jordá
* Davor Davidović
* Tomislav Piteša
* Enrique S. Quintana-Ortí
* Marin Sapunar


## Source Code

The source code is available at: [GitHub][Git]


[Git]: https://github.com/marin-sapunar/cis_nto
[F95]: https://software.intel.com/en-us/mkl-linux-developer-guide-fortran-95-interfaces-to-lapack-and-blas
[NTOpaper]: https://pubs.acs.org/doi/10.1021/acs.jctc.9b00235
[L2Mpaper]: https://dx.doi.org/10.1016/j.cpc.2020.107521


