# AMMeltpoolThermoFluidics 

### About
AMMeltpoolThermoFluidics is a C++ code developed to model heat conduction and fluid flow in a complex laser powder bed fusion process in Additive Manufacturing. The NonDimThermFluid has a finite-element code that models the relevant physics of the metlpool and has input parameters in terms of dimensionless numbers. The other branches have code suitable for different simplified problems which were used in the development of the model. The salient features of the computational implementation are: high-fidelity MPI-based parallel implementation, support for iterative solvers with preconditioning, and fully implicit time-stepping schemes. 

### Installation:

AMMeltpoolThermoFluidics code builds on top of the deal.II library.

1) Install CMake, PETSc, Trilinos, SLEPc, p4est, and deal.II (version 9.3.0 recommended)<br>

2) Clone the AMMeltpoolThermoFluidics GitHub repository <br>
```
$ git clone https://github.com/cmmg/AMMeltpoolThermoFluidics.git
$ cd AMMeltpoolThermoFluidics
$ git checkout master
$ cmake .
$ make -j nprocs
  ```
[here nprocs denotes the number of processors]

### Visualization:

  Output of the primary fields is in the vtk format (parallel:*.pvtu, serial:*.vtu files). These can be visualized with the following open-source applications:
  1. VisIt (https://visit-dav.github.io/visit-website/releases-as-tables/)
  2. Paraview (http://www.paraview.org/download/)


License
-------
GNU Lesser General Public License (LGPL). Please see the file LICENSE for details.
