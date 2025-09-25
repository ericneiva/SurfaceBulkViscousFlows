# SurfaceBulkViscousFlows

Repository that holds the Julia software with code demonstrators of the numerical examples in the paper _"Unfitted finite element modelling of surface-bulk viscous flows in animal cells"_. https://arxiv.org/abs/2505.05723

## List of examples

Referred to the paper sections

- _Section 4.1 Verification examples_
  - [Verification in fixed sphere](examples/Verification/VerificationInFixedSphere.jl)
  - [Verification with dynamic surface](examples/Verification/VerificationWithDynamicSurface.jl)
- _Section 4.2 Self-organised shape emergence_
  - [Shape emergence](examples/ShapeEmergence/ShapeEmergence.jl)
- _Section 4.3 Relaxation dynamics_
  - [Pear](examples/RelaxationDynamics/Pear.jl)
  - [Popcorn](examples/RelaxationDynamics/Popcorn.jl)
  - [Torus](examples/RelaxationDynamics/Torus.jl)
- _Section 4.4 Cell cleavage_
  - [2D axisymmetric division](examples/CellCleavage/2DAxisymmetricCleavage.jl)
  - [3D division](examples/CellCleavage/3DCleavage.jl)

## Instructions to run the examples

### Prerequisites

  - A Julia installation ([link to instructions](https://julialang.org/install/))
  - Install PETSc version 3.15.2 with mumps and let Julia find the library:

**Linux instructions:** On your folder of choice;

```shell
  wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.15.2.tar.gz
  tar -xvzf petsc-3.15.2.tar.gz
  cd petsc-3.15.2.tar.gz
  ./configure -download-mumps -download-scalapack -download-parmetis -download-metis --download-bison -download-ptscotch --download-make --with-debugging=0 --download-mpich --download-fblaslapack=1
  make PETSC_DIR=/your-folder/petsc-3.15.2 PETSC_ARCH=arch-linux-c-opt check
```

Note that for PETSc you might need to install missing required libraries.

Add `export JULIA_PETSC_LIBRARY="/your-folder/petsc-3.15.2/arch-linux-c-opt/lib/libpetsc.so"` to your `.bashrc` file.

### Running an example

  - Clone the repository at a path of your choice
  - Open a Julia interactive session on the repository directory

```
  git clone git@github.com:ericneiva/SurfaceBulkViscousFlows.git
  cd SurfaceBulkViscousFlows
  julia --project
```

  - Before _running an example for the first time_, set up the project dependencies

```julia
  julia> # press ] to enter Pkg mode...
  (SurfaceBulkViscousFlows) pkg> instantiate # Check that Julia points to your local PETSc installation
```

  - Run the example of your choice

```julia
  julia> include("examples/Verification/VerificationInFixedSphere.jl")
```

  - Output files are generated in the folder for inspection

### How to cite this code

In order to give credit to the contributors of this software, we simply ask you to cite the references below in any publication in which you have made use of the repository.



### Contact

Please, contact [Eric Neiva](mailto:eric.neiva@college-de-france.fr) or [Hervé Turlier](mailto:herve.turlier@college-de-france.fr) if you have any questions.
