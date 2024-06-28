Musubi
======

Octree based Lattice-Boltzmann solver with support for multiple species.
The actual sources are found in the musubi-source repository, which is
included here in the `mus` subdirectory.
Additionally all the other parts required for compilation of the solver
are included in this repository.

Use `git clone --recurse-submodules` when cloning this repository to fetch the
gathered subdirectories from the various repositories.

Prerequisite for building the solver is an installed Python, Fortran compiler
and MPI library. For compilation you need to point `FC` to the appropiate MPI
compiler wrapper. (Usually `export FC=mpif90`).

The solver can then be built with

```
bin/waf configure build
```

To install it, run:

```
bin/waf install
```

Run `bin/waf --help` to see all options.

Documentation
-------------

See the [documentation](https://geb.inf.tu-dresden.de/doxy/musubi/index.html)
for more details.

Developing Musubi
-----------------

The actual sources of musubi are found in the mus subdirectory, which is a git
submodule and thus has a repository (musubi-source) on its own.
To ease the work with this setup there is a `request` script that is meant to
take care of dealing with the tight coupling between this repository and
musubi-source.
