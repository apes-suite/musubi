Musubi
======

Octree based Lattice-Boltzmann solver with support for multi-physics.
The actual sources are found in the
[musubi-source repository](https://github.com/apes-suite/musubi-source), which is
included here in the `mus` subdirectory.
Additionally all the other parts required for compilation of the solver
are included in this repository.

Use `git clone --recurse-submodules` when cloning this repository to fetch the
gathered submodules from the various repositories.

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

We support two methods to conveniently install a development environment
with the APES tools:

* via the Spack manager with the packages in [apes-spack](https://github.com/apes-suite/apes-spack)
* via a Python virtual environment as provided in [apes-pyenv](https://github.com/apes-suite/apes-pyenv)

Please see the respective READMEs on instructions on how to
use one of these methods.

Documentation
-------------

See the [documentation](https://apes-suite.github.io/musubi/index.html)
for more details.

Developing Musubi
-----------------

The actual sources of musubi are found in the mus subdirectory, which is a git
submodule and, thus, has a repository
([musubi-source](https://github.com/apes-suite/musubi-source)) on its own.
To ease the work with this setup there is a `request` script in the `bin` subdirectory
that is meant to take care of dealing with the tight coupling between this repository
and musubi-source.

Testing
-------

The unit tests are automatically run during the compilation step, unless you
tell waf to only build musubi itself or explicitly switch off test runs.
The system tests are found in the `mus/examples` subdirectory and they can
be run with the help of pysys-test.
To obtain pysys-test and have it available, you can use the
[apes-pyenv](https://github.com/apes-suite/apes-pyenv).

While in the apes-pyenv environment install seeder and musubi as suggested
into the `$VIRTUAL_ENV` prefix and then run `pysys.py run` in the `mus/examples`
directory (or any of its subdirectories if you want to run only the
tests within that directory).
You can also list the available tests with `pysys.py ls` and see how to
run selected tests in its help (`pysys.py run --help`).

If you don't want to install Musubi into your `PATH`, you can also tell
pysys-test where to find it by setting the environment variable `APES_MUSUBI`
accordingly.
