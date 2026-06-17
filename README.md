# lava-flow
A finite element simulation code on adaptive quadtrees for lava flow dynamics, built on top of BIM++ (landslides branch). Includes a travelling vortex benchmark for validation.

## Dependencies

All dependencies are downloaded and built automatically by CMake:

| Library | Version / Source | Purpose |
|---|---|---|
| [LIS](https://www.ssisc.org/lis/) | 2.1.10 | Iterative linear solvers |
| [GNU Octave](https://www.gnu.org/software/octave/) | 6.2.0 | Scripting and post-processing |
| [octave_file_io](https://github.com/carlodefalco/octave_file_io) | git HEAD | Octave-based file I/O bridge |
| [BIM++](https://github.com/carlodefalco/bimpp) | `landslides` branch | Core finite-element/volume library |

## System Requirements

The following must be available on your system before building:

- CMake ≥ 3.17 (loaded: cmake/3.30.5)
- MPI compilers: `mpicc`, `mpicxx`, `mpif90` (loaded: openmpi/4.1.7)
- OpenMP-capable compilers — GCC recommended (loaded: gcc/12.2.0)
- Fortran compiler: `gfortran` (included with GCC)
- Standard build tools: `make`, `autoconf`, `automake`, `libtool`
- MUMPS and its dependencies (loaded: mumps):
  - `dmumps`, `mumps_common`, `pord`
  - `metis/5.1.0`, `parmetis/4.0.3`
  - `scotch/7.0.4`, `scotcherr`, `esmumps`
  - `scalapack` (loaded: netlib-scalapack/2.2.0)
  - `openblas/0.3.28`
- p4est/2.8 (adaptive mesh refinement)
- Lua/5.4.6 (if used for configuration/scripting)
- `ncurses/6.5`, `readline/8.2`, `pcre/8.45` (terminal/scripting dependencies)

## Building

```bash
# 1. Clone the repository
git clone <repo-url>
cd lava-flow

# 2. Configure
cmake -S . -B build

# 3. Build dependencies (first time only, this will take a while)
cd build
make deps -j$(nproc)

# 4. Build the lava-flow executable
make -j$(nproc)
```

After the first `make deps`, subsequent builds only require `make -j$(nproc)`.

## Running

The simulation is configured via JSON input files. Example input files are provided in the repository:

```bash
# Well-balancing test
mpirun -np NPROCS lava-flow ../glisX_input-lava-wb.json

# Lava pouring scenario
mpirun -np NPROCS lava-flow ../glisX_input-lava-pouring.json

# Travelling vortex benchmark
mpirun -np NPROCS lava-flow ../glisX_input-lava-travelling-vortex.json
```

Once runned the code, call from the build directory
```bash
../convertTovtk.sh
```

## Citation

If you use this code in your research, please cite:

```bibtex
@article{gatti2025second,
  title={Second-order Optimally Stable IMEX (pseudo-) staggered Galerkin discretization: application to lava flow modeling},
  author={Gatti, Federico and Orlando, Giuseppe},
  journal={arXiv preprint arXiv:2509.09460},
  year={2025}
}
```
