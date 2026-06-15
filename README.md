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

- CMake ≥ 3.17
- MPI compilers: `mpicc`, `mpicxx`, `mpif90`
- OpenMP-capable compilers (GCC recommended)
- MUMPS and its dependencies: `dmumps`, `mumps_common`, `metis`, `pord`, `esmumps`, `scotch`, `scotcherr`, `scalapack`, `openblas`
- Standard build tools: `make`, `autoconf`, `automake`, `libtool`
- Fortran compiler (`gfortran`)

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
# Lava flow simulation
./lava-flow ../glisX_input_lava.json

# Lava pouring scenario
./lava-flow ../glisX_input-lava-pouring.json

# Travelling vortex benchmark
./lava-flow ../glisX_input-lava-travelling-vortex.json
```

Once runned the code, call from the build directory
```bash
../convertTovtk.sh
```

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{lava-flow,
  author = {TODO},
  title  = {lava-flow},
  year   = {TODO},
  url    = {TODO}
}
```
