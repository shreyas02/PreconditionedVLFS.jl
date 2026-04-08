# PreconditionedVLFS.jl

`PreconditionedVLFS.jl` is a local Julia project for VLFS preconditioning experiments. This repository is used as a project workspace rather than a registered Julia package, with shell scripts in `run/` driving setup, sysimage generation, and the main cases.

## Prerequisites

- Julia `1.10.x`
- System OpenMPI available to Julia and Trilinos
- A Trilinos installation path exported through `run/env.sh`

## Setup

From the project root:

```bash
cp run/env.example.sh run/env.sh
```

Edit `run/env.sh` and set `TRILINOS_ROOT` for your machine. Add any required module loads or environment exports there as well.

Then prepare the Julia environment and MPI preferences:

```bash
bash run/warmup.sh
```

`run/warmup.sh` installs the project dependencies, selects the system MPI binary through `MPIPreferences`, and runs standard Julia precompilation.

## Build the cpp source files related to the Gridap Trilinos interface

From the project root:

```bash
bash run/buildcpp.sh
```

## Build Sysimage

To build the MPI-aware sysimage from the warmup traces, run:

```bash
bash run/build_sysimage.sh
```

This script:

- verifies that Julia is configured to use system OpenMPI
- traces the serial and MPI warmup cases
- merges the generated precompile statements
- builds a sysimage for `PreconditionedVLFS`

The output sysimage is written to:

```text
compile/PreconditionedVLFS.so
```

## Trilinos

Expected Trilinos setup:

- Trilinos `16.10`
- Required libraries: `Amesos2`, `Teuchos`, `Tpetra`, `Belos`, `ShyLU_DDFROSch`, `Xpetra`, `MueLu`
- SuiteSparse `7.8.2`

For the parallel cases, solver settings can be adjusted through the XML files in `data/`.

## Run Cases

Run the provided shell entrypoints from the project root:

```bash
bash run/toyrichardson.sh
bash run/wsi2d.sh
bash run/wsi3d.sh
bash run/scaling_2d_strong.sh
bash run/scaling_2d_weak.sh
bash run/scaling_3d_strong.sh
bash run/scaling_3d_weak.sh
```

Adjust the number of MPI processes by editing the corresponding shell script before running it.

The run scripts use the sysimage at `compile/PreconditionedVLFS.so`.