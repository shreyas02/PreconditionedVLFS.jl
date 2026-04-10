# PreconditionedVLFS.jl

`PreconditionedVLFS.jl` is a local Julia project for VLFS preconditioning experiments. This repository is used as a project workspace rather than a registered Julia package, with shell scripts in `run_local/` and `run_slurm/` driving setup, sysimage generation, and the main cases.

## Prerequisites

- Julia `1.10.x`
- System MPI available to Julia and Trilinos
- A Trilinos installation path exported through `run_local/env.sh` or `run_slurm/env.sh`

Trilinos must be built with the same MPI build used by the system binary configured through `MPIPreferences.jl`.

## Setup

### 1. Create the env.sh file

From the project root:

```bash
# For Local
cp run_local/env.example.sh run_local/env.sh
# For SLURM
cp run_slurm/env.example.sh run_slurm/env.sh
```

Edit `run_local/env.sh` or `run_slurm/env.sh` and set `TRILINOS_ROOT` for your machine. Add any required module loads or environment exports there as well.

### 2. Build the cpp source files related to the Gridap Trilinos interface

From the project root:

```bash
# For Local
bash run_local/buildcpp.sh
# For SLURM
sbatch slurm_jobs/<job-script-to-run-buildcpp.sh>.sh
```

### 3. Prepare the julia env

Then prepare the Julia environment and MPI preferences:

```bash
# For Local
bash run_local/warmup.sh
# For SLURM
bash run_slurm/warmup.sh
```

`run_local/warmup.sh` and `run_slurm/warmup.sh` install the project dependencies, select the system MPI binary through `MPIPreferences`, and run standard Julia precompilation.

### 4. Build Sysimage

To build the MPI-aware sysimage from the warmup traces, run:

```bash
# For Local
bash run_local/build_sysimage.sh
# For SLURM
sbatch slurm_jobs/<job-script-to-run-build_sysimage.sh>.sh
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

## Write your own jobs script for the SLRUM workflow

Use `example_job.sh` as a template and write your job script in the `slurm_jobs` dir. 

## Run Cases

Run the provided shell entrypoints from the project root:

```bash
# For Local
bash run_local/toyrichardson.sh
bash run_local/wsi2d.sh
bash run_local/wsi3d.sh
bash run_local/scaling_2d_strong.sh
bash run_local/scaling_2d_weak.sh
bash run_local/scaling_3d_strong.sh
bash run_local/scaling_3d_weak.sh
```

```bash
# For SLURM
sbatch slurm_jobs/<job-script>.sh
```

Adjust the number of MPI processes by editing the corresponding shell script before running it.

The run scripts use the sysimage at `compile/PreconditionedVLFS.so`.
