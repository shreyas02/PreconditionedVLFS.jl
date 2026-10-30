#!/bin/bash
# Run this script from anywhere inside the repository

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source ./run_slurm/env.sh

julia --project=. <<'JULIA'
using Pkg

patched_libs = [
    PackageSpec(url="https://github.com/shreyas02/GridapDistributed.jl", rev="preconditioner"),
    PackageSpec(url="https://github.com/shreyas02/GridapTrilinos.jl"),
    PackageSpec(url="https://github.com/Kyjor/FixedPointNumbers.jl",rev="167969b",),
]

Pkg.add(patched_libs)
Pkg.resolve()
Pkg.instantiate()

using MPIPreferences
MPIPreferences.use_system_binary()
JULIA

echo "Warmup complete: project instantiated and configured for system MPI."
