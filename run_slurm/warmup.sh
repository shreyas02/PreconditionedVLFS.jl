#!/bin/bash
# Run this script from anywhere inside the repository

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source ./run_slurm/env.sh

julia --project=. -e '
    using Pkg
    Pkg.resolve()
    using MPIPreferences
    MPIPreferences.use_system_binary()
    Pkg.precompile()
'

julia --project=. <<'JULIA'
using Pkg

patched_libs = [
    PackageSpec(url="https://github.com/shreyas02/Gridap.jl", rev="issue-1191"),
    PackageSpec(url="https://github.com/shreyas02/GridapDistributed.jl", rev="preconditioner"),
    PackageSpec(url="https://github.com/shreyas02/GridapSolvers.jl.git", rev="richardson_bugfix_v0.6.1"),
    PackageSpec(url="https://github.com/Kyjor/FixedPointNumbers.jl",rev="167969b",),
    PackageSpec(name="PackageCompiler"),
]

Pkg.add(patched_libs)
Pkg.instantiate()
Pkg.precompile()
JULIA

echo "Warmup complete: project instantiated, precompiled, and configured for system MPI."