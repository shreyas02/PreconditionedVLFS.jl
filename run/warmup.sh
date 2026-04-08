#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

source ./run/env.sh

set -euo pipefail

julia --project=. -e '
    using DrWatson
    @quickactivate "PreconditionedVLFS"

    using MPIPreferences
    MPIPreferences.use_system_binary()
'

julia --project=. <<'JULIA'
    using Pkg
    using DrWatson
    @quickactivate "PreconditionedVLFS"

    patched_libs = [PackageSpec(url="https://github.com/shreyas02/Gridap.jl", rev="issue-1191"),
                    PackageSpec(url="https://github.com/shreyas02/GridapDistributed.jl", rev="preconditioner"),
                    PackageSpec(url="https://github.com/shreyas02/GridapSolvers.jl.git", rev="richardson_bugfix_v0.6.1"),
                    PackageSpec(name="PackageCompiler")]

    Pkg.add(patched_libs)
    Pkg.instantiate()
    Pkg.precompile()
JULIA

echo "Warmup complete: project instantiated, precompiled, and configured for system MPI."
echo "Build the MPI-aware sysimage separately with: bash run/build_sysimage.sh"