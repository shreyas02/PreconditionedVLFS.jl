#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source ./run_slurm/env.sh

TRACE_DIR="${ROOT_DIR}/compile/traces"
SERIAL_TRACE="${TRACE_DIR}/serial_trace.jl"
MERGED_TRACE="${ROOT_DIR}/compile/PreconditionedVLFS_precompile_statements.jl"

mkdir -p "${TRACE_DIR}"
rm -f "${SERIAL_TRACE}" "${MERGED_TRACE}" "${TRACE_DIR}"/mpi_rank_*.jl

julia --project=. -e '
    using MPIPreferences
    MPIPreferences.use_system_binary()
'

julia --project=. -e '
    using MPIPreferences
    println("MPIPreferences binary=$(MPIPreferences.binary), abi=$(MPIPreferences.abi)")
    if string(MPIPreferences.binary) != "system"
        error("Expected system MPI before tracing and building the sysimage.")
    end
'

julia --project=. -e '
    using Pkg
    Pkg.instantiate()
    Pkg.precompile()
'

echo "Tracing serial warmup into ${SERIAL_TRACE}"
julia --project=. --trace-compile="${SERIAL_TRACE}" compile/trace_serial_warmup.jl

echo "Tracing MPI warmups with 2 ranks"
srun --mpi=pmix -n 2 bash -lc '
    set -euo pipefail
    rank="${SLURM_PROCID}"
    trace_file="compile/traces/mpi_rank_${rank}.jl"
    julia --project=. --trace-compile="${trace_file}" compile/trace_mpi_warmup.jl
'

awk '/^precompile/ { print }' "${SERIAL_TRACE}" "${TRACE_DIR}"/mpi_rank_*.jl | sort -u > "${MERGED_TRACE}"

echo "Building sysimage from ${MERGED_TRACE}"
julia --project=. compile/build_sysimage.jl

echo "Verifying sysimage can load MPI"

srun --mpi=pmix -n 1 julia --project=. -J compile/PreconditionedVLFS.so -e '
    using MPI
    MPI.Init()
    println("sysimage MPI ok on rank ", MPI.Comm_rank(MPI.COMM_WORLD))
    MPI.Finalize()
'

echo "Sysimage build complete: compile/PreconditionedVLFS.so"
