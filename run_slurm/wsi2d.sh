#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source ./run_slurm/env.sh

srun --mpi=pmix --ntasks="${SLURM_NTASKS}" julia --project=. \
    -J compile/PreconditionedVLFS.so \
    test/wsi2dtest.jl \
    > "${ROOT_DIR}/slurm_jobs/wsi_2d.log" 2>&1

echo "WSI 2D test completed."
