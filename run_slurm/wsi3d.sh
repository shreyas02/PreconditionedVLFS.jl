#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

set -euo pipefail

: "${SLURM_JOB_ID:?This script must be run inside a Slurm job allocation}"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source ./run_slurm/env.sh

echo "Running WSI 3D test (8 ranks)"

srun --mpi=pmix --ntasks=8 julia --project=. \
    -J compile/PreconditionedVLFS.so \
    test/wsi3dtest.jl \
    > "${ROOT_DIR}/slrum_jobs/wsi_3d.log" 2>&1

echo "WSI 3D test completed."
