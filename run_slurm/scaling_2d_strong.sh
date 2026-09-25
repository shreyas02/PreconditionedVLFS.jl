#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source ./run_slurm/env.sh

for n in 64 128 256 512; do
    nodes=$((n / 64))  # Assuming 64 cores per node
    echo "Running ${n}-rank case in ${nodes} nodes"
    srun --mpi=pmix --nodes="${nodes}" --ntasks="${n}" julia --project=. \
        test/periodic2dtest.jl strong_scaling \
        > "${ROOT_DIR}/slurm_jobs/periodic2d_test_strong_${n}cores.log" 2>&1
done

echo "Strong scaling test for 2D completed."
