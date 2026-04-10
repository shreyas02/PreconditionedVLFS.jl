#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source ./run_slurm/env.sh

for n in 8 6 4 2; do
    echo "Running ${n}-rank case"
    srun --mpi=pmix --ntasks="${n}" julia --project=. \
        -J compile/PreconditionedVLFS.so \
        test/periodic3d_test_strong.jl \
        > "${ROOT_DIR}/slurm_jobs/periodic3d_test_strong_${n}cores.log" 2>&1
done

echo "Strong scaling test for 3D completed."
