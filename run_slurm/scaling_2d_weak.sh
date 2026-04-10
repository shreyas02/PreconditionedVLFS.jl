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
        test/periodic2d_test_weak.jl \
        > "${ROOT_DIR}/slurm_jobs/periodic2d_test_weak_${n}cores.log" 2>&1
done

echo "Weak scaling test for 2D completed."
