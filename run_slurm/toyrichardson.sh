#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source ./run_slurm/env.sh

srun --ntasks=1 julia --project=. -J compile/PreconditionedVLFS.so test/toyrichardsontest.jl \
    > "${ROOT_DIR}/slurm_jobs/toy_richardson.log" 2>&1

echo "Toy Richardson test completed."
