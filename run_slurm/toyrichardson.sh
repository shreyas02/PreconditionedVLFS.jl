#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source ./run_slurm/env.sh

srun --ntasks=1 julia --project=. test/toyrichardsontest.jl all \
    > "${ROOT_DIR}/slurm_jobs/toy_richardson_all.log" 2>&1

echo "Toy Richardson test completed."
