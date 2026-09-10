#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

source ./run_local/env.sh

# mpiexecjl -n 8 julia --project=. test/wsi2dtest.jl case_1 &> output_wsi_2d_case_1.txt
# mpiexecjl -n 8 julia --project=. test/wsi2dtest.jl length_sweep &> output_wsi_2d_length_sweep.txt
# mpiexecjl -n 8 julia --project=. test/wsi2dtest.jl density_sweep &> output_wsi_2d_density_sweep.txt
mpiexecjl -n 8 julia --project=. test/wsi2dtest.jl all &> output_wsi_2d_all.txt

echo "WSI 2D test completed."
