#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

source ./run_local/env.sh

# mpiexecjl -n 8 julia --project=. test/wsi3dtest.jl case_1 &> output_wsi_3d_case_1.txt
# mpiexecjl -n 10 julia --project=. test/wsi3dtest.jl case_2 &> output_wsi_3d_case_2.txt
mpiexecjl -n 10 julia --project=. test/wsi3dtest.jl all &> output_wsi_3d_all.txt

echo "WSI 3D test completed."
