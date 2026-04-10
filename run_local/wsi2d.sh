#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

source ./run_local/env.sh

mpiexecjl -n 8 julia --project=. -J compile/PreconditionedVLFS.so test/wsi2dtest.jl &> output_wsi_2d.txt

echo "WSI 2D test completed."