#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

source ./run_local/env.sh

mpiexecjl -n 8 julia --project=. -J compile/PreconditionedVLFS.so test/wsi3dtest.jl &> output_wsi_3d.txt

echo "WSI 3D test completed."