#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

source ./run_local/env.sh

mpiexecjl -n 8 julia --project=. test/periodic2dtest.jl weak_scaling &> output_8cores.txt
mpiexecjl -n 6 julia --project=. test/periodic2dtest.jl weak_scaling &> output_6cores.txt
mpiexecjl -n 4 julia --project=. test/periodic2dtest.jl weak_scaling &> output_4cores.txt
mpiexecjl -n 2 julia --project=. test/periodic2dtest.jl weak_scaling &> output_2cores.txt

echo "Weak scaling test for 2D completed."
