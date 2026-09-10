#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

source ./run_local/env.sh

# julia --project=. test/toyrichardsontest.jl comparison &> output_toy_richardson_comparison.txt
julia --project=. test/toyrichardsontest.jl all &> output_toy_richardson_all.txt

echo "Toy Richardson test completed."
