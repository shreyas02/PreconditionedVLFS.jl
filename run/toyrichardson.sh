source ./run/env.sh

julia --project=. -J compile/PreconditionedVLFS.so test/toyrichardsontest.jl &> output_toy_richardson.txt

echo "Toy Richardson test completed."