source ./run/env.sh

source ./run/env.sh

mpiexecjl -n 8 julia --project=. -J compile/PreconditionedVLFS.so test/periodic2d_test_strong.jl &> output_8cores.txt
mpiexecjl -n 6 julia --project=. -J compile/PreconditionedVLFS.so test/periodic2d_test_strong.jl &> output_6cores.txt
mpiexecjl -n 4 julia --project=. -J compile/PreconditionedVLFS.so test/periodic2d_test_strong.jl &> output_4cores.txt
mpiexecjl -n 2 julia --project=. -J compile/PreconditionedVLFS.so test/periodic2d_test_strong.jl &> output_2cores.txt

echo "Strong scaling test for 2D completed."