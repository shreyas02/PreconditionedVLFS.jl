source ./run/env.sh

source ./run/env.sh

mpiexecjl -n 8 julia --project=. -J compile/PreconditionedVLFS.so test/periodic3d_test_weak.jl &> output_8cores.txt
mpiexecjl -n 6 julia --project=. -J compile/PreconditionedVLFS.so test/periodic3d_test_weak.jl &> output_6cores.txt
mpiexecjl -n 4 julia --project=. -J compile/PreconditionedVLFS.so test/periodic3d_test_weak.jl &> output_4cores.txt
mpiexecjl -n 2 julia --project=. -J compile/PreconditionedVLFS.so test/periodic3d_test_weak.jl &> output_2cores.txt

echo "Weak scaling test for 3D completed."