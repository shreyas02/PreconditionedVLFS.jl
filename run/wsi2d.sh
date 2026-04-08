source ./run/env.sh

mpiexecjl -n 8 julia --project=. -J compile/PreconditionedVLFS.so test/wsi2dtest.jl &> output_wsi_2d.txt

echo "WSI 2D test completed."