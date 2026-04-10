#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

RUN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd -- "$RUN_DIR/.." && pwd)
BUILD_DIR="$ROOT_DIR/src/cppcode"

source "$RUN_DIR/env.sh" # Place the env.sh file in the run directory

echo "The build directory is: $BUILD_DIR"
cd "$BUILD_DIR" || { echo "Failed to enter build directory"; exit 1; }
mkdir -p build
cd build
rm -rf * # Clean the build directory before configuring

# 2. Source the configure script
echo "Sourcing the configure script..."
source ../configure.sh

# 3. Detect available CPU cores and compile
cores=$(nproc)
echo "===================================================="
echo "Configuration finished. Starting compilation..."
echo "Running: make -j$cores"
echo "===================================================="

make -j"$cores"

cd "$ROOT_DIR"
echo "Configuration done!"