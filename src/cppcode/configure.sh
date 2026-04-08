#!/bin/bash

# Adding the path to Trilinos 
if [ -z "$TRILINOS_ROOT" ]; then
    echo "Error: TRILINOS_ROOT environment variable is not set."
    echo "Please set it before running this script. For example:"
    echo "export TRILINOS_ROOT=/path/to/TrilinosInstall"
    exit 1
fi

echo "Using Trilinos installation at: $TRILINOS_ROOT"

export LD_LIBRARY_PATH="${TRILINOS_ROOT}/lib:$LD_LIBRARY_PATH"
export CMAKE_PREFIX_PATH="${TRILINOS_ROOT}:$CMAKE_PREFIX_PATH"

SOURCE_DIR=../cppsrc
BUILD_DIR=`pwd`

cmake \
    ${SOURCE_DIR}