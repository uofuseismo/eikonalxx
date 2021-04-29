# Overview

This is a library for solving the eikonal equation in Cartesian coordinates using the fast sweeping method.  The implementation relies on SYCL thereby allowing the same code to be run on a CPU or GPU.

# Installation

## Prerequisites

The library requires

   1. CMake 3.10 or higher.
   1. A C++20 compliant compiler.
   2. A SYCL v1.2 compliant compiler - e.g., dpcpp or Intel.

Optionally, to build the Python bindings the following are needed

   1. Python3  
   2. pybind11

## Configuration

An example configuration script utilizing the Intel compiler is given

    #!/usr/bin/bash
    BUILD_DIR=intel_build
    if [ -d ${BUILD_DIR} ]; then
       rm -rf ${BUILD_DIR}
    fi
    mkdir ${BUILD_DIR}
    cd ${BUILD_DIR}
    cmake .. \
    -DCMAKE_CXX_COMPILER=/opt/intel/oneapi/compiler/2021.1.1/linux/bin/icpx \
    -DCMAKE_CXX_FLAGS="-Wall -tbb -fsycl -fsycl-unnamed-lambda" \
    -DWRAP_PYTHON=TRUE \
    -DPYTHON_EXECUTABLE=/home/bbaker/anaconda3/bin/python \
    -DPYTHON_INCLUDE_DIRS=/home/bbaker/anaconda3/include/python3.8 \
    -DPYTHON_LIBRARIES=/home/bbaker/anaconda3/lib

After configuration descend into the build directory and type

    make
    make test

