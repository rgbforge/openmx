#!/bin/bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
WORK_DIR="${SCRIPT_DIR}/work"

module load intel/2025.0.1/tbb/latest
module load intel/2025.0.1/umf/latest
module load intel/2025.0.1/compiler-rt/latest
module load intel/2025.0.1/compiler/latest
module load intel/2025.0.1/mpi/latest
module load intel/2025.0.1/mkl/latest

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake ..
make
make install
