#!/bin/bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build_openmx_intel"
mkdir -p "$BUILD_DIR"

module load intel/2025.0.1/tbb/latest
module load intel/2025.0.1/umf/latest
module load intel/2025.0.1/compiler-rt/latest
module load intel/2025.0.1/compiler/latest
module load intel/2025.0.1/mpi/latest
module load intel/2025.0.1/mkl/latest

export CC="mpiicx -O3 -xHOST -fiopenmp -fcommon -Wno-error=implicit-function-declaration -I${MKLROOT}/include -I${MKLROOT}/include/fftw"
export FC="mpiifx -O3 -xHOST -fiopenmp"
export LIB="-L${MKLROOT}/lib -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl -lifcore"
NUM_JOBS=${NUM_JOBS:-$(nproc || echo 1)}

OPENMX_VERSION="3.9"
OPENMX_PATCH="3.9.9"
OPENMX_SOURCE_DIR="$BUILD_DIR/openmx${OPENMX_VERSION}"
OPENMX_EXEC="$OPENMX_SOURCE_DIR/work/openmx"

build_openmx() {
    if [ ! -f "$OPENMX_EXEC" ]; then
        local OPENMX_TARBALL="$BUILD_DIR/openmx${OPENMX_VERSION}.tar.gz"
        local PATCH_TARBALL="$BUILD_DIR/patch${OPENMX_PATCH}.tar.gz"

        if [ ! -f "$OPENMX_TARBALL" ]; then
            echo "downloading openmx $OPENMX_VERSION..."
            wget -P "$BUILD_DIR" "https://www.openmx-square.org/openmx${OPENMX_VERSION}.tar.gz"
        fi

        if [ ! -f "$PATCH_TARBALL" ]; then
            echo "downloading patch $OPENMX_PATCH ..."
            wget -P "$BUILD_DIR" "https://www.openmx-square.org/bugfixed/21Oct17/patch${OPENMX_PATCH}.tar.gz"
        fi

        if [ -d "$OPENMX_SOURCE_DIR" ]; then
            echo "cleaning for rebuild..."
            rm -rf "$OPENMX_SOURCE_DIR"
        fi

        echo "extracting..."
        tar -C "$BUILD_DIR" -xzf "$OPENMX_TARBALL"

        echo "patching..."
        cd "$OPENMX_SOURCE_DIR/source"
        tar -xzf "$PATCH_TARBALL"
        mv kpoint.in ../work/

        echo "configuring..."
        sed -i "s#^CC\s*=.*#CC = ${CC}#" makefile
        sed -i "s#^FC\s*=.*#FC = ${FC}#" makefile
        sed -i "s#^LIB\s*=.*#LIB = ${LIB}#" makefile


        make clean

        echo "building elpa..."
        make mod_precision.o elpa_utilities.o elpa1_compute_real.o elpa1_compute_complex.o \
             aligned_mem.o elpa2_determine_workload.o mod_redist_band_real.o \
             mod_redist_band_complex.o mod_pack_unpack_cpu_real.o mod_pack_unpack_cpu_complex.o \
             real.o complex.o mod_single_hh_trafo_real.o mod_compute_hh_trafo_real.o \
             mod_compute_hh_trafo_complex.o elpa2_compute_real.o elpa2_compute_complex.o \
             elpa_solve_evp_real_2stage_double_impl.o elpa_solve_evp_complex_2stage_double_impl.o

        echo "building openmx..."
        make -j"$NUM_JOBS" all
        echo "installing openmx..."
        make install

        echo "build complete"
    else
        echo "openmx found at ${OPENMX_EXEC}, skipping build..."
    fi
}
clean_build() {
    rm -rf "$BUILD_DIR"
}

usage() {
    echo "usage: $0 -b [openmx|all|clean]"
    exit 1
}

if [ "$#" -eq 0 ]; then
    usage
fi

BUILD_TARGET=""
while getopts ":b:" opt; do
  case ${opt} in
    b) BUILD_TARGET="${OPTARG}" ;;
    \?|:) usage ;;
  esac
done

case "${BUILD_TARGET}" in
    openmx|all)
        build_openmx
        ;;
    clean)
        clean_build
        ;;
    *)
        usage
        ;;
esac
