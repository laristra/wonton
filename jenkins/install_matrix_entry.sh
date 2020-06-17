#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     $WORKSPACE/jenkins/build_matrix_entry.sh <compiler> <build_type>
#
# The exit code determines if the test succeeded or failed.
# Note that the environment variable WORKSPACE must be set (Jenkins
# will do this automatically).

# Exit on error
set -e
# Echo each command
set -x

compiler=$1
build_type=$2

echo "inside install_matrix entry"

# set modules and install paths

jali_version=1.1.4
lapack_version=3.8.0

export NGC=/usr/local/codes/ngc
ngc_include_dir=$NGC/private/include

wonton_install_dir=$NGC/private/wonton
wonton_version=dev

thrust_dir=${ngc_include_dir}

# compiler-specific settings
if [[ $compiler =~ "intel" ]]; then

    compiler_version=18.0.1
    cxxmodule=intel/${compiler_version}
    compiler_suffix="-intel-${compiler_version}"

    openmpi_version=2.1.2
    mpi_module=openmpi/2.1.2
    mpi_suffix="-openmpi-${openmpi_version}"
    
elif [[ $compiler =~ "gcc" ]]; then

    openmpi_version=2.1.2
    if [[ $compiler == "gcc6" ]]; then
	compiler_version=6.4.0
    elif [[ $compiler == "gcc7" ]]; then
	compiler_version=7.3.0
    elif [[ $compiler == "gcc8" ]]; then
	compiler_ersion=8.2.0
	openmpi_version=3.1.3
    fi
    
    cxxmodule=gcc/${compiler_version}
    compiler_suffix="-gcc-${compiler_version}"
    
    mpi_module=openmpi/${openmpi_version}
    mpi_suffix="-openmpi-${openmpi_version}"

fi

jali_install_dir=$NGC/private/jali/${jali_version}${compiler_suffix}${mpi_suffix}
jali_flags="-D WONTON_ENABLE_Jali:BOOL=True -D Jali_ROOT:PATH=$jali_install_dir"

lapacke_dir=$NGC/private/lapack/${lapack_version}-patched${compiler_suffix}
lapacke_flags="-D WONTON_ENABLE_LAPACKE:BOOL=True -D LAPACKE_ROOT:PATH=$lapacke_dir"

if [[ $compiler == "gcc6" ]]; then
    flecsi_flags="-D WONTON_ENABLE_FleCSI:BOOL=True"
fi


mpi_flags="-D WONTON_ENABLE_MPI=True"
if [[ $build_type == "serial" ]]; then
    mpi_flags="-D WONTON_ENABLE_MPI=False"
    mpi_suffix=
    jali_flags=
    flecsi_flags=
fi


thrust_flags=
thrust_suffix=
thrust_version=1.8.3
thrust_install_dir=$NGC/private/thrust/${thrust_version}
if [[ $build_type == "thrust" ]]; then
    thrust_flags="-D WONTON_ENABLE_THRUST=True -D THRUST_ROOT=${thrust_install_dir}"
    thrust_suffix="-thrust"
fi

kokkos_flags=
kokkos_suffix=
kokkos_version=3.1.01
kokkos_install_dir=$NGC/private/kokkos/${kokkos_version}${compiler_suffix}
if [[ $build_type == "kokkos" ]]; then
    kokkos_flags="-D WONTON_ENABLE_Kokkos=True -D Kokkos_ROOT=${kokkos_install_dir}"
    kokkos_suffix="-kokkos"
fi

wonton_install_dir=$wonton_install_dir/${wonton_version}${compiler_suffix}${mpi_suffix}${thrust_suffix}${kokkos_suffix}

export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load $cxxmodule
if [[ -n "$mpi_flags" ]]; then
    module load ${mpi_module}
fi
module load cmake/3.14.0 # 3.13 or higher is required

module load git
branch=`git rev-parse --abbrev-ref HEAD`

if [[ branch == "master" && build_type == "kokkos" ]]; then
    exit
fi

if [[ branch == "kokkos" && build_type != "kokkos" ]]; then
    exit
fi

cmake_build_type=Release

echo "JENKINS WORKSPACE = $WORKSPACE"
cd $WORKSPACE

rm -fr build
mkdir -p build
cd build

cmake \
    -D CMAKE_BUILD_TYPE=$cmake_build_type \
    -D CMAKE_CXX_FLAGS="-Wall -Werror" \
    -D CMAKE_INSTALL_PREFIX=$wonton_install_dir \
    -D CMAKE_PREFIX_PATH=$ngc_include_dir \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_JENKINS_OUTPUT=True \
    $mpi_flags \
    $jali_flags \
    $lapacke_flags \
    $thrust_flags \
    $kokkos_flags \
    .. &&
make -j2 &&
ctest --output-on-failure &&
make install
