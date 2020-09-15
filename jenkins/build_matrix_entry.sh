#!/usr/bin/env bash
# This script is executed on Jenkins using
#
# $WORKSPACE/jenkins/build_matrix_entry.sh build_type config_type compiler <VER>
#
# build_type  -  pr, nightly, install
# config_type -  base, debug, serial, readme, thrust, kokkos
# compiler    -  intel, gcc6, gcc7
#
# if VER is absent, the HEAD of the master branch will be taken
# (except for kokkos config_type which takes the HEAD of the kokkos
# branch)
#
# The exit code determines if the test succeeded or failed.
# Note that the environment variable WORKSPACE must be set (Jenkins
# will do this automatically).

# Exit on error
set -e
# Echo each command
set -x

build_type=$1
config_type=$2
compiler=$3
version=$4
if [[ $version == "" ]]; then
    version=dev
fi

# Don't build kokkos config from master branch and general configs
# from kokkos branch
if [[ $BRANCH_NAME == "master" && config_type == "kokkos" ]]; then
    exit
fi
if [[ $BRANCH_NAME == "kokkos" && config_type != "kokkos" ]]; then
    exit
fi


# versions of dependencies
jali_version=1.1.4
thrust_version=1.8.3
kokkos_version=3.1.01
lapack_version=3.8.0


echo "inside build_matrix with build_type=$build_type config_type=$config_type compiler=$compiler"

jali_version=1.1.4
thrust_version=1.8.3
kokkos_version=3.1.01
lapack_version=3.8.0

# special case for README builds
if [[ $build_type != "install" && config_type == "readme" ]]; then

    # Put a couple of settings in place to generate test output even if
    # the README doesn't ask for it.
    export CTEST_OUTPUT_ON_FAILURE=1
    CACHE_OPTIONS="-D ENABLE_JENKINS_OUTPUT=True"
    sed "s/^ *cmake/& $CACHE_OPTIONS/g" $WORKSPACE/README.md >$WORKSPACE/README.md.1
    python2 $WORKSPACE/jenkins/parseREADME.py \
	    $WORKSPACE/README.md.1 \
	    $WORKSPACE \
	    varan
    exit

fi

# set modules and install paths

export NGC=/usr/local/codes/ngc
ngc_include_dir=$NGC/private/include


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
    fi  
    cxxmodule=gcc/${compiler_version}
    compiler_suffix="-gcc-${compiler_version}"
    
    mpi_module=openmpi/${openmpi_version}
    mpi_suffix="-openmpi-${openmpi_version}"

fi

# Jali
jali_install_dir=$NGC/private/jali/${jali_version}${compiler_suffix}${mpi_suffix}
jali_flags="-D WONTON_ENABLE_Jali:BOOL=True -D Jali_ROOT:PATH=$jali_install_dir"

# LAPACKE
lapacke_dir=$NGC/private/lapack/${lapack_version}-patched${compiler_suffix}
lapacke_flags="-D WONTON_ENABLE_LAPACKE:BOOL=True -D LAPACKE_ROOT:PATH=$lapacke_dir"

# Flecsi
flecsi_flags="-D WONTON_ENABLE_FleCSI:BOOL=False"
if [[ $compiler == "gcc6" ]]; then
    flecsi_install_dir=$NGC/private/flecsi/374b56b-gcc-6.4.0
    flecsisp_install_dir=$NGC/private/flecsi-sp/e78c594-gcc-6.4.0
    flecsi_flags="-D WONTON_ENABLE_FleCSI:BOOL=True -D FleCSI_ROOT:PATH=$flecsi_install_dir -D FleCSISP_ROOT:PATH=$flecsisp_install_dir"
fi

# THRUST
thrust_flags=
thrust_suffix=
thrust_install_dir=$NGC/private/thrust/${thrust_version}
if [[ $config_type == "thrust" ]]; then
    thrust_flags="-D WONTON_ENABLE_THRUST=True -D THRUST_ROOT=${thrust_install_dir}"
    thrust_suffix="-thrust"
fi

# Kokkos
kokkos_flags=
kokkos_suffix=
kokkos_install_dir=$NGC/private/kokkos/${kokkos_version}${compiler_suffix}
if [[ $config_type == "kokkos" ]]; then
    kokkos_flags="-D WONTON_ENABLE_Kokkos=True -D Kokkos_ROOT=${kokkos_install_dir}"
    kokkos_suffix="-kokkos"
fi

# MPI or not
mpi_flags="-D WONTON_ENABLE_MPI=True"
if [[ $config_type == "serial" ]]; then
    mpi_flags="-D WONTON_ENABLE_MPI=False"
    mpi_suffix=
    jali_flags=
    flecsi_flags=
fi

# Debug or Optimized build
cmake_build_type=Release
debug_suffix=
if [[ $config_type == "debug" ]]; then
    cmake_build_type=Debug
    debug_suffix="-debug"
fi

# Build up an install dir name
wonton_install_dir=$NGC/private/wonton/${version}${compiler_suffix}${mpi_suffix}${thrust_suffix}${kokkos_suffix}${debug_suffix}


export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load $cxxmodule
if [[ -n "$mpi_flags" ]]; then
    module load ${mpi_module}
fi
module load cmake/3.14.0 # 3.13 or higher is required

module load git

echo "JENKINS WORKSPACE = $WORKSPACE"
cd $WORKSPACE

rm -fr build
mkdir build
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
    $flecsi_flags \
    $lapacke_flags \
    $thrust_flags \
    $kokkos_flags \
    ..
make -j2
ctest -j2 --output-on-failure

if [[ $build_type == "install" ]]; then
    make install
fi
