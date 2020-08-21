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

echo "inside build_matrix entry"
# special case for README builds
if [[ $build_type == "readme" ]]; then
  python2 $WORKSPACE/jenkins/parseREADME.py $WORKSPACE/README.md $WORKSPACE
  exit
fi

# set modules and install paths

jali_version=1.1.4
lapack_version=3.8.0

export NGC=/usr/projects/ngc
ngc_include_dir=$NGC/private/include

thrust_dir=${ngc_include_dir}


# compiler-specific settings
if [[ $compiler == "intel18" ]]; then
  compiler_version=18.0.5
  cxxmodule=intel/${compiler_version}
  compiler_suffix="-intel-${compiler_version}"

  openmpi_version=2.1.2
  mpi_module=openmpi/2.1.2
  mpi_suffix="-openmpi-${openmpi_version}"

elif [[ $compiler =~ "gcc" ]]; then
    
  if [[ $compiler == "gcc6" ]]; then
      compiler_version=6.4.0
  elif [[ $compiler == "gcc7" ]]; then
      compiler_version=7.4.0
  fi  
  cxxmodule=gcc/${compiler_version}
  compiler_suffix="-gcc-${compiler_version}"
  
  openmpi_version=2.1.2
  mpi_module=openmpi/${openmpi_version}
  mpi_suffix="-openmpi-${openmpi_version}"
  
fi
jali_install_dir=$NGC/private/jali/${jali_version}${compiler_suffix}${mpi_suffix}
lapacke_dir=$NGC/private/lapack/${lapack_version}-patched${compiler_suffix}


# build-type-specific settings
cmake_build_type=Release
mpi_flags="-D WONTON_ENABLE_MPI=True"
extra_flags=
thrust_flags=
jali_flags="-D WONTON_ENABLE_Jali=True -D Jali_ROOT:FILEPATH=$jali_install_dir"
lapacke_flags="-D WONTON_ENABLE_LAPACKE=True -D LAPACKE_ROOT:FILEPATH=$lapacke_dir"
if [[ $build_type == "debug" ]]; then
    cmake_build_type=Debug
elif [[ $build_type == "serial" ]]; then
    mpi_flags=
    jali_flags=         # jali is not available in serial
elif [[ $build_type == "thrust" ]]; then
    thrust_flags="-D WONTON_ENABLE_THRUST=True -DTHRUST_ROOT:FILEPATH=${thrust_dir}"
fi

# We don't seem to have a suitable (read as "ancient") flecsi install on HPC 
#if [[ $compiler == "gcc6" && $build_type != "serial" ]]; then
#    flecsi_install_dir=$NGC/private/flecsi/374b56b-gcc-6.4.0
#    flecsisp_install_dir=$NGC/private/flecsi-sp/e78c594-gcc-6.4.0
#    flecsi_flags="-D WONTON_ENABLE_FleCSI:BOOL=True -D FleCSI_ROOT:PATH=$flecsi_install_dir -D FleCSISP_ROOT:PATH=$flecsisp_install_dir"
#fi

export SHELL=/bin/sh

export MODULEPATH=""
. /usr/share/lmod/lmod/init/sh

if ! module avail -t 2>&1 | grep -wq "^${cxxmodule}" ; then 
    echo "Module ${cxxmodule} not found"
    exit 2
fi
module load $cxxmodule

if [[ -n "$mpi_flags" ]]; then
    if ! module avail -t 2>&1 | grep -wq "^${mpi_module}" ; then 
        echo "Module ${mpi_module} not found"
        exit 2
    fi
    module load ${mpi_module}
fi

if ! module avail -t 2>&1 | grep -wq "cmake/3.14.0" ; then 
    echo "Module cmake/3.14.0 not found"
    exit 2
fi
module load cmake/3.14.0 # 3.13 or higher is required

echo $WORKSPACE
cd $WORKSPACE

mkdir build
cd build

cmake \
    -D CMAKE_BUILD_TYPE=$cmake_build_type \
    -D CMAKE_CXX_FLAGS="-Wall -Werror" \
    -D CMAKE_PREFIX_PATH=$ngc_include_dir \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_JENKINS_OUTPUT=True \
    $mpi_flags \
    $jali_flags \
    $flecsi_flags \
    $lapacke_flags \
    $thrust_flags \
    ..
make -j2
ctest -j8 --output-on-failure
