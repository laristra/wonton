#!/usr/bin/env bash
#
# Script to build the doxygen docs.  Here, we don't really care about which
# compiler or build options are set.

# Exit on error
set -e

# Echo each command
set -x

# Set umask so installations have correct permissions
umask 002

JALI_VERSION=1.0.0
openmpi_version=2.1.2

# location on XLAN
NGC_DIR=/usr/local/codes/ngc

JALI_INST=${NGC_DIR}/private/jali/${JALI_VERSION}-gcc-6.4.0-openmpi-${openmpi_version}


git config user.email ""
git config user.name "Jenkins"
git merge origin/master

export SHELL=/bin/sh
export MODULEPATH=""

# Setup modules
. /opt/local/packages/Modules/default/init/sh
module load gcc/6.4.0 openmpi/${openmpi_version} cmake

# the system doxygen and LaTeX are too old; use these instead
export PATH=/usr/local/codes/ngc/home/cmalone/texlive/2016/bin/x86_64-linux/:$PATH
DOXY_EXE=/home/cmalone/code/doxygen/build_1_8_13/bin/doxygen

echo $WORKSPACE
cd $WORKSPACE

mkdir build
cd build

cmake \
    -D Jali_DIR:FILEPATH=${JALI_INST}/lib \
    -D ENABLE_DOXYGEN=True \
    -D DOXYGEN_EXECUTABLE=$DOXY_EXE \
    ..
make doxygen

# build a PDF copy of the docs
pushd doc/doxygen/latex
make pdf
popd
