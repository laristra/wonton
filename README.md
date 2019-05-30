
# Wonton 

Wonton library provides low level infrastructure for wrappers 
to various unstructured mesh and data (state) managers.
A wrapper is an implementation of an interface to the native
interface of the mesh and state manager and can be used by 
libraries such as Portage and Tangram without copying data
from the mesh and state manager. 

The primary reason for Wonton's existence is to enable applications
to interface to the Portage and Tangram libraries. Wonton also provides
a common repository for these interface implementations as well
as common definitions and methods. 

## Getting Started 
To obtain a copy of wonton and its submodules from GitHub, clone 
recursively:
```sh
git clone --recursive https://github.com/laristra/wonton
```

If you are familiar with Docker, take a look at
our
[Dockerfile](https://github.com/laristra/wonton/blob/master/docker/Dockerfile) for
a working build environment.  In particular, the Dockerfile builds off
of
the [wonton-buildenv](https://github.com/laristra/wonton-buildenv)
Dockerfile, and uses
our
[travis.yml](https://github.com/laristra/wonton/blob/master/.travis.yml) file
with Travis CI.

### Prerequisites
Wonton uses standard C++11 features, so a fairly modern compiler 
is needed. We regularly test with Intel 18.0.1, GCC 6.4.0, and GCC 7.3.0. 
The build system _requires_ CMake version 3.0+. 

The following libraries are also _required_:

- LAPACKE (3.8.0+)

- **__Either__** Boost (1.58.0+) **__or__** Thrust (1.6.0+):
  We wrap some features of either one of these packages.  If you would
  like to run with OpenMP or TBB threads, then you _must_ use Thrust.

Wonton provides wrappers for a few third-party mesh types.  Building
support for these is _optional_:

- [Jali](http://github.com/lanl/jali):

  We regularly test with verison 1.0.4.  You will need to set the
  `Jali_Dir` CMake variable if you wish to build support for Jali and
  its tests (see examples below).

- [FleCSI Burton Specialization](http://github.com/laristra/flecsi-sp):

  The Burton specialization in the `flecsi-sp` repository is built on
  top of [FleCSI](http://github.com/laristra/flecsi).  You will need
  _both_ projects to build support for the Burton mesh specialization
  and its tests.  You will need to set `ENABLE_FleCSI=True` and add
  the FleCSI and FleCSI-sp install paths to the `CMAKE_PREFIX_PATH`;
  see examples below.  Both FleCSI packages are under constant
  development.  This version of wonton is known to work with hash
  `bd29de5d` of the FleCSI _stable_ branch, and hash `e78c594` of the
  FleCSI-SP _stable_ branch.

The [documentation](http://wonton.lanl.gov) is built using doxygen (1.8+). 

### Installing

In the simplest case where you have the appropriate versions mentioned
above and Boost and LAPACKE are in the usual locations that CMake
searches, then the build step is:

```sh
wonton $ mkdir build
wonton $ cd build
wonton/build $ cmake -DENABLE_UNIT_TESTS=True ..
wonton/build $ make
```

This compiles the serial code and about a dozen unit tests.  To
run the tests, execute

```sh
wonton/build $ ctest
```

If you wish to install the code into the `CMAKE_INSTALL_PREFIX` then
execute
```sh
wonton/build $ make install
```

To build the documentation, one would configure with the
`-DENABLE_DOXYGEN=True` flag, and then `make doxygen`.

# License

This project is licensed under a modified 3-clause BSD license - see
the [LICENSE](https://github.com/laristra/wonton/blob/master/LICENSE)
file for details.

# Release

This software has been approved for open source release and has been
assigned **LA-CC-18-019**.

----
----

# Example builds

Below we list copy & paste instructions for several local machines; we
have a script that parses this README file to execute the examples
below to ensure they build.

## Varan

Execute the following from the Wonton root directory:

```c++
# machine=varan
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load intel/18.0.1 openmpi/2.1.2 cmake/3.10.2
NGC_INCLUDE_DIR=/usr/local/codes/ngc/private/include
TPL_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali-tpl/1.2.0-intel-18.0.1-openmpi-2.1.2
JALI_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali/1.0.4-intel-18.0.1-openmpi-2.1.2
LAPACKE_INSTALL_PREFIX=/usr/local/codes/ngc/private/lapack/3.8.0-patched-intel-18.0.1
mkdir build
cd build
cmake \
  -D CMAKE_C_COMPILER=`which mpicc` \
  -D CMAKE_CXX_COMPILER=`which mpiCC` \
  -D CMAKE_BUILD_TYPE=Release \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D LAPACKE_DIR:FILEPATH=$LAPACKE_INSTALL_PREFIX \
  -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
  -D NGC_INCLUDE_DIR:FILEPATH=$NGC_INCLUDE_DIR \
  ..
make -j2
ctest --output-on-failure
```
