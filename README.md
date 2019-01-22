
# Wonton 

Wonton library provides low level infrastructure for wrappers 
to various unstructured mesh and data (state) managers as 
well as various math and geometry algorithms and data structures.

The primary reason for Wonton's existence is to maintain a single
and common source of such low level types for any other dependent 
libraries such as Portage and Tangram.  

## Getting Started 
To obtain a copy of wonton and its submodules from GitHub, clone 
recursively:
```sh
git clone --recursive https://github.com/laristra/wonton
```

### Prerequisites
Wonton uses standard C++11 features, so a fairly modern compiler 
is needed. We regularly test with Intel 17+ or GCC 5.3+. The build 
system _requires_ CMake version 3.0+. 

The following libraries are also _required_:

- LAPACKE (3.8.0+)

- **__Either__** Boost (1.68.0+) **__or__** Thrust (1.6.0+):
  We wrap some features of either one of these packages.  If you would
  like to run with OpenMP or TBB threads, then you _must_ use Thrust.

Wonton provides wrappers for a few third-party mesh types.  Building
support for these is _optional_:

- [Jali](http://github.com/lanl/jali):

  We regularly test with verison 1.0.0.  You will need to set the
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
  `374b56b` of the FleCSI _stable_ branch, and hash `e78c594` of the
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

This compiles the serial code and about a dozen application tests.  To
run the tests, simply execute

```sh
wonton/build $ ctest
```

If you wish to install the code into the `CMAKE_INSTALL_PREFIX` then
simply execute
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

## Darwin

Execute the following from the wonton root directory:

```c++
# machine=darwin-fe
module load openmpi/2.1.2-intel_17.0.6 boost/1.58.0 cmake
setenv JALI_INSTALL_PREFIX /usr/projects/ngc/private/jali/1.0.0-intel-17.0.6-openmpi-2.1.2
setenv TPL_INSTALL_PREFIX /usr/projects/ngc/private/jali-tpl/1.1.0-intel-17.0.6-openmpi-2.1.2
export LAPACKE_DIR=/usr/projects/ngc/private/lapack/3.8.0-patched-intel-17.0.6
mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE=Debug \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_MPI=True \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    -D Boost_INCLUDE_DIR:PATH=$TPL_INSTALL_PREFIX/include \
    -D LAPACKE_DIR=$LAPACKE_DIR \
    ..
make -j16
ctest -j16 --output-on-failure
```

## Varan

Execute the following from the wonton root directory:

```c++
# machine=varan
echo "This is a test on varan"
```
