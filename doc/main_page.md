# Welcome to Wonton!   {#mainpage}

Wonton is primarily a library of interfaces for accessing mesh and state (mesh field) 
data of an application without explicitly copying or converting that data for use 
in another application. The mesh and state wrappers in Wonton allow packages such
as the [Portage](https://github.com/laristra/portage) remapping package and the 
[Tangram](https://github.com/laristra/tangram) interface reconstruction package 
to easily interface with any application with minimal effort and for the data 
to be read and written efficiently.  

To summarize, Wonton is a low-level library that provides:
- wrappers for mesh and state data managers
- geometric and algebraic type definitions 
- standard math algorithms such as svd, least-squares fittings, etc. 

Wonton also provides access to the [R3D](https://github.com/devonmpowell/r3d.git) 
library for polyhedral intersection methods as a git submodule. 

For further high-level discussion of the methods and data structures provided
within Wonton, please see the [Concepts](@ref concepts) page. 

---

# Details and Requirements

At a minimum, wonton requires:
- A C++-11 compatible compiler; regular testing is performed with GCC
  5.3+ and Intel 17+.
- CMake 3.0+
- LAPACKE (3.7.1+)
- Boost (1.53.0+) **or** Thrust (1.6.0+)

By design, Wonton wrappers do not depend on MPI. However, testing for 
some wrappers such as Jali requires MPI. For such cases, MPI can be 
enabled in Wonton by setting the CMake varible 
`ENABLE_MPI=True`.

On-node parallelism is exposed through
the [Thrust](https://thrust.github.io) library.  Enabling Thrust
within Wonton requires setting at least two CMake variables:
`ENABLE_THRUST=True` and `THRUST_DIR=<path_to_thrust_directory>`.
Additionally, one can specify the Thrust backend to utilize, with the
default being the OpenMP backend
`THRUST_BACKEND=THRUST_DEVICE_SYSTEM_OMP`.  Wonton also supports the
`THRUST_DEVICE_SYSTEM_TBB` backend.  Regular testing happens with
Thrust 1.8. 

## Obtaining wonton

The latest release of [wonton](https://github.com/laristra/wonton)
lives on GitHub.  Wonton makes use of git submodules, so it must be
cloned recursively:

```sh
git clone --recursive https://github.com/laristra/wonton
```

## Building

Wonton uses the CMake build system. A build with MPI, Thrust (for on-node
parallelism), unit test support, documentation support, and support for 
both [Jali](https://github.com/lanl/jali) and the Burton 
[specialization](https://github.com/laristra/flecsi-sp) of 
the [FleCSI](https://github.com/laristra/flecsi) library, assuming that CMake knows where to
find your Boost and LAPACKE installations,  would look like:

~~~sh
wonton/ $ mkdir build
wonton/ $ cd build
wonton/build/ $ cmake -D ENABLE_UNIT_TESTS=True \
					   -D ENABLE_MPI=True \
					   -D ENABLE_THRUST=True 
             				   -D THRUST_DIR=/path/to/thrust/include/directory \
					   -D Jali_DIR=path/ to/Jali/lib \
					   -D ENABLE_FleCSI=True \
					   -D CMAKE_PREFIX_PATH="/path/to/FleCSI/install;/path/to/FleCSI-sp/install" \
					   -D ENABLE_DOXYGEN=True \
					   -D LAPACKE_DIR=/path/to/LAPACKE/ \
					   ..
wonton/build/ $ make           # builds the library and tests
wonton/build/ $ ctest          # runs the tests
wonton/build/ $ make doxygen   # builds this HTML and a PDF form of the documentation
wonton/build/ $ make install   # installs the wonton library and headers into CMAKE_INSTALL_PREFIX
~~~

## Useful CMake Flags
Below is a non-exhaustive list of useful CMake flags for building
Wonton.

| CMake flag:type | Description | Default |
|:----------|:------------|:--------|
| `CMAKE_BUILD_TYPE:STRING`| `Debug` or optimized `Release` build | `Debug` |
| `CMAKE_INSTALL_PREFIX:PATH` | Location for the wonton library and headers to be installed | `/usr/local` |
| `CMAKE_PREFIX_PATH:PATH` | Locations where CMake can look for packages; needs to be set to the FleCSI and FleCSI-SP locations if using FleCSI | "" |
| `ENABLE_APP_TESTS:BOOL` | Turn on compilation and test harness of application tests | `False` |
| `ENABLE_DOXYGEN:BOOL` | Create a target to build this documentation | `False` |
| `ENABLE_MPI:BOOL` | Build with support for MPI | `False` |
| `ENABLE_TCMALLOC:BOOL` | Build with support for TCMalloc | `False` |
| `ENABLE_THRUST:BOOL` | Turn on Thrust support for on-node parallelism | `False` |
| `ENABLE_UNIT_TESTS:BOOL` | Turn on compilation and test harness of unit tests | `False` |
| `Jali_DIR:PATH` | Hint location for CMake to find Jali.  This version of wonton works with version 1.0.0 of Jali | "" |
| `ENABLE_FleCSI:BOOL` | Turn on support for the FleCSI Burton specialization; must set `CMAKE_PREFIX_PATH` to a location where _both_ FleCSI and FleCSI-SP can be found. Both FleCSI packages are under constant development.  This version of wonton is known to work with hash `374b56b` of the FleCSI _stable_ branch, and hash `e78c594` of the FleCSI-SP _stable_ branch. | `False` |
| `LAPACKE_DIR:PATH` | Hint location for CMake to find LAPACKE include files if `pkg_config` fails. | "" |
| `TCMALLOC_LIB:PATH` | The TCMalloc library to use | `${HOME}` |
| `THRUST_DIR:PATH` | Directory of the Thrust install | "" |
| `THRUST_BACKEND:STRING` | Backend to use for Thrust | `"THRUST_DEVICE_SYSTEM_OMP"` |
