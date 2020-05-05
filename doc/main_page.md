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
- A C++-14 compatible compiler; regular testing is performed with GCC
  6.3+ and Intel 17+.
- CMake 3.0+
- LAPACKE (3.8.0+)
- Boost (1.58.0+) **or** Thrust (1.6.0+)

By design, Wonton wrappers do not depend on MPI. However, testing for 
some wrappers such as [Jali](https://github.com/lanl/jali) requires MPI. 
For such cases, MPI can be enabled in Wonton by setting the CMake varible 
`WONTON_ENABLE_MPI=True`.

On-node parallelism in Portage and Tangram are exposed through aliases 
of the transform operator in the [Thrust](https://thrust.github.io) library
which is defined in Wonton. In order to use 
[Thrust](https://thrust.github.io) based transform operators, Wonton
should be build with at least two CMake variables:
`WONTON_ENABLE_THRUST=True` and `THRUST_ROOT=<path_to_thrust_directory>` (alternatively, the path to Thrust can be included in `CMAKE_PREFIX_PATH`).
Additionally, one can specify the Thrust backend to utilize, with the
default being the OpenMP backend
`THRUST_BACKEND=THRUST_DEVICE_SYSTEM_OMP`.  Wonton also supports the
`THRUST_DEVICE_SYSTEM_TBB` backend.  Regular testing happens with
Thrust 1.8.

**If you turn on Thrust for multi-threading-enabled executables, the team strongly
recommends linking to the TCMalloc library available in 
[Google Performance Tools](https://github.com/gperftools/gperftools) 
to see the expected scaling.**

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
					   -D CMAKE_PREFIX_PATH="/path/to/FleCSI/install;/path/to/FleCSI-sp/install" \
					   -D WONTON_ENABLE_MPI=True \
					   -D WONTON_ENABLE_THRUST=True 
             		   -D THRUST_ROOT=/path/to/thrust/include/directory \
					   -D WONTON_ENABLE_Jali=True \
					   -D Jali_ROOT=path/ to/Jali/lib \
					   -D WONTON_ENABLE_FleCSI=True \
					   -D ENABLE_DOXYGEN=True \
					   -D WONTON_ENABLE_LAPACKE=/path/to/LAPACKE/ \
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
| `ENABLE_DOXYGEN:BOOL` | Create a target to build this documentation | `False` |
| `WONTON_ENABLE_MPI:BOOL` | Build with support for MPI | `False` |
| `ENABLE_TCMALLOC:BOOL` | Build with support for TCMalloc | `False` |
| `WONTON_ENABLE_THRUST:BOOL` | Turn on Thrust support for on-node parallelism | `False` |
| `THRUST_ROOT:PATH` | Directory of the Thrust install | "" |
| `THRUST_BACKEND:STRING` | Backend to use for Thrust | `"THRUST_DEVICE_SYSTEM_OMP"` |
| `ENABLE_UNIT_TESTS:BOOL` | Turn on compilation and test harness of unit tests | `False` |
| `WONTON_ENABLE_Jali:BOOL` | Turn on Jali mesh/state wrappers | `False` |
| `Jali_ROOT:PATH` | Where to find Jali. This version of Wonton works with version 1.1.4+ of Jali | "" |
| `WONTON_ENABLE_FleCSI:BOOL` | Turn on support for the FleCSI Burton specialization; must set `CMAKE_PREFIX_PATH` to a location where _both_ FleCSI and FleCSI-SP can be found. Both FleCSI packages are under constant development.  This version of wonton is known to work with hash `374b56b` of the FleCSI _stable_ branch, and hash `e78c594` of the FleCSI-SP _stable_ branch. | `False` |
| `WONTON_ENABLE_LAPACKE:BOOL` | Turn on LAPACKE solvers | `False` |
| `LAPACKE_ROOT:PATH` | Where to find LAPACKE. | "" |
