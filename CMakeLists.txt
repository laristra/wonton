#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

cmake_minimum_required(VERSION 3.13)

project(wonton CXX C)  # need C for FortranC interface in case LAPACKE is ON

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_VERSION_MAJOR GREATER_EQUAL 3.13)
  CMAKE_POLICY(SET CMP0079 NEW)  # allow target_link_libraries to reference
                                 # targets from other directories
endif()

cmake_policy(SET CMP0074 NEW)  # Don't ignore Pkg_ROOT variables

#
# Also when CMake 3.15 is the minimum version we can find_package to
# prefer config files and fall back to module files
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.15)
  set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
endif ()


# SEMANTIC VERSION NUMBERS - UPDATE DILIGENTLY
# As soon as a change with a new version number is merged into the master,
# tag the central repository.

set(WONTON_VERSION_MAJOR 1)
set(WONTON_VERSION_MINOR 2)
set(WONTON_VERSION_PATCH 8)


# Top level target
add_library(wonton INTERFACE)

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(wonton::wonton ALIAS wonton)
set(WONTON_LIBRARIES wonton::wonton CACHE STRING "Wonton library target")


# Find our modules first
if (CMAKE_VERSION GREATER_EQUAL 3.15)
  list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else ()
  set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
endif ()

#------------------------------------------------------------------------------#
# Enable documentation generation
#------------------------------------------------------------------------------#
option(ENABLE_DOXYGEN "Enable Doxygen documentation" OFF)
option(ENABLE_DOXYGEN_WARN "Enable Doxygen warnings" OFF)

if (ENABLE_DOXYGEN)
  include(cmake/doxygen.cmake)
  add_doxygen()
endif()

#------------------------------------------------------------------------------#
# Discover packages here but set them as dependencies of
# wonton_support target in the wonton/support subdirectory. Since
# wonton_support is a dependency of the top level target
# wonton::wonton, projects using wonton will get the transitive
# dependencies
# ------------------------------------------------------------------------------


#------------------------------------------------------------------------------#
# Set up MPI builds
#------------------------------------------------------------------------------#
set(WONTON_ENABLE_MPI False CACHE BOOL "Is MPI enabled in Wonton?") # create cache entry
if (WONTON_ENABLE_MPI)  # if overwritten by user input or environment var
  find_package(MPI REQUIRED)
endif ()


#-----------------------------------------------------------------------------
# FleCSI and FleCSI-SP location
#-----------------------------------------------------------------------------

set(WONTON_ENABLE_FleCSI FALSE CACHE BOOL "FleCSI interface enabled?")
if (WONTON_ENABLE_FleCSI AND NOT FleCSI_LIBRARIES)
  find_package(FleCSI REQUIRED)
  find_package(FleCSISP REQUIRED)

  message(STATUS "FleCSI_LIBRARIES=${FleCSI_LIBRARIES}" )
  message(STATUS "FleCSISP_LIBRARIES=${FleCSISP_LIBRARIES}" )

  target_include_directories(wonton SYSTEM INTERFACE ${FleCSI_INCLUDE_DIR})
  target_include_directories(wonton SYSTEM INTERFACE ${FleCSISP_INCLUDE_DIR})

  target_link_libraries(wonton INTERFACE ${FleCSI_LIBRARIES})
  target_link_libraries(wonton INTERFACE ${FleCSISP_LIBRARIES})

  target_compile_definitions(wonton INTERFACE WONTON_ENABLE_FleCSI)

  if (NOT FleCSI_ROOT)
    get_filename_component(FleCSI_ROOT ${FleCSI_INCLUDE_DIR} DIRECTORY CACHE)
  endif ()
  if (NOT FleCSISP_ROOT)
    get_filename_component(FleCSISP_ROOT ${FleCSISP_INCLUDE_DIR} DIRECTORY CACHE)
  endif ()
endif()


#------------------------------------------------------------------------------#
# Locate Jali
#------------------------------------------------------------------------------#

set(WONTON_ENABLE_Jali False CACHE BOOL "Jali Interface enabled?")
if (WONTON_ENABLE_Jali AND WONTON_ENABLE_MPI AND NOT Jali_LIBRARIES)
  # Look for the Jali package
  
  find_package(Jali REQUIRED)  # specify in Jali_ROOT or CMAKE_PREFIX_PATH

  message(STATUS "Located Jali")
  message(STATUS "Jali_LIBRARIES ${Jali_LIBRARIES}")

  target_link_libraries(wonton INTERFACE ${Jali_LIBRARIES})
  
  if (NOT Jali_ROOT)
    # Jali_CONFIG should be the full path to where the config file was found
    # which is typically SOMEDIR/lib/cmake - back out Jali_ROOT from that
    # until we fix JaliConfig
    get_filename_component(Jali_CONFIGLOC ${Jali_CONFIG} DIRECTORY)
    get_filename_component(Jali_CONFIGLOC_UP1 ${Jali_CONFIGLOC} DIRECTORY)
    get_filename_component(Jali_CONFIGLOC_UP2 ${Jali_CONFIGLOC_UP1} DIRECTORY)
    get_filename_component(Jali_ROOT ${Jali_CONFIGLOC_UP2} DIRECTORY CACHE "Where Jali lives")
  endif ()
  message(STATUS "Jali_ROOT ${Jali_ROOT}")
  
  target_compile_definitions(wonton INTERFACE WONTON_ENABLE_Jali)
endif ()


#-----------------------------------------------------------------------------
# Thrust information
#-----------------------------------------------------------------------------

set(WONTON_ENABLE_THRUST False CACHE BOOL "Is the Thrust library being used?")
if (WONTON_ENABLE_THRUST)

  # Allow for swapping backends
  set(THRUST_HOST_BACKEND "THRUST_HOST_SYSTEM_CPP"  CACHE STRING "Thrust host backend")
  set(THRUST_DEVICE_BACKEND "THRUST_DEVICE_SYSTEM_OMP" CACHE STRING "Thrust device backend")

  if ((${THRUST_HOST_BACKEND} STREQUAL "THRUST_HOST_SYSTEM_OMP") OR
      (${THRUST_DEVICE_BACKEND} STREQUAL "THRUST_DEVICE_SYSTEM_OMP"))
    list(APPEND _components OpenMP)
  endif ()
  if (${THRUST_DEVICE_BACKEND} STREQUAL "THRUST_DEVICE_SYSTEM_CUDA")
    list(APPEND _components CUDA)
  endif ()

  find_package(THRUST COMPONENTS ${_components} REQUIRED MODULE)

  message(STATUS "Enabling compilation with Thrust")
  message(STATUS "Using THRUST_ROOT=${THRUST_ROOT}")
  message(STATUS "THRUST_INCLUDE_DIRS ${THRUST_INCLUDE_DIRS}")

  message(STATUS "Using ${THRUST_HOST_BACKEND} as Thrust host backend")
  message(STATUS "Using ${THRUST_DEVICE_BACKEND} as Thrust device backend")

else ()

  #-----------------------------------------------------------------------------
  # Include Boost libraries (for counting_iterator etc)
  #-----------------------------------------------------------------------------

  find_package(Boost REQUIRED)
  target_include_directories(wonton SYSTEM INTERFACE ${Boost_INCLUDE_DIRS})
  message(STATUS "Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
  #  message(STATUS "Boost_LIBRARIES?=${Boost_LIBRARIES}")

  if (NOT BOOST_ROOT)
    if (Boost_DIR)
      set(BOOST_ROOT ${Boost_DIR} CACHE PATH "Boost installation directories")
    elseif (Boost_INCLUDE_DIRS)
      get_filename_component(BOOST_ROOT ${Boost_INCLUDE_DIRS} DIRECTORY CACHE)
    endif ()
  endif ()
  message(STATUS "BOOST_ROOT ${BOOST_ROOT}")

endif()

set(WONTON_ENABLE_TCMALLOC False CACHE BOOL "Link in tcmalloc?")
if (WONTON_ENABLE_TCMALLOC)  # if overridden by command line
  find_library(TCMALLOC_LIBRARIES NAMES tcmalloc PATHS ${TCMALLOC_ROOT})
  if (NOT TCMALLOC_LIBRARIES)
    message(FATAL_ERROR "TCMALLOC not found")
  endif ()
endif ()



#-----------------------------------------------------------------------------
# Use Kokkos
#-----------------------------------------------------------------------------
set(WONTON_ENABLE_Kokkos FALSE CACHE BOOL "Use Kokkos")
if (WONTON_ENABLE_Kokkos)
  # find package and link against it.
  # conflict with the one packed in Trilinos.
  # ugly hack: override headers and library locations

  find_package(Kokkos REQUIRED)

  # Sadly Kokkos CMake config does not set Kokkos_LIBRARIES so we will set it
  set(Kokkos_LIBRARIES "Kokkos::kokkos" CACHE STRING "Kokkos top level target")

  list(FIND Kokkos_DEVICES OPENMP INDEX)
  set(Kokkos_HAS_OPENMP False)
  if (INDEX GREATER -1)
    set(Kokkos_HAS_OPENMP True)
  endif ()

  list(FIND Kokkos_DEVICES CUDA INDEX)
  set(Kokkos_HAS_CUDA False)
  if (INDEX GREATER -1)
    set(Kokkos_HAS_CUDA True)
  endif ()
  
  message(STATUS "Enabling Kokkos")
  message(STATUS "Kokkos_LIBRARIES = ${Kokkos_LIBRARIES}")

  if (WONTON_ENABLE_Kokkos_CUDA)
    if (NOT Kokkos_HAS_CUDA)
      message(FATAL_ERROR "CUDA enabled for Wonton but Kokkos not compiled with CUDA support")
    endif ()
    
    enable_language(CUDA)

    # Do we need to do this? Should Kokkos take care of this?
    # Also, is there a CMake option to do this? - RVG
    
    # Looks like CUDAFLAGS environment variables or CMAKE_CUDA_FLAGS
    # cmake variable can be used instead of CUDA_NVCC_FLAGS?
    
    # additional options for nvcc:
    # allow '__host__',' __device__' annotations in lambdas.
    # allow host code to invoke '__device__ constexpr' functions,
    # and device code to invoke '__host__ constexpr' functions.
    # disable warning on __global__ defaulted functions.
    string(APPEND CUDA_NVCC_FLAGS " --expt-extended-lambda")
    string(APPEND CUDA_NVCC_FLAGS " --expt-relaxed-constexpr")
    string(APPEND CUDA_NVCC_FLAGS " -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored")

    # nvidia graphic cards architectures:
    # _30 kepler: tesla k40|80, geforce 700, gt-730, support for unified memory.
    # _35 kepler: add support for dynamic parallelism for tesla k40.
    # _37 kepler: add a few more registers for tesla k80.
    # _50 maxwell: tesla/quadro M series.
    # _52 maxwell: quadro m6000, geforce 900, gtx-970|80, gtx titan x.
    # _60 pascal: gp100, tesla p100.
    # _61 pascal: gtx 1030|50|60|70|80, titan xp, tesla P4|40.
    # _70 volta: tesla v100, gtx 1180 (gv104).
    string(APPEND CUDA_NVCC_FLAGS " -arch=sm_30")
    string(APPEND CUDA_NVCC_FLAGS " -gencode=arch=compute_30,code=sm_30")
    string(APPEND CUDA_NVCC_FLAGS " -gencode=arch=compute_35,code=sm_35")
    string(APPEND CUDA_NVCC_FLAGS " -gencode=arch=compute_37,code=sm_37")
    string(APPEND CUDA_NVCC_FLAGS " -gencode=arch=compute_50,code=sm_50")
    string(APPEND CUDA_NVCC_FLAGS " -gencode=arch=compute_52,code=sm_52")
    string(APPEND CUDA_NVCC_FLAGS " -gencode=arch=compute_60,code=sm_60")
    string(APPEND CUDA_NVCC_FLAGS " -gencode=arch=compute_61,code=sm_61")
    string(APPEND CUDA_NVCC_FLAGS " -gencode=arch=compute_70,code=sm_70")
    string(APPEND CUDA_NVCC_FLAGS " -gencode=arch=compute_70,code=compute_70")
    string(APPEND CUDA_NVCC_FLAGS " ")
  endif ()

  # use OpenMP backend on CPU
  if (WONTON_ENABLE_Kokkos_OpenMP)
    if (NOT Kokkos_HAS_OPENMP)
      message(FATAL_ERROR "OPENMP enabled for Wonton but Kokkos not compiled with OPENMP support")
    endif ()
    
    find_package(OpenMP)
  endif ()

  if (NOT Kokkos_ROOT)
    get_property(Kokkos_INCDIR TARGET ${Kokkos_LIBRARIES} PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    if (Kokkos_INCDIR)
      get_filename_component(Kokkos_ROOT ${Kokkos_INCDIR} DIRECTORY CACHE)
    endif ()
  endif ()
  message(STATUS "Kokkos_ROOT ${Kokkos_ROOT}")
endif ()


#------------------------------------------------------------------------------#
# Configure LAPACKE
#------------------------------------------------------------------------------#
set(WONTON_ENABLE_LAPACKE False CACHE BOOL "Use LAPACKE Solvers")
if (WONTON_ENABLE_LAPACKE)  # if overridden by command line or environment
  # Depending on the version we may find LAPACK config file (3.5.0 or older)
  # or LAPACKE config file (new)
  set(LAPACKE_CONFIG_FOUND False CACHE STRING "Was a LAPACKE config file found?")
  set(LAPACK_CONFIG_FOUND False CACHE STRING "Was a LAPACK config file found?")
  
  # Find lapacke-config.cmake or lapackeconfig.cmake
  find_package(LAPACKE QUIET
    CONFIG
    NAMES LAPACKE lapacke)

  if (LAPACKE_FOUND)
    set(LAPACKE_CONFIG_FOUND True)
  else ()
    # Maybe it's an older version (3.5.0 or older) that only puts out
    # lapack-config.cmake. When looking for this config file, also add
    # LAPACKE_ROOT to the list of directories that are searched. This
    # version sets the transitive dependencies for the target
    # "lapacke" correctly so we don't have to explicitly tell it to
    # link in lapack or blas

    find_package(LAPACK
      CONFIG
      NAMES LAPACK lapack
      PATHS ${LAPACKE_ROOT})

    if (LAPACK_FOUND)
      # is there a target named lapacke?
      if (TARGET lapacke)
	set(LAPACKE_LIBRARIES lapacke)
	set(LAPACKE_FOUND True)
	set(LAPACK_CONFIG_FOUND True)
      endif ()
    endif ()
    
  endif ()
endif ()

if (WONTON_ENABLE_LAPACKE)
  if (LAPACKE_FOUND)
    enable_language(Fortran)
    include(FortranCInterface)  # will ensure the fortran library is linked in
    
    message(STATUS "LAPACKE_FOUND ${LAPACKE_FOUND}")
    message(STATUS "LAPACKE_LIBRARIES  ${LAPACKE_LIBRARIES}")
    
    if (NOT LAPACKE_ROOT)
      get_property(LAPACKE_LOC TARGET ${LAPACKE_LIBRARIES} PROPERTY LOCATION)
      get_filename_component(LAPACKE_LIBDIR ${LAPACKE_LOC} DIRECTORY)
      get_filename_component(LAPACKE_ROOT ${LAPACKE_LIBDIR} DIRECTORY CACHE "Where LAPACKE lives")
    endif ()
    message(STATUS "LAPACKE_ROOT ${LAPACKE_ROOT}")
  else ()
    message(FATAL_ERROR "LAPACKE enabled but not found.")
  endif ()
endif ()


if (ENABLE_UNIT_TESTS)

  # This needs to be set in the root directory for tests to be run by
  # 'make test' or ctest
  enable_testing()

  include(cmake/unittest.cmake)
  
endif ()



#-----------------------------------------------------------------------------
# Recurse down the source directories building up dependencies
#-----------------------------------------------------------------------------

add_subdirectory(wonton)

# In addition to the include directories of the source, we need to
# include the build or directory to get the autogenerated
# wonton-config.h (The first of these is needed if Wonton is included
# as a submodule, the second is needed for the auto-generated config
# file if Wonton is included as a submodule, the third is to get the
# autogenerated config header if Wonton is being compiled separately
# and the last is for dependencies in installations)

target_include_directories(wonton INTERFACE
  $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
  $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

# Wonton targets

install(TARGETS wonton
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )

#-----------------------------------------------------------------------------
# Add any applications built upon wonton
#-----------------------------------------------------------------------------

add_subdirectory(app)

#-----------------------------------------------------------------------------
# Prepare output for configuration files to be used by projects importing Wonton
#-----------------------------------------------------------------------------

# Write a configuration file from template replacing only variables enclosed
# by the @ sign.

configure_file(${PROJECT_SOURCE_DIR}/cmake/wontonConfig.cmake.in 
  wontonConfig.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/wontonConfig.cmake
  DESTINATION lib/cmake/wonton)


# write out a version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(wontonConfigVersion.cmake
  VERSION "${WONTON_MAJOR_VERSION}.${WONTON_MINOR_VERSION}.${WONTON_PATCH_VERSION}"
  COMPATIBILITY SameMajorVersion)
install(FILES ${PROJECT_BINARY_DIR}/wontonConfigVersion.cmake
  DESTINATION lib/cmake/wonton)


# export targets

install(EXPORT wonton_LIBRARIES
  FILE wontonTargets.cmake
  NAMESPACE wonton::
  EXPORT_LINK_INTERFACE_LIBRARIES
  DESTINATION lib/cmake/wonton)



# Dynamically configured header files that contains defines like
# WONTON_ENABLE_MPI etc. if enabled.

configure_file(${PROJECT_SOURCE_DIR}/config/wonton-config.h.in
  ${PROJECT_BINARY_DIR}/wonton-config.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/wonton-config.h
  DESTINATION include)


# Install the FindTHRUST module needed by downstream packages when
# processing wontonConfig.cmake
install(FILES ${PROJECT_SOURCE_DIR}/cmake/FindTHRUST.cmake
  DESTINATION lib/cmake/wonton/modules)
