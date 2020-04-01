#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

cmake_minimum_required(VERSION 3.13)

project(wonton CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

cinch_minimum_required(VERSION 1.0)

if (CMAKE_VERSION_MAJOR GREATER_EQUAL 3.13)
  CMAKE_POLICY(SET CMP0079 NEW)  # allow target_link_libraries to reference
                                 # targets from other directories
endif()

cmake_policy(SET CMP0074 NEW)  # Don't ignore Pkg_ROOT variables



# SEMANTIC VERSION NUMBERS - UPDATE DILIGENTLY
# As soon as a change with a new version number is merged into the master,
# tag the central repository.

set(WONTON_VERSION_MAJOR 1)
set(WONTON_VERSION_MINOR 1)
set(WONTON_VERSION_PATCH 3)


# Top level target
add_library(wonton INTERFACE)

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(wonton::wonton ALIAS wonton)
set(WONTON_LIBRARIES wonton::wonton CACHE STRING "Wonton library target")


# cinch extras

cinch_load_extras()

set(CINCH_HEADER_SUFFIXES "\\.h")

# Explicitly turn off module paths inherited from Cinch. When we get rid of
# Cinch we can delete this line
set(CMAKE_MODULE_PATH "")

# Find our modules first
if (CMAKE_VERSION GREATER_EQUAL 3.15)
  list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else ()
  set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
endif ()


#------------------------------------------------------------------------------#
# Discover packages here but set them as dependencies of
# wonton_support target in the wonton/support subdirectory. Since
# wonton_support is a dependency of the top level target
# wonton::wonton, projects using wonton will get the transitive
# dependencies
# ------------------------------------------------------------------------------


#------------------------------------------------------------------------------#
# Set up MPI builds
# (eventually most of this should be pushed down into cinch)
#------------------------------------------------------------------------------#
set(WONTON_ENABLE_MPI False CACHE BOOL "Is MPI enabled in Wonton?")
if (ENABLE_MPI)
  find_package(MPI REQUIRED)
  set(WONTON_ENABLE_MPI True)
endif (ENABLE_MPI)


#-----------------------------------------------------------------------------
# FleCSI and FleCSI-SP location
#-----------------------------------------------------------------------------

set(ENABLE_FleCSI FALSE CACHE BOOL "Use FleCSI")
if (ENABLE_FleCSI AND NOT FleCSI_LIBRARIES)
  find_package(FleCSI REQUIRED)
  find_package(FleCSISP REQUIRED)

  set(WONTON_ENABLE_FleCSI True CACHE BOOL "FleCSI interface enabled?")

  message(STATUS "FleCSI_LIBRARIES=${FleCSI_LIBRARIES}" )
  message(STATUS "FleCSISP_LIBRARIES=${FleCSISP_LIBRARIES}" )

  target_include_directories(wonton INTERFACE ${FleCSI_INCLUDE_DIR})
  target_include_directories(wonton INTERFACE ${FleCSISP_INCLUDE_DIR})

  target_link_libaries(wonton INTERFACE ${FleCSI_LIBRARIES})
  target_link_libraries(wonton INTERFACE ${FleCSISP_LIBRARIES})

  target_compile_definitions(wonton INTERFACE WONTON_ENABLE_FleCSI)
endif(ENABLE_FleCSI AND NOT FleCSI_LIBRARIES)


#------------------------------------------------------------------------------#
# Configure Jali
# (this includes the TPLs that Jali will need)
#------------------------------------------------------------------------------#

if (ENABLE_Jali AND ENABLE_MPI AND NOT Jali_LIBRARIES)
  # Look for the Jali package
  
  find_package(Jali REQUIRED)  # specify in Jali_ROOT or CMAKE_PREFIX_PATH

  set(WONTON_ENABLE_Jali True CACHE BOOL "Jali interface enabled?")

  message(STATUS "Located Jali")
  message(STATUS "Jali_LIBRARIES ${Jali_LIBRARIES}")

  target_link_libraries(wonton INTERFACE ${Jali_LIBRARIES})
  
  if (NOT Jali_ROOT)
    # Jali_CONFIG should be the full path to where the config file was found
    # which is typically SOMEDIR/lib/cmake - back out Jali_ROOT from that
    # until we fix JaliConfig
    set(Jali_ROOT ${Jali_CONFIG}/../.. CACHE FILEPATH "Where Jali lives")
  endif ()
  
  target_compile_definitions(wonton INTERFACE WONTON_ENABLE_Jali)
endif ()


#-----------------------------------------------------------------------------
# Thrust information
#-----------------------------------------------------------------------------
if (ENABLE_THRUST)   # if it is overridden by the command line

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

  message(STATUS "Using ${THRUST_HOST_BACKEND} as Thrust host backend")
  message(STATUS "Using ${THRUST_DEVICE_BACKEND} as Thrust device backend")


  set(WONTON_ENABLE_THRUST True CACHE BOOL "Is the Thrust library being used?" FORCE)

else ()

  #-----------------------------------------------------------------------------
  # Include Boost libraries (for counting_iterator etc)
  #-----------------------------------------------------------------------------

  find_package(Boost REQUIRED)
  target_include_directories(wonton SYSTEM PUBLIC ${Boost_INCLUDE_DIR})
  message(STATUS "Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIR}")

endif()

if (ENABLE_TCMALLOC)
  find_library(TCMALLOC_LIBRARIES NAMES tcmalloc PATHS ${TCMALLOC_ROOT})
  if (TCMALLOC_LIBRARIES)
    set(WONTON_ENABLE_TCMALLOC True CACHE BOOL "Link in tcmalloc?")
  endif ()
endif ()



#-----------------------------------------------------------------------------
# Use Kokkos
#-----------------------------------------------------------------------------
set(ENABLE_Kokkos FALSE CACHE BOOL "Use Kokkos")
if (ENABLE_Kokkos)
  # find package and link against it.
  # conflict with the one packed in Trilinos.
  # ugly hack: override headers and library locations
  find_path(KOKKOS_INCLUDE_DIR Kokkos_Core.hpp HINTS "${Kokkos_DIR}/include")
  find_library(KOKKOS_CORE_LIBRARY NAMES kokkoscore
      HINTS "${Kokkos_DIR}/lib" "${Kokkos_DIR}/lib64")
  find_library(KOKKOS_CONTAINERS_LIBRARY NAMES kokkoscontainers
      HINTS "${Kokkos_DIR}/lib" "${Kokkos_DIR}/lib64")
  set(KOKKOS_LIBRARIES ${KOKKOS_CORE_LIBRARY} ${KOKKOS_CONTAINERS_LIBRARY})
  #find_package(Kokkos REQUIRED HINTS "${Kokkos_DIR}")

  message(STATUS "Enabling Kokkos")
  set(WONTON_ENABLE_KOKKOS True CACHE BOOL "Whether Kokkos is enabled")
  include_directories("${KOKKOS_INCLUDE_DIR}")
  list(APPEND WONTON_EXTRA_LIBRARIES "${KOKKOS_LIBRARIES}")

  if (USE_GPU)
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

    # add to nvcc flags
    string(APPEND CMAKE_CXX_FLAGS "${CUDA_NVCC_FLAGS}")
    message(STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}")
  endif (USE_GPU)

  # use OpenMP backend on CPU
  find_package(OpenMP)
  if (OPENMP_FOUND)
    string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
    string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " ${OpenMP_EXE_LINKER_FLAGS}")
  endif (OPENMP_FOUND)

  # use topology.
  find_package(Hwloc)
  if (HWLOC_FOUND)
    include_directories("${HWLOC_INCLUDE_DIR}")
    list(APPEND WONTON_EXTRA_LIBRARIES "${HWLOC_LIBRARY}")
  endif (HWLOC_FOUND)
endif (ENABLE_Kokkos)


#------------------------------------------------------------------------------#
# Configure LAPACKE
#------------------------------------------------------------------------------#

if (ENABLE_LAPACKE)
  # Depending on the version we may find LAPACK config file (3.5.0 or older)
  # or LAPACKE config file (new)
  set(LAPACKE_CONFIG_FOUND False CACHE STRING "Was a LAPACKE config file found?")
  set(LAPACK_CONFIG_FOUND False CACHE STRING "Was a LAPACK config file found?")
  
  # Find lapacke-config.cmake or lapackeconfig.cmake
  find_package(LAPACKE
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

if (LAPACKE_FOUND)
  enable_language(Fortran)
  include(FortranCInterface)  # will ensure the fortran library is linked in
  
  set(WONTON_ENABLE_LAPACKE True CACHE BOOL "LAPACKE libraries linked in?")

  target_include_directories(wonton INTERFACE ${LAPACKE_INCLUDE_DIRS})
  target_compile_definitions(wonton INTERFACE WONTON_HAS_LAPACKE)

  target_link_libraries(wonton INTERFACE ${LAPACKE_LIBRARIES})

  message(STATUS "LAPACKE_FOUND ${LAPACKE_FOUND}")
  message(STATUS "LAPACKE_LIBRARIES  ${LAPACKE_LIBRARIES}")
else ()
  message(FATAL_ERROR "LAPACKE enabled but not found.")
endif (LAPACKE_FOUND)


#-----------------------------------------------------------------------------
# Thrust information
#-----------------------------------------------------------------------------
set(ENABLE_THRUST FALSE CACHE BOOL "Use Thrust")
if(ENABLE_THRUST)
  message(STATUS "Enabling compilation with Thrust")
  # allow the user to specify a THRUST_DIR, otherwise use ${NGC_INCLUDE_DIR}
  # NOTE: thrust internally uses include paths from the 'root' directory, e.g.
  #
  #       #include "thrust/device_vector.h"
  #
  #       so the path here should point to the directory that has thrust as
  #       a subdirectory.
  # Use THRUST_DIR directly if specified, otherwise try to build from NGC
  set(THRUST_DIR "${NGC_INCLUDE_DIR}" CACHE PATH "Thrust directory")
  message(STATUS "Using THRUST_DIR=${THRUST_DIR}")

  # Allow for swapping backends
  set(THRUST_BACKEND "THRUST_DEVICE_SYSTEM_OMP" CACHE STRING "Thrust backend")
  message(STATUS "Using ${THRUST_BACKEND} as Thrust backend.")
  target_include_directories(wonton SYSTEM INTERFACE ${THRUST_DIR})
  target_compile_definitions(wonton INTERFACE THRUST_DEVICE_SYSTEM=${THRUST_BACKEND})

  set(WONTON_ENABLE_THRUST True CACHE BOOL "Is the Thrust library being used?")

  if("${THRUST_BACKEND}" STREQUAL "THRUST_DEVICE_SYSTEM_OMP")
    FIND_PACKAGE( OpenMP REQUIRED)
    if(OPENMP_FOUND)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif(OPENMP_FOUND)
  endif ()

  if("${THRUST_BACKEND}" STREQUAL "THRUST_DEVICE_SYSTEM_TBB")
    FIND_PACKAGE(TBB REQUIRED)
    if(TBB_FOUND)
      target_include_directories(wonton SYSTEM INTERFACE ${TBB_INCLUDE_DIRS})
      target_link_libraries(wonton SYSTEM INTERFACE ${TBB_LIBRARIES})
      target_compile_definitions(wonton INTERFACE ${TBB_DEFINITIONS})
    endif(TBB_FOUND)
  endif()

endif(ENABLE_THRUST)

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
# WONTON_ENABLE_MPI etc. if enabled. We write to a temporary name
# (with a .gen suffix) so that it is not seen during the build process
# and then rename it properly during the install step

configure_file(${PROJECT_SOURCE_DIR}/config/wonton-config.h.in
  ${PROJECT_BINARY_DIR}/wonton-config.h.gen @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/wonton-config.h.gen
  DESTINATION include RENAME wonton-config.h)


# Install the FindTHRUST module needed by downstream packages when
# processing wontonConfig.cmake
install(FILES ${PROJECT_SOURCE_DIR}/cmake/FindTHRUST.cmake
  DESTINATION lib/cmake/wonton/modules)
