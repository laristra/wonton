#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

project(wonton CXX)

cinch_minimum_required(VERSION 1.0)

cmake_minimum_required(VERSION 3.13)

if (CMAKE_VERSION_MAJOR GREATER_EQUAL 3.13)
  CMAKE_POLICY(SET CMP0079 NEW)  # allow target_link_libraries to reference
                                 # targets from other directories
endif()
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
  cmake_policy(SET CMP0074 NEW)  # Don't ignore Pkg_ROOT variables
endif()



# SEMANTIC VERSION NUMBERS - UPDATE DILIGENTLY
# As soon as a change with a new version number is merged into the master,
# tag the central repository.

set(WONTON_VERSION_MAJOR 1)
set(WONTON_VERSION_MINOR 1)
set(WONTON_VERSION_PATCH 3)

# If a C++14 compiler is available, then set the appropriate flags
include(cxx14)
check_for_cxx14_compiler(CXX14_COMPILER)
if(CXX14_COMPILER)
  enable_cxx14()
else()
  message(STATUS "C++14 compatible compiler not found")
endif()

# If we couldn't find a C++14 compiler, try to see if a C++11
# compiler is available, then set the appropriate flags
if (NOT CXX14_COMPILER)
  include(cxx11)
  check_for_cxx11_compiler(CXX11_COMPILER)
  if(CXX11_COMPILER)
    enable_cxx11()
  else()
    message(FATAL_ERROR "C++11 compatible compiler not found")
  endif()
endif()


# cinch extras

cinch_load_extras()

set(CINCH_HEADER_SUFFIXES "\\.h")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")


# Top level target
add_library(wonton INTERFACE)

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(wonton::wonton ALIAS wonton)


#------------------------------------------------------------------------------#
# Set up MPI builds
# (eventually most of this should be pushed down into cinch)
#------------------------------------------------------------------------------#
set(ENABLE_MPI OFF CACHE BOOL "")
if (ENABLE_MPI)
  find_package(MPI REQUIRED)
  set(WONTON_ENABLE_MPI True CACHE BOOL "Whether MPI is enabled")
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE FILEPATH "C compiler to use" FORCE)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "C++ compiler to use" FORCE)
endif (ENABLE_MPI)

#-----------------------------------------------------------------------------
# FleCSI and FleCSI-SP location
#-----------------------------------------------------------------------------

set(ENABLE_FleCSI FALSE CACHE BOOL "Use FleCSI")
if (ENABLE_FleCSI AND NOT FleCSI_LIBRARIES)

  find_package(FleCSI REQUIRED)

  message(STATUS "FleCSI_LIBRARIES=${FleCSI_LIBRARIES}" )
  target_include_directories(wonton PUBLIC ${FleCSI_INCLUDE_DIR})

  message(STATUS "FleCSI_INCLUDE_DIRS=${FleCSI_INCLUDE_DIR}")
  target_link_libaries(wonton PUBLIC ${FleCSI_LIBRARIES})
  
  find_package(FleCSISP REQUIRED)

  message(STATUS "FleCSISP_LIBRARIES=${FleCSISP_LIBRARIES}" )
  target_include_directories(wonton PUBLIC ${FleCSISP_INCLUDE_DIR})

  message(STATUS "FleCSISP_INCLUDE_DIRS=${FleCSISP_INCLUDE_DIR}")
  target_link_libraries(wonton INTERFACE ${FleCSISP_LIBRARIES})

endif(ENABLE_FleCSI AND NOT FleCSI_LIBRARIES)


#------------------------------------------------------------------------------#
# Configure Jali
# (this includes the TPLs that Jali will need)
#------------------------------------------------------------------------------#

if (ENABLE_Jali AND NOT Jali_LIBRARIES)

   # Look for the Jali package

   find_package(Jali REQUIRED)  # specify in Jali_ROOT or CMAKE_PREFIX_PATH

   message(STATUS "Located Jali")
   message(STATUS "Jali_ROOT=${Jali_ROOT}")

   set(ENABLE_Jali ON)

   message(STATUS "Jali_LIBRARIES ${Jali_LIBRARIES}")
   target_link_libraries(wonton INTERFACE ${Jali_LIBRARIES})
  
endif (Jali_DIR AND NOT Jali_LIBRARIES)


#-----------------------------------------------------------------------------
# Include Boost libraries for builds without Jali
#-----------------------------------------------------------------------------

if (NOT ENABLE_Jali)
  find_package(Boost REQUIRED)
  target_include_directories(wonton SYSTEM PUBLIC ${Boost_INCLUDE_DIR})
  message(STATUS "Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIR}")

  target_link_libraries(wonton INTERFACE ${Boost_LIBRARIES})
endif( NOT Jali_DIR)


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

if (LAPACKE_DIR)
  # Directly look for cmake config file in LAPACKE_DIR and below
  file(GLOB_RECURSE LAPACKE_CONFIG_FILE ${LAPACKE_DIR}/lapacke-config.cmake)

  if (NOT LAPACKE_CONFIG_FILE)
    message(FATAL_ERROR " LAPACKE CMAKE config file not found under LAPACKE_DIR (${LAPACKE_DIR})")
  endif (NOT LAPACKE_CONFIG_FILE)

  message(STATUS "LAPACKE_CONFIG_FILE ${LAPACKE_CONFIG_FILE}")

  get_filename_component(LAPACKE_CONFIG_PATH ${LAPACKE_CONFIG_FILE} DIRECTORY)
  message(status " LAPACKE_CONFIG_PATH ${LAPACKE_CONFIG_PATH}")

  # If successful, the config file will set LAPACKE_LIBRARIES,
  # LAPACKE_lapack_LIBRARIES and LAPACKE_blas_LIBRARIES

  find_package(LAPACKE NO_MODULE NO_DEFAULT_PATH HINTS ${LAPACKE_CONFIG_PATH})

  if (LAPACKE_LIBRARIES STREQUAL "lapacke")

    # LAPACKE config file does not set the library path but it does set the
    # LAPACKE_INCLUDE_DIRS path. Try to back out the library path using this
    # and the top level directory as starting points for a find_library command

    find_library(LAPACKE_LIBRARY NAMES lapacke
                 NO_CMAKE_SYSTEM_PATH NO_DEFAULT_PATH
                 HINTS ${LAPACKE_DIR} ${LAPACKE_INCLUDE_DIRS}/..
	         PATH_SUFFIXES lib lib64)


    # Extract path of directory in which library files live to pass as a lib
    # search directory for the linker to find lapacke, lapack and blas libs

    get_filename_component(LAPACKE_LIBRARY_DIR ${LAPACKE_LIBRARY} DIRECTORY)

    set(LAPACKE_LIBRARIES "-Wl,-rpath,${LAPACKE_LIBRARY_DIR} -L${LAPACKE_LIBRARY_DIR} -l${LAPACKE_LIBRARIES} -l${LAPACK_lapack_LIBRARIES} -l${LAPACK_blas_LIBRARIES}")

    # If we don't want to link with Fortran then we have to tell it to link
    # with the Fortran libraries because LAPACK is written/compiled in Fortran
    #
    # NEEDED FOR STATIC LAPACK LIBS

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(LAPACKE_LIBRARIES "${LAPACKE_LIBRARIES} -lgfortran")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      set(LAPACKE_LIBRARIES "${LAPACKE_LIBRARIES} -lifcore")
    endif()

  endif(LAPACKE_LIBRARIES STREQUAL "lapacke")

else (LAPACKE_DIR)

  # Use FindLAPACKE.cmake provided by cinch or cmake to find it
  # FindLAPACKE.cmake provided by cinch requires PC_LAPACKE_INCLUDE_DIRS and
  # PC_LAPACKE_LIBRARY to be able to find LAPACKE

  find_package(LAPACKE)

endif (LAPACKE_DIR)

if (LAPACKE_FOUND)
  enable_language(Fortran)
  include(FortranCInterface)  # will ensure the fortran library is linked in
  
  target_include_directories(wonton INTERFACE ${LAPACKE_INCLUDE_DIRS})
  target_compile_definitions(wonton INTERFACE HAVE_LAPACKE)

  target_link_libraries(wonton INTERFACE ${LAPACKE_LIBRARIES})

  message(STATUS "LAPACKE_FOUND ${LAPACKE_FOUND}")
  message(STATUS "LAPACKE_LIBRARIES  ${LAPACKE_LIBRARIES}")
else (LAPACKE_FOUND)
   unset(LAPACKE_LIBRARIES)  # otherwise it will be LAPACKE-NOTFOUND or something
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
  target_include_directories(wonton SYSTEM PUBLIC ${THRUST_DIR})
  target_compile_definitions(wonton PUBLIC THRUST_DEVICE_SYSTEM=${THRUST_BACKEND})

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
      target_include_directories(wonton SYSTEM PUBLIC ${TBB_INCLUDE_DIRS})
      target_link_libraries(wonton PUBLIC ${TBB_LIBRARIES})
      target_compile_definitions(wonton PUBLIC ${TBB_DEFINITIONS})
    endif(TBB_FOUND)
  endif()

endif(ENABLE_THRUST)

#-----------------------------------------------------------------------------
# Recurse down the source directories building up dependencies
#-----------------------------------------------------------------------------

add_subdirectory(wonton)

# In addition to the include directories of the source, we need to
# include the build or directory to get the autogenerated
# wonton-config.h

target_include_directories(wonton INTERFACE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIRECTORY}>
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

# set the name of the WONTON library

set(WONTON_LIBRARY "wonton" CACHE STRING "Name of the wonton library")

# Write a configuration file from template replacing only variables enclosed
# by the @ sign.

configure_file(${PROJECT_SOURCE_DIR}/cmake/wontonConfig.cmake.in 
  wontonConfig.cmake @ONLY)
install(FILES wontonConfig.cmake DESTINATION lib/cmake/wonton)


# write out a version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(wontonConfigVersion.cmake
  VERSION "${WONTON_MAJOR_VERSION}.${WONTON_MINOR_VERSION}.${WONTON_PATCH_VERSION}"
  COMPATIBILITY SameMajorVersion)
install(FILES wontonConfigVersion.cmake DESTINATION lib/cmake/wonton)


# export targets

install(EXPORT wonton_LIBRARIES
  FILE wontonTargets.cmake
  NAMESPACE wonton::
  EXPORT_LINK_INTERFACE_LIBRARIES
  DESTINATION lib/cmake)



# Dynamically configured header files that contains defines like
# WONTON_ENABLE_MPI etc. if enabled

configure_file(${PROJECT_SOURCE_DIR}/config/wonton-config.h.in
  ${PROJECT_BINARY_DIR}/wonton-config.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/wonton-config.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
