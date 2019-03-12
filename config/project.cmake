#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

project(wonton CXX)

# SEMANTIC VERSION NUMBERS - UPDATE DILIGENTLY
# As soon as a change with a new version number is merged into the master,
# tag the central repository.

set(WONTON_VERSION_MAJOR 0)
set(WONTON_VERSION_MINOR 0)
set(WONTON_VERSION_PATCH b2fa9daa213)



cinch_minimum_required(VERSION 1.0)

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

# set the name of the WONTON library

set(WONTON_LIBRARY "wonton" CACHE STRING "Name of the wonton library")

#-----------------------------------------------------------------------------
# Gather all the third party libraries needed for Wonton
#-----------------------------------------------------------------------------

set(WONTON_EXTRA_LIBRARIES)

#------------------------------------------------------------------------------#
# Set up MPI builds
# (eventually most of this should be pushed down into cinch)
#------------------------------------------------------------------------------#
set(ENABLE_MPI OFF CACHE BOOL "")
if (ENABLE_MPI)
  find_package(MPI REQUIRED)
  set(WONTON_ENABLE_MPI True CACHE BOOL "Whether MPI is enabled")
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE FILEPATH "C compiler to use" FORCE)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "C compiler to use" FORCE)
endif (ENABLE_MPI)


set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

#-----------------------------------------------------------------------------
# FleCSI and FleCSI-SP location
#-----------------------------------------------------------------------------

set(ENABLE_FleCSI FALSE CACHE BOOL "Use FleCSI")
if (ENABLE_FleCSI AND NOT FleCSI_LIBRARIES)

 find_package(FleCSI REQUIRED)
 message(STATUS "FleCSI_LIBRARIES=${FleCSI_LIBRARIES}" )
 include_directories(${FleCSI_INCLUDE_DIR})
 message(STATUS "FleCSI_INCLUDE_DIRS=${FleCSI_INCLUDE_DIR}")
 list(APPEND WONTON_EXTRA_LIBRARIES ${FleCSI_LIBRARIES})

 find_package(FleCSISP REQUIRED)
 message(STATUS "FleCSISP_LIBRARIES=${FleCSISP_LIBRARIES}" )
 include_directories(${FleCSISP_INCLUDE_DIR})
 message(STATUS "FleCSISP_INCLUDE_DIRS=${FleCSISP_INCLUDE_DIR}")
 list(APPEND WONTON_EXTRA_LIBRARIES ${FleCSISP_LIBRARIES})

endif(ENABLE_FleCSI AND NOT FleCSI_LIBRARIES)


#------------------------------------------------------------------------------#
# Configure Jali
# (this includes the TPLs that Jali will need)
#------------------------------------------------------------------------------#

if (Jali_DIR AND NOT Jali_LIBRARIES)

   # Look for the Jali package

   find_package(Jali REQUIRED
                HINTS ${Jali_DIR}/lib)

   message(STATUS "Located Jali")
   message(STATUS "Jali_DIR=${Jali_DIR}")

   # add full path to jali libs
   unset(_LIBS)
   foreach (_lib ${Jali_LIBRARIES})
      set(_LIBS ${_LIBS};${Jali_LIBRARY_DIRS}/lib${_lib}.a)
   endforeach()
   set(Jali_LIBRARIES ${_LIBS})

   include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})

   list(APPEND WONTON_EXTRA_LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})
endif (Jali_DIR AND NOT Jali_LIBRARIES)

#-----------------------------------------------------------------------------
# Include Boost libraries for builds without Jali
#-----------------------------------------------------------------------------

if ( NOT Jali_DIR)

 find_package(Boost REQUIRED)
 include_directories(${Boost_INCLUDE_DIR})
 message(STATUS "Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIR}")

endif( NOT Jali_DIR)

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
  include_directories(${LAPACKE_INCLUDE_DIRS})
  add_definitions("-DHAVE_LAPACKE")

  list(APPEND WONTON_EXTRA_LIBRARIES ${LAPACKE_LIBRARIES})

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
  include_directories(${THRUST_DIR})
  add_definitions(-DTHRUST_DEVICE_SYSTEM=${THRUST_BACKEND})

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
      include_directories(${TBB_INCLUDE_DIRS})
      link_directories(${TBB_LIBRARY_DIRS})
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltbb")
    endif(TBB_FOUND)
  endif()

endif(ENABLE_THRUST)

#-----------------------------------------------------------------------------
# Now add the source directories
#-----------------------------------------------------------------------------

# In addition to the include directories of the source set by cinch,
# we need to include the build directory to get the autogenerated
# wonton-config.h

include_directories(${CMAKE_BINARY_DIRECTORY})

# Libraries

cinch_add_library_target(wonton wonton)
# TODO - merge LAPACKE_LIBRARIES into WONTON_LIBRARIES
cinch_target_link_libraries(wonton ${WONTON_LIBRARIES})
cinch_target_link_libraries(wonton ${LAPACKE_LIBRARIES})


# build the WONTON_LIBRARIES variable
set(WONTON_LIBRARIES ${WONTON_LIBRARY} ${WONTON_EXTRA_LIBRARIES} CACHE STRING "List of libraries to link with wonton")

# retrieve all the definitions we added for compiling
get_directory_property(WONTON_COMPILE_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)

############################################################################## 
# Write a configuration file from template replacing only variables enclosed
# by the @ sign. This will let other programs build on WONTON discover how
# WONTON was built and which TPLs it used
#############################################################################

configure_file(${PROJECT_SOURCE_DIR}/cmake/wonton_config.cmake.in 
  ${PROJECT_BINARY_DIR}/wonton_config.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/wonton_config.cmake 
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/cmake/)


configure_file(${PROJECT_SOURCE_DIR}/config/wonton-config.h.in
  ${PROJECT_BINARY_DIR}/wonton-config.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/wonton-config.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
