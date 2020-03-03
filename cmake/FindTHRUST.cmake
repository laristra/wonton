#------------------------------------------------------------------------------#
# Copyright (c) 2020 Triad National Security, LLC
# All rights reserved.
#
# Find the Thrust library
# This module will set the following variables
#
# THRUST_ROOT         - Where Thrust is installed/found
# THRUST_INCLUDE_DIRS - Where Thrust includes are
# THRUST_LIBRARIES    - Name of Thrust library target (INTERFACE target)
# THRUST_COMPONENTS   - Which backends were found (OpenMP, TBB, CUDA)
#------------------------------------------------------------------------------#

# Create an INTERFACE library

set(THRUST_LIBRARIES THRUST::THRUST)
add_library(${THRUST_LIBRARIES} INTERFACE IMPORTED)

#------------------------------------------------------------------------------#
# Find the header file
#------------------------------------------------------------------------------#

if (THRUST_ROOT)
  find_path(THRUST_INCLUDE_DIRS thrust/version.h PATHS ${THRUST_ROOT})
else ()
  find_path(THRUST_INCLUDE_DIRS thrust/version.h)
  if (THRUST_INCLUDE_DIRS)
    set(THRUST_ROOT ${THRUST_INCLUDE_DIRS}/.. CACHE FILEPATH "Where THRUST is found")
  endif ()
endif ()

target_include_directories(${THRUST_LIBRARIES} INTERFACE ${THRUST_INCLUDE_DIR})

#------------------------------------------------------------------------------#
# Find the backend components
#------------------------------------------------------------------------------#

set(THRUST_COMPONENTS "")
foreach (_component IN LISTS THRUST_FIND_COMPONENTS)
  if (${_component} STREQUAL "OpenMP")
    find_package(OpenMP)
    if (OpenMP_FOUND)
      target_link_libraries(${THRUST_LIBRARIES} INTERFACE OpenMP::OpenMP_CXX)
      list(APPEND THRUST_COMPONENTS OpenMP)
    endif ()
  elseif (${_component} STREQUAL "CUDA")
    find_package(CUDA)
    if (CUDA_FOUND)
      target_link_libraries(${THRUST_LIBRARIES} INTERFACE "${CUDA_LIBRARIES}")
      list(APPEND THRUST_COMPONENTS CUDA)
    endif ()
  elseif ()
    message(FATAL_ERROR "Unknown component ${_component) requested for THRUST")
  endif ()
endforeach ()


message(STATUS "THRUST components found: ${THRUST_COMPONENTS}")

#------------------------------------------------------------------------------#
# Set standard args stuff
#------------------------------------------------------------------------------#

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(THRUST
  REQUIRED_VARS THRUST_ROOT THRUST_LIBRARIES THRUST_INCLUDE_DIRS)


