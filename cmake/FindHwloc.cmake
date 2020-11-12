#------------------------------------------------------------------------------#
# Copyright (c) 2020 Triad National Security, LLC
# All rights reserved.
#
# Find the HWLOC library
# This module will set the following variables
#
# HWLOC_ROOT         - Where HWLOC is installed/found
# HWLOC_INCLUDE_DIRS - Where HWLOC includes are
# HWLOC_LIBRARIES    - Name of HWLOC library target(s)
#------------------------------------------------------------------------------#

# First use pkg-config to parse an installed .pc file to find the
# library although we cannot rely on it

find_package(PkgConfig)
pkg_check_modules(PC_HWLOC QUIET hwloc)


#------------------------------------------------------------------------------#
# Find the header file
#------------------------------------------------------------------------------#

find_path(HWLOC_INCLUDE_DIR
  NAMES hwloc.h
  HINTS ${PC_HWLOC_INCLUDE_DIRS}
  PATHS ${HWLOC_ROOT}
  )
if (NOT HWLOC_INCLUDE_DIR)
  if (HWLOC_FIND_REQUIRED)
    message(FATAL_ERROR "Cannot find hwloc.h")
  else ()
    if (NOT HWLOC_FIND_QUIETLY)
      message(WARNING "Cannot find hwloc.h")
    endif ()
  endif ()
endif ()

if (HWLOC_INCLUDE_DIR)
  set(HWLOC_INCLUDE_DIRS ${HWLOC_INCLUDE_DIR} CACHE FILEPATH "Where hwloc.h is located")
  if (NOT HWLOC_ROOT)
    get_filename_component(hwloc_root "${HWLOC_INCLUDE_DIR}/.." ABSOLUTE)
    set(HWLOC_ROOT ${hwloc_root} CACHE FILEPATH "Where hwloc is found")
  endif ()
endif ()

#------------------------------------------------------------------------------#
# Find the library
#------------------------------------------------------------------------------#

find_library(HWLOC_LIBRARY
  NAMES hwloc
  HINTS ${PC_HWLOC_LIBRARY_DIRS}
  PATHS ${HWLOC_ROOT})

if (NOT HWLOC_LIBRARY)
  if (HWLOC_FIND_REQUIRED)
    message(FATAL_ERROR "Cannot locate HWLOC library")
  else ()
    if (NOT HWLOC_FIND_QUIETLY)
      message(WARNING "Cannot locate HWLOC library")
    endif ()
  endif ()
endif ()

set(HWLOC_VERSION PC_HWLOC_VERSION)  # No guarantee

# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HWLOC
  DEFAULT_MSG
  HWLOC_INCLUDE_DIRS
  HWLOC_LIBRARY)

# Create HWLOC target

if (HWLOC_FOUND AND NOT TARGET hwloc::hwloc)
  set(HWLOC_LIBRARIES hwloc:hwloc)
  add_library(${HWLOC_LIBRARIES} UNKNOWN IMPORTED)
  set_target_properties(${HWLOC_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${HWLOC_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${HWLOC_INCLUDE_DIR}")
endif ()


mark_as_advanced(HWLOC_INCLUDE_DIR HWLOC_LIBRARY)
