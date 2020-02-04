#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#
project(support)

# Add header files
set(support_HEADERS
  lsfits.h
  CoordinateSystem.h
  Matrix.h
  moment_index.h
  Point.h
  svd.h
  Vector.h
  Polytope.h
  prime_factors.h
  equifactor.h
  structured_partitioner.h
  wonton.h
)

# Add source files
set(support_SOURCES
  Matrix.cc
  svd.cc
)

add_library(wonton_support ${support_SOURCES})

target_include_directories(wonton_support PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
  $<INSTALL_INTERFACE:include/wonton/support>)

set_target_properties(wonton_support PROPERTIES PUBLIC_HEADER ${support_HEADERS})


#-----------------------------------------------------------------------------
# Thrust information
#-----------------------------------------------------------------------------
set(ENABLE_THRUST FALSE CACHE BOOL "Use Thrust")  # create a cache var with default
if (ENABLE_THRUST)   # if it is overridden by the command line

  message(STATUS "Enabling compilation with Thrust")
  message(STATUS "Using THRUST_ROOT=${THRUST_ROOT}")

  # Allow for swapping backends
  set(THRUST_HOST_BACKEND "THRUST_HOST_SYSTEM_CPP"  CACHE STRING "Thrust host backend")
  set(THRUST_DEVICE_BACKEND "THRUST_DEVICE_SYSTEM_OMP" CACHE STRING "Thrust device backend")
  message(STATUS "Using ${THRUST_HOST_BACKEND} as Thrust host backend and ${THRUST_DEVICE_BACKEND} as Thrust device backend.")

  target_include_directories(wonton_support PUBLIC ${THRUST_ROOT})
  target_compile_definitions(wonton_support PUBLIC THRUST)
  target_compile_definitions(wonton_support PUBLIC THRUST_HOST_SYSTEM=${THRUST_HOST_BACKEND})
  target_compile_definitions(wonton_support PUBLIC THRUST_DEVICE_SYSTEM=${THRUST_DEVICE_BACKEND})

  set(WONTON_ENABLE_THRUST True CACHE BOOL "Is the Thrust library being used?" FORCE)

  if ("${THRUST_DEVICE_BACKEND}" STREQUAL "THRUST_DEVICE_SYSTEM_OMP")
    find_package(OpenMP REQUIRED)
    target_link_libraries(wonton_support PUBLIC ${OpenMP_LIBRARIES})
  endif ()

endif()

if (ENABLE_TCMALLOC)
  if (NOT TCMALLOC_LIB)
    message(WARNING "TCMALLOC enabled but TCMALLOC_LIB not specified")
  else ()
    target_link_libraries(wonton_support PUBLIC ${TCMALLOC_LIB})
  endif ()
endif ()



target_link_libraries(wonton INTERFACE wonton_support)

install(TARGETS wonton_support
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_point
    SOURCES test/test_point.cc
    POLICY SERIAL)

  cinch_add_unit(test_vector
      SOURCES test/test_vector.cc
      POLICY SERIAL)

  cinch_add_unit(test_matfuncs
    SOURCES test/test_matfuncs.cc
    LIBRARIES wonton_support
    POLICY MPI
    THREADS 1)

  cinch_add_unit(test_coord_sys
    SOURCES test/test_coord_sys.cc
    LIBRARIES wonton_support
    POLICY SERIAL
    )

  cinch_add_unit(test_moment_index
    SOURCES test/test_moment_index.cc
    LIBRARIES wonton_support
    POLICY SERIAL
    )

  cinch_add_unit(test_Polytope2
    SOURCES test/test_Polytope2.cc
    LIBRARIES wonton_support
    POLICY SERIAL)

  cinch_add_unit(test_Polytope3
    SOURCES test/test_Polytope3.cc
    LIBRARIES wonton_support
    POLICY SERIAL) 

  cinch_add_unit(test_prime_factors
    SOURCES test/test_prime_factors.cc
    POLICY SERIAL)

  cinch_add_unit(test_equifactor
    SOURCES test/test_equifactor.cc
    POLICY SERIAL)
  
  cinch_add_unit(test_structured_partitioner
    SOURCES test/test_structured_partitioner.cc
    POLICY SERIAL)

  cinch_add_unit(test_structured_partitioner_par
    SOURCES test/test_structured_partitioner_par.cc
    POLICY MPI
    THREADS 4)
  
endif(ENABLE_UNIT_TESTS)
