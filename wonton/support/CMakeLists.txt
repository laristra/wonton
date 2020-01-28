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
    POLICY SERIAL
    )

  cinch_add_unit(test_moment_index
    SOURCES test/test_moment_index.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_polytope_2
    SOURCES test/test_polytope_2.cc
    POLICY SERIAL)

  cinch_add_unit(test_polytope_3
    SOURCES test/test_polytope_3.cc
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
