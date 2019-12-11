#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

# Add header files
set(support_HEADERS
  lsfits.h
  CoordinateSystem.h
  Matrix.h
  Point.h
  svd.h
  Vector.h
  Polytope.h
  prime_factors.h
  equifactor.h
  structured_partitioner.h
  wonton.h
  PARENT_SCOPE
)

# Add source files
set(support_SOURCES
  Matrix.cc
  svd.cc
  PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_matfuncs
    SOURCES test/test_matfuncs.cc
    LIBRARIES wonton
    POLICY MPI
    THREADS 1)

  cinch_add_unit(test_Polytope2
    SOURCES test/test_Polytope2.cc
    POLICY SERIAL)

  cinch_add_unit(test_Polytope3
    SOURCES test/test_Polytope3.cc
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
