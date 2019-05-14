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
  moment_index.h
  Point.h
  svd.h
  Vector.h
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

  cinch_add_unit(moment_index
    SOURCES test/test_moment_index.cc
    POLICY SERIAL
    )

endif(ENABLE_UNIT_TESTS)
