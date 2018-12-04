#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

# Add header files
set(intersect_HEADERS
    r3d/r2d.h
    r3d/r3d.h
    PARENT_SCOPE
)

# Add source files
set(intersect_SOURCES
    r3d/r2d.c
    r3d/r3d.c
    PARENT_SCOPE
)

# Add unit test
#if(ENABLE_UNIT_TESTS)
#endif(ENABLE_UNIT_TESTS)

