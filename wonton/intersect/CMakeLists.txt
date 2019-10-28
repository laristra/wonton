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
#    PARENT_SCOPE
)

add_library(r3d_intersect ${intersect_SOURCES})

target_include_directories(support PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(support PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(support PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIRECTORY}>
  $<INSTALL_INTERFACE:include>)  #<cmake_intall_prefix>/include


target_link_libraries(wonton PRIVATE intersect)



# Add unit test
#if(ENABLE_UNIT_TESTS)
#endif(ENABLE_UNIT_TESTS)

