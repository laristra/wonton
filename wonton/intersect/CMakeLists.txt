#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#
project(intersect)

add_library(wonton_r3d_intersect)
target_sources(wonton_r3d_intersect PRIVATE r3d/r3d.c r3d/r2d.c)

target_include_directories(wonton_r3d_intersect PUBLIC
  $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
  $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

set(R3D_HEADERS r3d/r3d.h r3d/r2d.h)
set_target_properties(wonton_r3d_intersect PROPERTIES
  PUBLIC_HEADER "${R3D_HEADERS}"
  )

target_link_libraries(wonton INTERFACE wonton_r3d_intersect)

install(TARGETS wonton_r3d_intersect
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/intersect/r3d
  INCLUDES DESTINATION include/wonton/intersect/r3d
  )


