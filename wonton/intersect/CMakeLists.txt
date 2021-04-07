#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#
project(intersect)

add_subdirectory(r3d)
set(r3d_ROOT ${CMAKE_INSTALL_PREFIX} CACHE PATH "Path to R3D Installation")

add_library(wonton_r3d_intersect INTERFACE)

target_link_libraries(wonton_r3d_intersect INTERFACE r3d::r3d)

target_include_directories(wonton_r3d_intersect INTERFACE
  $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
  $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

target_link_libraries(wonton INTERFACE wonton_r3d_intersect)

install(TARGETS wonton_r3d_intersect
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  )


