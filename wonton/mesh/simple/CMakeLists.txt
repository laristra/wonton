#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

set(simple_mesh_wapper_HEADERS
  simple_mesh_wrapper.h
  simple_mesh.h
  )

add_library(simple_mesh_wrapper INTERFACE)

target_include_directories(simple_mesh_wrapper INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

# add this header-only library as a dependency of upper level mesh library

target_link_libraries(wonton_mesh INTERFACE simple_mesh_wrapper)

install(TARGETS simple_mesh_wrapper
  EXPORT wonton-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh/simple
  INCLUDES DESTINATION include/wonton/mesh/simple
  )
      
