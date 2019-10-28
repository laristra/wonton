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

add_library(simple_mesh_wrapper INTERFACE ${simple_mesh_wrapper_HEADERS})
target_include_directories(simple_mesh_wrapper INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}>)

# add this header-only library as a dependency of upper level mesh library

target_link_libraries(mesh INTERFACE simple_mesh_wrapper)
