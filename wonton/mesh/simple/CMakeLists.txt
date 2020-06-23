#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

add_library(simple_mesh_wrapper INTERFACE)

target_include_directories(simple_mesh_wrapper INTERFACE
  $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
  $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

set(simple_mesh_HEADERS
  simple_mesh.h
  simple_mesh_wrapper.h
  )

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(simple_mesh_wrapper PROPERTIES
#  PUBLIC_HEADER "${simple_mesh_wrapper_HEADERS}")
#
# Directly install files instead
install(FILES ${simple_mesh_HEADERS} DESTINATION include/wonton/mesh/simple)

target_link_libraries(simple_mesh_wrapper INTERFACE wonton_support)

# add this header-only library as a dependency of upper level mesh library

target_link_libraries(wonton_mesh INTERFACE simple_mesh_wrapper)

install(TARGETS simple_mesh_wrapper
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh/simple
  INCLUDES DESTINATION include/wonton/mesh/simple
  )
      
if (ENABLE_UNIT_TESTS)
  wonton_add_unittest(test_simple_mesh
    SOURCES test/test_simple_mesh.cc
    LIBRARIES simple_mesh_wrapper
    POLICY SERIAL
    )

  # Simple Mesh wrapper
  wonton_add_unittest(test_simple_mesh_wrapper
    SOURCES test/test_simple_mesh_wrapper.cc
    LIBRARIES simple_mesh_wrapper
    POLICY SERIAL)  
endif ()
