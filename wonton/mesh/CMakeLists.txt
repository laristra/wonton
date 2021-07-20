#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#
project(mesh)

add_library(wonton_mesh INTERFACE)

target_include_directories(wonton_mesh INTERFACE
  $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
  $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

set(wonton_mesh_HEADERS AuxMeshTopology.h)

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(wonton_mesh PROPERTIES
#   PUBLIC_HEADER "${wonton_mesh_HEADERS}")
#
# Directly install files instead
install(FILES ${wonton_mesh_HEADERS} DESTINATION include/wonton/mesh)

target_link_libraries(wonton_mesh INTERFACE wonton_support)

# Recurse down in subdirectories
if (WONTON_ENABLE_Jali)
  add_subdirectory(jali)
endif ()

if (WONTON_ENABLE_FleCSI)
  add_subdirectory(flecsi)
endif ()

if (WONTON_ENABLE_Conduit)
  add_subdirectory(conduit)
endif ()

add_subdirectory(simple)

add_subdirectory(adaptive_refinement)

add_subdirectory(direct_product)

add_subdirectory(flat)

# Make this a dependency of the main target
target_link_libraries(wonton INTERFACE wonton_mesh)

install(TARGETS wonton_mesh
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh
  INCLUDES DESTINATION include/wonton/mesh
  )

if (ENABLE_UNIT_TESTS)
  wonton_add_unittest(test_mesh_api
    SOURCES test/test_mesh_api.cc
    LIBRARIES wonton_mesh
    POLICY SERIAL
    )
endif ()
