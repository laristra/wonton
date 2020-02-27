#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#
project(mesh)

add_library(wonton_mesh INTERFACE)

target_include_directories(wonton_mesh INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
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
if (ENABLE_Jali)
  add_subdirectory(jali)
endif (ENABLE_Jali)

if (FleCSI_DIR)
  add_subdirectory(flecsi)
endif (FleCSI_DIR)

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
  if (ENABLE_FleCSI)
    cinch_add_unit(test_flecsi_mesh_wrapper
      SOURCES flecsi/test/test_flecsi_mesh_wrapper.cc
      LIBRARIES wonton_mesh
      POLICY MPI
      THREADS 1)
  endif (ENABLE_FleCSI)
endif(ENABLE_UNIT_TESTS)
