#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#
project(mesh)


# Only header files in this directory and it's subdirectories
set(mesh_HEADERS
  AuxMeshTopology.h
  )


add_library(wonton_mesh INTERFACE)

target_include_directories(wonton_mesh INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)


if (ENABLE_Jali)
  add_subdirectory(jali)
endif (ENABLE_Jali)

if (FleCSI_DIR)
  add_subdirectory(flecsi)
endif (FleCSI_DIR)

add_subdirectory(simple)

add_subdirectory(adaptive_refinement)

add_subdirectory(direct_product)

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
