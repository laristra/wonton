#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

set(adaptive_refinement_mesh_wapper_HEADERS
  adaptive_refinement_mesh_wrapper.h
  adaptive_refinement_mesh.h
  )

add_library(adaptive_refinement_mesh_wrapper INTERFACE)

target_include_directories(adaptive_refinement_mesh_wrapper INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(wonton_mesh INTERFACE adaptive_refinement_mesh_wrapper)

install(TARGETS adaptive_refinement_mesh_wrapper
  EXPORT wonton-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh/adaptive_refinement
  INCLUDES DESTINATION include/wonton/mesh/adaptive_refinement
  )
      

if (ENABLE_UNIT_TESTS)
  
  # Adaptive Refinement Mesh
  cinch_add_unit(test_adaptive_refinement_mesh
    SOURCES test/test_adaptive_refinement_mesh.cc
    POLICY SERIAL
    )

  # Adaptive Refinement Mesh Wrapper
  cinch_add_unit(test_adaptive_refinement_mesh_wrapper
    SOURCES test/test_adaptive_refinement_mesh_wrapper.cc
    POLICY SERIAL
    )

endif (ENABLE_UNIT_TESTS)
