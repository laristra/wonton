#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

add_library(adaptive_refinement_mesh_wrapper INTERFACE)

target_include_directories(adaptive_refinement_mesh_wrapper INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

set(armesh_HEADERS
  adaptive_refinement_mesh.h
  adaptive_refinement_mesh_wrapper.h
  )

# Not yet allowed for INTERFACE libraries
#
#set_target_properties(adaptive_refinement_mesh_wrapper PROPERTIES
#  PUBLIC_HEADER "${wonton_ar_HEADERS}")
#
# Directly install files instead
install(FILES ${armesh_HEADERS}
  DESTINATION include/wonton/mesh/adaptive_refinement)

target_link_libraries(adaptive_refinement_mesh_wrapper INTERFACE wonton_support)

target_link_libraries(wonton_mesh INTERFACE adaptive_refinement_mesh_wrapper)

install(TARGETS adaptive_refinement_mesh_wrapper
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh/adaptive_refinement
  INCLUDES DESTINATION include/wonton/mesh/adaptive_refinement
  )
      

if (ENABLE_UNIT_TESTS)
  
  # Adaptive Refinement Mesh
  wonton_add_unit(test_adaptive_refinement_mesh
    SOURCES test/test_adaptive_refinement_mesh.cc
    LIBRARIES adaptive_refinement_mesh_wrapper
    POLICY SERIAL
    )

  # Adaptive Refinement Mesh Wrapper
  wonton_add_unit(test_adaptive_refinement_mesh_wrapper
    SOURCES test/test_adaptive_refinement_mesh_wrapper.cc
    LIBRARIES adaptive_refinement_mesh_wrapper
    POLICY SERIAL
    )

endif (ENABLE_UNIT_TESTS)
