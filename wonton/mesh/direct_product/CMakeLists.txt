#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

set(direct_product_mesh_wapper_HEADERS
  direct_product_mesh_wrapper.h
  direct_mesh.h
  )

add_library(direct_product_mesh_wrapper INTERFACE)

target_include_directories(direct_product_mesh_wrapper INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(wonton_mesh INTERFACE direct_product_mesh_wrapper)

install(TARGETS direct_product_mesh_wrapper
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh/direct_product
  INCLUDES DESTINATION include/wonton/mesh/direct_product
  )
      

if (ENABLE_UNIT_TESTING)

  # Direct Product Mesh
  cinch_add_unit(test_direct_product_mesh
    SOURCES test/test_direct_product_mesh.cc
    LIBRARIES direct_product_mesh_wrapper
    POLICY SERIAL
    )
  
  # Direct Product Mesh Wrapper
  cinch_add_unit(test_direct_product_mesh_wrapper
    SOURCES test/test_direct_product_mesh_wrapper.cc
    LIBRARIES direct_product_mesh_wrapper
    POLICY SERIAL
    )

endif (ENABLE_UNIT_TESTING)
