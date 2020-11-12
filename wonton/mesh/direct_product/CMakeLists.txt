#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

add_library(direct_product_mesh_wrapper INTERFACE)

target_include_directories(direct_product_mesh_wrapper INTERFACE
  $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
  $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

set(dpmesh_wrapper_HEADERS
  direct_product_mesh.h
  direct_product_mesh_wrapper.h
  )

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(direct_product_mesh_wrapper PROPERTIES
#  PUBLIC_HEADER "${dpmesh_wrapper_HEADERS}")
#
# Directly install files instead
install(FILES ${dpmesh_wrapper_HEADERS}
  DESTINATION include/wonton/mesh/direct_product)

target_link_libraries(direct_product_mesh_wrapper INTERFACE wonton_support)

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
  wonton_add_unittest(test_direct_product_mesh
    SOURCES test/test_direct_product_mesh.cc
    LIBRARIES direct_product_mesh_wrapper
    POLICY SERIAL
    )
  
  # Direct Product Mesh Wrapper
  wonton_add_unittest(test_direct_product_mesh_wrapper
    SOURCES test/test_direct_product_mesh_wrapper.cc
    LIBRARIES direct_product_mesh_wrapper
    POLICY SERIAL
    )

endif (ENABLE_UNIT_TESTING)
