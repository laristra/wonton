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

if (ENABLE_UNIT_TESTING)

  # Direct Product Mesh
  cinch_add_unit(test_direct_product_mesh
    SOURCES test/test_direct_product_mesh.cc
    POLICY SERIAL
    )
  
  # Direct Product Mesh Wrapper
  cinch_add_unit(test_direct_product_mesh_wrapper
    SOURCES test/test_direct_product_mesh_wrapper.cc
    POLICY SERIAL
    )

endif (ENABLE_UNIT_TESTING)
