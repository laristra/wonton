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

target_sources(mesh INTERFACE ${adaptive_refinement_mesh_wrapper_HEADERS})

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
