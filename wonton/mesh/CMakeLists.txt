#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

if (Jali_DIR)
  set(jali_mesh_wrapper_headers
    jali/jali_mesh_wrapper.h
    )
endif (Jali_DIR)

if (ENABLE_FleCSI)
  set(flecsi_mesh_wrapper_headers
    flecsi/flecsi_mesh_wrapper.h
    )
endif ()

set(mesh_HEADERS
  AuxMeshTopology.h
  adaptive_refinement/adaptive_refinement_mesh.h
  adaptive_refinement/adaptive_refinement_mesh_wrapper.h
  direct_product/direct_product_mesh.h
  direct_product/direct_product_mesh_wrapper.h
  simple/simple_mesh.h
  simple/simple_mesh_wrapper.h
  flat/flat_mesh_wrapper.h
  ${jali_mesh_wrapper_headers}
  ${flecsi_mesh_wrapper_headers}
  PARENT_SCOPE
  )

if (ENABLE_UNIT_TESTS)

  # Adaptive Refinement Mesh
  cinch_add_unit(test_adaptive_refinement_mesh
    SOURCES adaptive_refinement/test/test_adaptive_refinement_mesh.cc
    POLICY SERIAL
    )

  # Adaptive Refinement Mesh Wrapper
  cinch_add_unit(test_adaptive_refinement_mesh_wrapper
    SOURCES adaptive_refinement/test/test_adaptive_refinement_mesh_wrapper.cc
    POLICY SERIAL
    )

  # Direct Product Mesh
  cinch_add_unit(test_direct_product_mesh
    SOURCES direct_product/test/test_direct_product_mesh.cc
    POLICY SERIAL
    )

  # Direct Product Mesh Wrapper
  cinch_add_unit(test_direct_product_mesh_wrapper
    SOURCES direct_product/test/test_direct_product_mesh_wrapper.cc
    POLICY SERIAL
    )

  # Simple Mesh
  cinch_add_unit(test_simple_mesh
    SOURCES simple/test/test_simple_mesh.cc
    POLICY SERIAL
    )

  # Simple Mesh wrapper
  cinch_add_unit(test_simple_mesh_wrapper
    SOURCES simple/test/test_simple_mesh_wrapper.cc
    POLICY SERIAL)

  if (Jali_DIR)
    cinch_add_unit(test_jali_mesh_wrapper
      SOURCES jali/test/test_jali_mesh_wrapper.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      # this is relative to the wrapper directory
      INPUTS ../../test_data/shotshell.exo
      ../../test_data/shotshell-v.exo
      POLICY MPI
      THREADS 1)

    # Flat mesh wrapper (has to be paired with a real MPI capable mesh
    # wrapper)
    cinch_add_unit(test_flat_mesh_wrapper
      SOURCES flat/test/test_flat_mesh_wrapper.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      # this is relative to the wrapper directory
      INPUTS ../../test_data/shotshell.exo
      ../../test_data/shotshell-v.exo
      POLICY MPI
      THREADS 1)
  endif (Jali_DIR)

  if (ENABLE_FleCSI)
    cinch_add_unit(test_flecsi_mesh_wrapper
      SOURCES flecsi/test/test_flecsi_mesh_wrapper.cc
      LIBRARIES ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES}
      POLICY MPI
      THREADS 1)
  endif (ENABLE_FleCSI)
endif(ENABLE_UNIT_TESTS)
