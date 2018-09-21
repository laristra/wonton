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
  simple/simple_mesh.h
  simple/simple_mesh_wrapper.h
  flat/flat_mesh_wrapper.h
  ${jali_mesh_wrapper_headers}
  ${flecsi_mesh_wrapper_headers}
  PARENT_SCOPE
  )

if (ENABLE_UNIT_TESTS)

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
    cinch_add_unit(flecsi_2d
      SOURCES flecsi/test/flecsi_2d.cc
      LIBRARIES ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES}
      POLICY MPI
      THREADS 1)
    
    cinch_add_unit(flecsi_3d
      SOURCES flecsi/test/flecsi_3d.cc
      LIBRARIES ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES}
      POLICY MPI
      THREADS 1)
    
    #cinch_add_unit(test_flecsi_mesh_wrapper
    #               SOURCES mesh/flecsi/test/test_flecsi_mesh_wrapper.cc
    #               LIBRARIES wonton ${FLECSI_LIBRARIES})
  endif ()
  
endif(ENABLE_UNIT_TESTS)

