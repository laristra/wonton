#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

if (ENABLE_FleCSI)
  set(flecsi_mesh_wrapper_HEADERS
    flecsi_mesh_wrapper.h
    )
  
  target_sources(mesh INTERFACE ${flecsi_mesh_wrapper_HEADERS})

  # FleCSI dependencies
  target_include_directories(mesh INTERFACE ${FleCSI_INCLUDE_DIRS} ${FleCSISP_INCLUDE_DIRS})
  target_link_libraries(mesh INTERFACE ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES})

  cinch_add_unit(test_flecsi_mesh_wrapper
    SOURCES flecsi/test/test_flecsi_mesh_wrapper.cc
    LIBRARIES ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES}
    POLICY MPI
    THREADS 1)

  # unit test needs FleCSI includes
  target_include_directories(test_flecsi_mesh_wrapper
    PRIVATE ${FleCSI_INCLUDE_DIRS} ${FleCSISP_INCLUDE_DIRS})
  
endif ()
