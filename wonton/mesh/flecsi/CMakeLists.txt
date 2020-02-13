#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

if (ENABLE_FleCSI)

  add_library(flecsi_mesh_wrapper INTERFACE)
  
  target_include_directories(flecsi_mesh_wrapper INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>)
  
  target_link_libraries(flecsi_mesh_wrapper INTERFACE wonton_support)

  # FleCSI dependencies
  target_include_directories(flecsi_mesh_wrapper INTERFACE
    ${FleCSI_INCLUDE_DIRS} ${FleCSISP_INCLUDE_DIRS})

  target_link_libraries(flecsi_mesh_wrapper INTERFACE
    ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES})

  # add this header-only library as a dependency of upper level mesh library
  
  target_link_libraries(wonton_mesh INTERFACE flecsi_mesh_wrapper)
  
  install(TARGETS flecsi_mesh_wrapper
    EXPORT wonton-targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/wonton/mesh/flecsi
    INCLUDES DESTINATION include/wonton/mesh/flecsi
    )


  cinch_add_unit(test_flecsi_mesh_wrapper
    SOURCES flecsi/test/test_flecsi_mesh_wrapper.cc
    LIBRARIES ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES}
    POLICY MPI
    THREADS 1)

  # unit test needs FleCSI includes
  target_include_directories(test_flecsi_mesh_wrapper
    PRIVATE ${FleCSI_INCLUDE_DIRS} ${FleCSISP_INCLUDE_DIRS})
  
endif ()
