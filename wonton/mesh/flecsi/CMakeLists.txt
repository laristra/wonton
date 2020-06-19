#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

if (WONTON_ENABLE_FleCSI)

  add_library(flecsi_mesh_wrapper INTERFACE)
  
  target_include_directories(flecsi_mesh_wrapper INTERFACE
    $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
    $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>)

  set(flecsi_mesh_wrapper_HEADERS
    flecsi_mesh_wrapper.h
    )

  # Not yet allowed for INTERFACE libraries
  #
  #  set_target_properties(flecsi_mesh_wrapper PROPERTIES
  #    PUBLIC_HEADER "${flecsi_mesh_wrapper_HEADERS}")
  #
  # Directly install files instead
  install(FILES ${flecsi_mesh_wrapper_HEADERS}
    DESTINATION include/wonton/mesh/flecsi)
  
  target_link_libraries(flecsi_mesh_wrapper INTERFACE wonton_support)

  # FleCSI dependencies
  target_include_directories(flecsi_mesh_wrapper INTERFACE
    ${FleCSI_INCLUDE_DIR} ${FleCSISP_INCLUDE_DIR})

  target_link_libraries(flecsi_mesh_wrapper INTERFACE
    ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES})

  # add this header-only library as a dependency of upper level mesh library
  
  target_link_libraries(wonton_mesh INTERFACE flecsi_mesh_wrapper)
  
  install(TARGETS flecsi_mesh_wrapper
    EXPORT wonton_LIBRARIES
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/wonton/mesh/flecsi
    INCLUDES DESTINATION include/wonton/mesh/flecsi
    )


  wonton_add_unittest(test_flecsi_mesh_wrapper
    SOURCES test/test_flecsi_mesh_wrapper.cc
    LIBRARIES flecsi_mesh_wrapper
    POLICY MPI
    THREADS 1)

endif ()
