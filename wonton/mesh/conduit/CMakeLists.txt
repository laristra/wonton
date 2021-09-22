#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

if (WONTON_ENABLE_Conduit)

  add_library(conduit_mesh_wrapper INTERFACE)

  target_include_directories(conduit_mesh_wrapper INTERFACE
    $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
    $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>)

  set(conduit_mesh_wrapper_HEADERS
    conduit_mesh_wrapper.h
    )

  # Not yet allowed for INTERFACE libraries
  #
  # set_target_properties(conduit_mesh_wrapper install(DIRECTORY ${wonton_SOURCE_DIR}/include/SI DESTINATION include)PROPERTIES
  #   PUBLIC_HEADER "${conduit_mesh_wrapper_HEADERS}")
  #
  # Directly install files instead
  install(FILES ${conduit_mesh_wrapper_HEADERS}
    DESTINATION include/wonton/mesh/conduit)
  
  target_link_libraries(conduit_mesh_wrapper INTERFACE wonton_support)


  # Conduit dependencies
  target_link_libraries(conduit_mesh_wrapper INTERFACE ${Conduit_LIBRARIES})

  target_link_libraries(wonton_mesh INTERFACE conduit_mesh_wrapper)

  install(TARGETS conduit_mesh_wrapper
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh/conduit
  INCLUDES DESTINATION include/wonton/mesh/conduit
  )

  if (ENABLE_UNIT_TESTS)
    wonton_add_unittest(test_conduit_mesh_wrapper
      SOURCES test/test_conduit_mesh_wrapper.cc
      LIBRARIES conduit_mesh_wrapper
      POLICY SERIAL)  
  endif ()

endif()
