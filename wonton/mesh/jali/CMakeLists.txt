#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

if (WONTON_ENABLE_Jali)

  add_library(jali_mesh_wrapper INTERFACE)

  target_include_directories(jali_mesh_wrapper INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>)

  set(jali_mesh_wrapper_HEADERS
    jali_mesh_wrapper.h
    )

  # Not yet allowed for INTERFACE libraries
  #
  # set_target_properties(jali_mesh_wrapper install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/SI DESTINATION include)PROPERTIES
  #   PUBLIC_HEADER "${jali_mesh_wrapper_HEADERS}")
  #
  # Directly install files instead
  install(FILES ${jali_mesh_wrapper_HEADERS}
    DESTINATION include/wonton/mesh/jali)
  
  target_link_libraries(jali_mesh_wrapper INTERFACE wonton_support)


  # Jali dependencies
  target_include_directories(jali_mesh_wrapper INTERFACE ${Jali_INCLUDE_DIRS})
  target_link_libraries(jali_mesh_wrapper INTERFACE ${Jali_LIBRARIES})

  target_link_libraries(wonton_mesh INTERFACE jali_mesh_wrapper)

  install(TARGETS jali_mesh_wrapper
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh/jali
  INCLUDES DESTINATION include/wonton/mesh/jali
  )

  if (ENABLE_UNIT_TESTS)
    wonton_add_unit(test_jali_mesh_wrapper
      SOURCES test/test_jali_mesh_wrapper.cc
      LIBRARIES jali_mesh_wrapper
      # this is relative to the wrapper directory
      INPUTS ${wonton_SOURCE_DIR}/test_data/shotshell.exo
      ${wonton_SOURCE_DIR}/test_data/shotshell-v.exo
      POLICY MPI
      THREADS 1)
  endif ()

endif()
