#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

if (ENABLE_Jali)
  set(jali_mesh_wrapper_HEADERS
    jali_mesh_wrapper.h
    )

  add_library(jali_mesh_wrapper INTERFACE)

  target_include_directories(jali_mesh_wrapper INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>)

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
    cinch_add_unit(test_jali_mesh_wrapper
      SOURCES test/test_jali_mesh_wrapper.cc
      LIBRARIES jali_mesh_wrapper
      # this is relative to the wrapper directory
      INPUTS ${CMAKE_SOURCE_DIR}/test_data/shotshell.exo
      ${CMAKE_SOURCE_DIR}/test_data/shotshell-v.exo
      POLICY MPI
      THREADS 1)
  endif ()

endif (ENABLE_Jali)
