#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]


# Add header files

if (WONTON_ENABLE_Jali)
  add_library(jali_state_wrapper INTERFACE)
  
  target_include_directories(jali_state_wrapper INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
    )

  set(jali_state_wrapper_HEADERS
    jali_state_wrapper.h
    )

  # Not yet allowed for INTERFACE targets
  #
  #  set_target_properties(jali_state_wrapper PROPERTIES
  #    PUBLIC_HEADER "${jali_state_wrapper_HEADERS}")
  #
  # Directly install files instead
  install(FILES ${jali_state_wrapper_HEADERS}
    DESTINATION include/wonton/state/jali)

  
  target_link_libraries(jali_state_wrapper INTERFACE jali_mesh_wrapper)
  target_link_libraries(jali_state_wrapper INTERFACE wonton_support)

  target_link_libraries(jali_state_wrapper INTERFACE ${Jali_LIBRARIES})

  target_link_libraries(wonton_state INTERFACE jali_state_wrapper)

  install(TARGETS jali_state_wrapper
    EXPORT wonton_LIBRARIES
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/wonton/state/jali
    INCLUDES DESTINATION include/wonton/state/jali
    )
      
  if (ENABLE_UNIT_TESTS)
    cinch_add_unit(test_jali_state_wrapper
      SOURCES test/test_jali_state_wrapper.cc
      LIBRARIES jali_state_wrapper
      POLICY MPI
      THREADS 1)
  endif (ENABLE_UNIT_TESTS)
endif (ENABLE_Jali)
