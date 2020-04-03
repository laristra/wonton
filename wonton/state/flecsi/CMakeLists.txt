#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]


# Add header files

if (WONTON_ENABLE_FleCSI)
  add_library(flecsi_state_wrapper INTERFACE)
  
  target_include_directories(flecsi_state_wrapper INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
  )

  set(flecsi_state_wrapper_HEADERS
    flecsi_state_wrapper.h
    )

  # Not yet allowed for INTERFACE targets
  #
  #  set_target_properties(flecsi_state_wrapper PROPERTIES
  #    PUBLIC_HEADER "${flecsi_state_wrapper_HEADERS}")
  #
  # Directly install files instead
  install(FILES ${flecsi_state_wrapper_HEADERS}
    DESTINATION include/wonton/state/flecsi)
    
  target_link_libraries(flecsi_state_wrapper INTERFACE wonton_support)

  # Add flecsi dependencies
  target_link_libraries(wonton_state INTERFACE flecsi_state_wrapper)

  install(TARGETS flecsi_state_wrapper
    EXPORT wonton_LIBRARIES
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/wonton/state/flecsi
    INCLUDES DESTINATION include/wonton/state/flecsi
    )
  
  # Have to add FleCSI state wrapper tests
endif ()
