#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]


# Add header files

if (FleCSI_DIR)
  add_library(flecsi_state_wrapper INTERFACE)
  
  target_include_directories(flecsi_state_wrapper PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
    $<INSTALL_INTERFACE:include)

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
endif (FleCSI_DIR)
