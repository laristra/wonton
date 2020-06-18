#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

add_library(simple_state_wrapper INTERFACE)

target_include_directories(simple_state_wrapper INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

set(simple_state_HEADERS
  simple_state.h
  simple_state_wrapper.h
  simple_state_mm_wrapper.h
  )

# Not yet allowed for INTERFACE targets
#
# set_target_properties(simple_state_wrapper PROPERTIES
#   PUBLIC_HEADER "${simple_state_HEADERS}")
#
# Directly install files instead
install(FILES ${simple_state_HEADERS} DESTINATION include/wonton/state/simple)


target_link_libraries(simple_state_wrapper INTERFACE simple_mesh_wrapper)

target_link_libraries(simple_state_wrapper INTERFACE wonton_support)

target_link_libraries(wonton_state INTERFACE simple_state_wrapper)

install(TARGETS simple_state_wrapper
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/state/simple
  INCLUDES DESTINATION include/wonton/state/simple
  )

if (ENABLE_UNIT_TESTS)
  wonton_add_unit(test_simple_state
    SOURCES test/test_simple_state.cc
    LIBRARIES simple_state_wrapper
    POLICY SERIAL
    )
  
  wonton_add_unit(test_simple_state_wrapper
    SOURCES test/test_simple_state_wrapper.cc
    LIBRARIES simple_state_wrapper
    POLICY SERIAL)  
  
  wonton_add_unit(test_simple_state_mm_wrapper
    SOURCES test/test_simple_state_mm_wrapper.cc
    LIBRARIES simple_state_wrapper
    POLICY SERIAL)
endif (ENABLE_UNIT_TESTS)
  
  
