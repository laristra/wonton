#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

project(swarm)

add_library(wonton_swarm INTERFACE)

target_include_directories(wonton_swarm INTERFACE
    $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
    $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
    $<INSTALL_INTERFACE:include/wonton/swarm>)

# Add header files
set (wonton_swarm_HEADERS
    swarm.h
    swarm_state.h)

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(wonton_swarm PROPERTIES
#   PUBLIC_HEADER "${wonton_swarm_HEADERS}")
#
# Directly install files instead
#
install(FILES ${wonton_swarm_HEADERS} DESTINATION include/wonton/swarm)

target_link_libraries(wonton_swarm INTERFACE wonton_support)

target_link_libraries(wonton INTERFACE wonton_swarm)

install(TARGETS wonton_swarm
    EXPORT wonton_LIBRARIES
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/wonton/swarm
    INCLUDES DESTINATION include/wonton/swarm)


# Unit tests
if (ENABLE_UNIT_TESTS)
  wonton_add_unittest(test_swarm
      SOURCES test/test_swarm.cc
      LIBRARIES wonton_swarm
      POLICY SERIAL)

  wonton_add_unittest(test_swarm_state
      SOURCES test/test_swarm_state.cc
      LIBRARIES wonton_swarm
      POLICY SERIAL)

  wonton_add_unittest(test_swarm_state_multi
      SOURCES test/test_swarm_state_multi.cc
      LIBRARIES wonton_swarm
      POLICY SERIAL)
endif (ENABLE_UNIT_TESTS)
