#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
project(state)

add_library(wonton_state INTERFACE)

target_include_directories(wonton_state INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/wonton/state>)

set(wonton_state_HEADERS
  state_manager.h
  state_vector_base.h
  state_vector_multi.h
  state_vector_multi_raw.h
  state_vector_uni.h
  state_vector_uni_raw.h
  )

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(wonton_state PROPERTIES
#   PUBLIC_HEADER "${wonton_state_HEADERS}")
#
# Directly install files instead
install(FILES ${wonton_state_HEADERS} DESTINATION include/wonton/state)

target_link_libraries(wonton_state INTERFACE wonton_support)

if (ENABLE_Jali)
  add_subdirectory(jali)
endif (ENABLE_Jali)

if (ENABLE_FleCSI)
  add_subdirectory(flecsi)
endif ()

add_subdirectory(simple)

add_subdirectory(flat)


target_link_libraries(wonton INTERFACE wonton_state)


install(TARGETS wonton_state
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/state
  INCLUDES DESTINATION include/wonton/state
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_state_vector_uni
    SOURCES test/test_state_vector_uni.cc
    LIBRARIES wonton_state
    POLICY SERIAL
    )

  cinch_add_unit(test_state_vector_multi
    SOURCES test/test_state_vector_multi.cc
    LIBRARIES wonton_state
    POLICY SERIAL
    )

  cinch_add_unit(test_state_manager
    SOURCES test/test_state_manager.cc
    LIBRARIES wonton_state
    POLICY SERIAL
    )

  cinch_add_unit(test_state_vector_uni_raw
    SOURCES test/test_state_vector_uni_raw.cc
    LIBRARIES wonton_state
    POLICY SERIAL
    )

  cinch_add_unit(test_state_vector_multi_raw
    SOURCES test/test_state_vector_multi_raw.cc
    LIBRARIES wonton_state
    POLICY SERIAL
    )

  cinch_add_unit(test_state_manager_raw
    SOURCES test/test_state_manager_raw.cc
    LIBRARIES wonton_state
    POLICY SERIAL
    )
  
  # HAVE TO ADD FLECSI STATE WRAPPER TESTS
  
endif (ENABLE_UNIT_TESTS)
