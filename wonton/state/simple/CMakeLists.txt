#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

add_library(simple_state_wrapper INTERFACE)

set(simple_state_wrapper_HEADERS
  simple_state_mm_wrapper.h
  )

target_sources(simple_state_wrapper state INTERFACE ${simple_state_wrapper_HEADERS})

target_link_libraries(wonton INTERFACE simple_state_wrapper)

if (ENABLE_UNIT_TESTS)
  cinch_add_unit(test_simple_state
    SOURCES test/test_simple_state.cc
    LIBRARIES simple_state_wrapper
    POLICY SERIAL
    )
  
  cinch_add_unit(test_simple_state_wrapper
    SOURCES test/test_simple_state_wrapper.cc
    LIBRARIES simple_state_wrapper
    POLICY SERIAL)  
  
  cinch_add_unit(test_simple_state_mm_wrapper
    SOURCES test/test_simple_state_mm_wrapper.cc
    LIBRARIES simple_state_wrapper
    POLICY SERIAL)
endif (ENABLE_UNIT_TESTS)
  
  
