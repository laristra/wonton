#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]


# Add header files

if (Jali_DIR)
  add_library(jali_state_wrapper INTERFACE)
  
  set(jali_state_wrapper_HEADERS
    jali/jali_state_wrapper.h
    )

  target_sources(jali_state_wrapper INTERFACE ${jali_state_wrapper_HEADERS})
  target_link_libraries(jali_state_wrapper INTERFACE ${Jali_LIBRARIES})

  if (ENABLE_UNIT_TESTS)
    cinch_add_unit(test_jali_state_wrapper
      SOURCES test/test_jali_state_wrapper.cc
      LIBRARIES jali_state_wrapper)
  endif (ENABLE_UNIT_TESTS)
endif (Jali_DIR)
