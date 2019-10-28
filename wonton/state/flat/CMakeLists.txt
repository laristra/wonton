#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

add_library(flat_state_wrapper INTERFACE)

set(flat_state_wrapper_HEADERS
  flat_state_mm_wrapper.h
  )

target_link_libraries(state INTERFACE flat_state_wrapper)

if (Jali_DIR)  
  if (ENABLE_UNIT_TESTS)

    # has to be processed after directory jali_state_wrapper is processed
    
    cinch_add_unit(test_flat_state_wrapper
      SOURCES test/test_flat_state_wrapper.cc
      LIBRARIES flat_state_wrapper jali_state_wrapper
      POLICY MPI
      THREADS 1)
  endif (ENABLE_UNIT_TESTS)
  
endif (Jali_DIR)
