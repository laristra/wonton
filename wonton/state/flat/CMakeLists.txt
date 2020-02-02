#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

set(flat_state_wrapper_HEADERS
  flat_state_mm_wrapper.h
  )

add_library(flat_state_wrapper INTERFACE)

target_include_directories(flat_state_wrapper INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(wonton_state INTERFACE flat_state_wrapper)

install(TARGETS flat_state_wrapper
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/state/flat
  INCLUDES DESTINATION include/wonton/state/flat
  )
      
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
