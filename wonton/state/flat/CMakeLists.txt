#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

add_library(flat_state_wrapper INTERFACE)

target_include_directories(flat_state_wrapper INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  )

set(flat_state_HEADERS
  flat_state_wrapper.h
  flat_state_mm_wrapper.h
  )

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(flat_state_wrapper PROPERTIES
#   PUBLIC_HEADER "${flat_state_HEADERS}")
#
# Directly install files instead
install(FILES ${flat_state_HEADERS} DESTINATION include/wonton/state/flat)


target_link_libraries(flat_state_wrapper INTERFACE wonton_support)

target_link_libraries(wonton_state INTERFACE flat_state_wrapper)

install(TARGETS flat_state_wrapper
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/state/flat
  INCLUDES DESTINATION include/wonton/state/flat
  )
      
if (ENABLE_Jali)  
  if (ENABLE_UNIT_TESTS)

    # has to be processed after directory jali_state_wrapper is processed
    
    wonton_add_unit(test_flat_state_wrapper
      SOURCES test/test_flat_state_wrapper.cc
      LIBRARIES flat_state_wrapper jali_state_wrapper
      POLICY MPI
      THREADS 1)
  endif (ENABLE_UNIT_TESTS)
  
endif (ENABLE_Jali)
